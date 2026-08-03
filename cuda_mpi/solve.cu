/*
solve.cu (CUDA + MPI)
*/

#include "ofd.h"
#include "ofd_cuda.h"
#include "ofd_prototype.h"
#include "finc.h"      // TPA 検証用の透過率測定で入射波形 finc() を使う (ofd.h の後)

#include "hdf5.h"
#include "ofd_hdf5.h"
#include <mpi.h>
#include <limits.h>
#define FILE_NAME "time_series_data.h5"

static void setup_cuda_mpi();
static void copy_to_host();

void solve(int io, double *tdft, FILE *fp)
{
    // HDF5ファイルの作成
    // 関数から?(fp の入替え?)
    // local

    double fmax[] = {0, 0};
    char   str[BUFSIZ];
    int    converged = 0;

    // setup boundary index (MPI)
    setup_mpi();

    // setup host memory
    setup_host();

    // setup (GPU)
    if (GPU) {
        setup_gpu();
        setup_cuda_mpi();
    }

    // initial field
    initfield();

    // TPA (二光子吸収) : material id -> β テーブル作成 (cuda/updateTpa.cu)
    if (NTpa) {
        setupTpa();
    }

    // TPA 検証用の透過率測定 : CW 波源 (waveamp) + 平面波 + point がある場合、
    // 最終 1 周期の point #1 位置の全電界振幅 |E_tot| から T = (|E_t|/E0)^2 を求める
    // (MPI では point を含むプロセスだけが測定し、最後に全プロセスで最大値を取る)
    const int tpaMon = (NTpa && IPlanewave && (WaveAmp > 0) && (NPoint > 0));
    double tpaEmax = 0;
    int tpaStart = INT_MAX;
    if (tpaMon) {
        const int nper = (int)(2 * PI / (WaveOmega * Dt)) + 1;  // 1 周期のステップ数
        tpaStart = Solver.maxiter - nper;
    }

    // HDF5 ファイルの作成 (rank 0 のみ / 直列アクセス)
    //
    // このファイルに書かれるのはすべて rank 0 が持つ値であり、グループ・
    // データセットの生成も書き込みも rank 0 だけが行っている (下記の
    // if (io) / if (commRank == 0) ブロック)。
    // ところが以前は H5Pset_fapl_mpio で MPI-IO ドライバを使って開いていた。
    // 並列 HDF5 では H5Gcreate / H5Dcreate / H5Dclose / H5Gclose / H5Fclose は
    // 全ランクが同じ順序で呼ぶ必要がある集団操作なので、rank 0 だけが呼ぶと
    // 他のランクと足並みが揃わずデッドロックする。
    // CPU+MPI 版 (mpi/solve.c) と同じ修正で、出力されるファイルの内容・構造は
    // 従来と同一。→ CUDA+MPI ビルドにも並列 HDF5 は不要。
    if (commRank == 0) {
        hdf5_open(1);
    }

    // time step iteration
    int itime;
    double t = 0;
    for (itime = 0; itime <= Solver.maxiter; itime++) {
        // update H
        t += 0.5 * Dt;
        updateHx(t);
        updateHy(t);
        updateHz(t);

        // ABC H
        if      (iABC == 0) {
            murH(numMurHx, (GPU ? d_fMurHx : fMurHx), Hx);
            murH(numMurHy, (GPU ? d_fMurHy : fMurHy), Hy);
            murH(numMurHz, (GPU ? d_fMurHz : fMurHz), Hz);
        }
        else if (iABC == 1) {
            pmlHx();
            pmlHy();
            pmlHz();
        }

        // PBC H
        if (PBCx) {
            if (Npx > 1) {
                comm_cuda_X(1);
            }
            else {
                pbcx();
            }
        }
        if (PBCy) {
            if (Npy > 1) {
                comm_cuda_Y(1);
            }
            else {
                pbcy();
            }
        }
        if (PBCz) {
            if (Npz > 1) {
                comm_cuda_Z(1);
            }
            else {
                pbcz();
            }
        }

        // share boundary H (MPI)
        if (Npx > 1) {
            comm_cuda_X(0);
        }
        if (Npy > 1) {
            comm_cuda_Y(0);
        }
        if (Npz > 1) {
            comm_cuda_Z(0);
        }

        // update E
        t += 0.5 * Dt;
        updateEx(t);
        updateEy(t);
        updateEz(t);

        // dispersion E
        if (numDispersionEx) {
            dispersionEx(t);
        }
        if (numDispersionEy) {
            dispersionEy(t);
        }
        if (numDispersionEz) {
            dispersionEz(t);
        }

        // ABC E
        if      (iABC == 1) {
            pmlEx();
            pmlEy();
            pmlEz();
        }

        // feed
        if (NFeed) {
            efeed(itime);
        }

        // inductor
        if (NInductor) {
            eload();
        }

        // TPA (二光子吸収) 非線形減衰 (cuda/updateTpa.cu)
        // 領域境界の E ハローを先に交換する : updateTpa は |E|^2 の colocated
        // 近似のため隣接セルの E 成分を読む (E 更新自体は H しか読まないので
        // 既存の comm_cuda_X/Y/Z(0) は H しか交換していない)
        if (NTpa) {
            if (Npx > 1) comm_cuda_E_X();
            if (Npy > 1) comm_cuda_E_Y();
            if (Npz > 1) comm_cuda_E_Z();
            updateTpa(t);
        }

        // TPA 検証用 : 最終 1 周期の全電界振幅を測定 (point を含むプロセスのみ)
        // (GPU 実行時は E が device 側にあるため、UM でない場合は測定できない)
        if (tpaMon && (itime >= tpaStart)) {
            if (GPU) cudaDeviceSynchronize();
            const int pi = Point[0].i;
            const int pj = Point[0].j;
            const int pk = Point[0].k;
            if (comm_inproc(pi, pj, pk)) {
                real_t fi = 0, dfi = 0;
                double e = 0;
                if      (Point[0].dir == 'X') {
                    finc(h_Xc[pi], h_Yn[pj], h_Zn[pk], t, Planewave.r0, Planewave.ri, Planewave.ei[0], Planewave.ai, Dt, &fi, &dfi);
                    e = EX(pi, pj, pk) + fi;
                }
                else if (Point[0].dir == 'Y') {
                    finc(h_Xn[pi], h_Yc[pj], h_Zn[pk], t, Planewave.r0, Planewave.ri, Planewave.ei[1], Planewave.ai, Dt, &fi, &dfi);
                    e = EY(pi, pj, pk) + fi;
                }
                else if (Point[0].dir == 'Z') {
                    finc(h_Xn[pi], h_Yn[pj], h_Zc[pk], t, Planewave.r0, Planewave.ri, Planewave.ei[2], Planewave.ai, Dt, &fi, &dfi);
                    e = EZ(pi, pj, pk) + fi;
                }
                tpaEmax = MAX(tpaEmax, fabs(e));
            }
        }

        // point
        if (NPoint) {
            vpoint(itime);
        }

        // DFT
        if (GPU) cudaDeviceSynchronize();
        const double t0 = comm_cputime();
        dftNear3d(itime);
        if (GPU) cudaDeviceSynchronize();
        *tdft += comm_cputime() - t0;

        // average and convergence
        if ((itime % Solver.nout == 0) || (itime == Solver.maxiter)) {
            // average
            double fsum[2];
            average(fsum);

            // allreduce average (MPI)
            if (commSize > 1) {
                comm_average(fsum);
            }

            // average (post)
            if (commRank == 0) {
                Eiter[Niter] = fsum[0];
                Hiter[Niter] = fsum[1];
                Niter++;
            }

            // monitor
            if (io) {
                sprintf(str, "%7d %.6f %.6f", itime, fsum[0], fsum[1]);
                fprintf(fp,     "%s\n", str);
                fprintf(stdout, "%s\n", str);
                fflush(fp);
                fflush(stdout);
                
                // copy near3d from device to host
                memcopy3_gpu();

            }

            // HDF5 : 瞬時値スナップショット
            // comm_snapshot() は全ランクが参加する集団操作なので、
            // io (= rank 0 か) の条件の外で呼ぶこと (mpi/solve.c と同じ)。
            if (Hdf5Output) {
                if (GPU) cudaDeviceSynchronize();
                comm_snapshot();
                if (commRank == 0) {
                    setupSize(1, 1, 1, 0);
                    hdf5_write_snapshot(itime, t, g_Ex, g_Ey, g_Ez, g_Hx, g_Hy, g_Hz);
                    setupSize(Npx, Npy, Npz, commRank);
                }
            }

            // check convergence
            fmax[0] = MAX(fmax[0], fsum[0]);
            fmax[1] = MAX(fmax[1], fsum[1]);
            if ((fsum[0] < fmax[0] * Solver.converg) &&
                (fsum[1] < fmax[1] * Solver.converg)) {
                converged = 1;
                break;
            }
        }
    }

    // result
    if (io) {
        sprintf(str, "    --- %s ---", (converged ? "converged" : "max steps"));
        fprintf(fp,     "%s\n", str);
        fprintf(stdout, "%s\n", str);
        fflush(fp);
        fflush(stdout);
    }

    // TPA 検証用 : 透過率を出力 (CI が ofd.log のこの行を判定に使う)
    // point を含むプロセスだけが測定しているので、全プロセスで最大値を取って集約する
    // (非所有プロセスの tpaEmax は 0 なので MAX で正しい値が得られる)
    if (tpaMon) {
        double gmax = 0;
        MPI_Allreduce(&tpaEmax, &gmax, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        tpaEmax = gmax;
    }
    if (io && tpaMon) {
        // I0 = (1/2) ε0 c E0^2 : 入射平面波 (真空中) の強度
        const double i0 = 0.5 * EPS0 * C * WaveAmp * WaveAmp;
        const double trans = (tpaEmax / WaveAmp) * (tpaEmax / WaveAmp);
        sprintf(str, "TPA: transmission = %.6f (I0=%.6e W/m^2)", trans, i0);
        fprintf(fp,     "%s\n", str);
        fprintf(stdout, "%s\n", str);
        fflush(fp);
        fflush(stdout);
    }

    // time steps
    Ntime = itime + converged;

    // copy point from device to host
    if (GPU) {
        copy_to_host();
    }

    // グループの作成前に同期
    //MPI_Barrier(MPI_COMM_WORLD);

    // HDF5 : 周波数領域の結果は comm_near3d() の集約後に書く
    // (この関数の末尾を参照)

    // free
    memfree2_gpu();

    // copy near3d from device to host
    memcopy3_gpu();

    // free
    memfree3_gpu();

    // MPI : send to root
    if (commSize > 1) {
        // feed waveform
        if (NFeed) {
            comm_feed();
        }

        // point waveform
        if (NPoint) {
            comm_point();
        }

        // near3d
        if (NFreq2) {
            comm_near3d();
        }
    }

    // HDF5 : 周波数領域の最終結果・収束履歴・メタデータ
    //
    // 1 プロセスのときは上の集約ブロックごとスキップされるので、この書き出しを
    // 中に置くと 1 プロセスだけ /freqdomain が出ない。ブロックの外に置き、
    // 集約の有無で参照する配列を切り替える (mpi/solve.c と同じ)。
    if (commRank == 0) {
        const int gathered = (commSize > 1) && NFreq2;
        hdf5_write_freqdomain(
            gathered ? g_cEx_r : cEx_r, gathered ? g_cEx_i : cEx_i,
            gathered ? g_cEy_r : cEy_r, gathered ? g_cEy_i : cEy_i,
            gathered ? g_cEz_r : cEz_r, gathered ? g_cEz_i : cEz_i,
            gathered ? g_cHx_r : cHx_r, gathered ? g_cHx_i : cHx_i,
            gathered ? g_cHy_r : cHy_r, gathered ? g_cHy_i : cHy_i,
            gathered ? g_cHz_r : cHz_r, gathered ? g_cHz_i : cHz_i);
        hdf5_write_convergence(Niter, Eiter, Hiter);
        hdf5_close();
    }
}


// setup
static void setup_cuda_mpi()
{
    size_t size;
    //printf("%d %d %d %d %d %d %d %d\n", commSize, commRank, bid.numhy_x, bid.numhz_x, bid.numhz_y, bid.numhx_y, bid.numhx_z, bid.numhy_z); fflush(stdout);

    // X boundary
    size = Bid.numhy_x * sizeof(real_t);
    cuda_malloc(GPU, UM, (void **)&d_Sendbuf_hy_x, size);
    cuda_malloc(GPU, UM, (void **)&d_Recvbuf_hy_x, size);

    size = Bid.numhz_x * sizeof(real_t);
    cuda_malloc(GPU, UM, (void **)&d_Sendbuf_hz_x, size);
    cuda_malloc(GPU, UM, (void **)&d_Recvbuf_hz_x, size);

    // Y boundary
    size = Bid.numhz_y * sizeof(real_t);
    cuda_malloc(GPU, UM, (void **)&d_Sendbuf_hz_y, size);
    cuda_malloc(GPU, UM, (void **)&d_Recvbuf_hz_y, size);

    size = Bid.numhx_y * sizeof(real_t);
    cuda_malloc(GPU, UM, (void **)&d_Sendbuf_hx_y, size);
    cuda_malloc(GPU, UM, (void **)&d_Recvbuf_hx_y, size);

    // Z boundary
    size = Bid.numhx_z * sizeof(real_t);
    cuda_malloc(GPU, UM, (void **)&d_Sendbuf_hx_z, size);
    cuda_malloc(GPU, UM, (void **)&d_Recvbuf_hx_z, size);

    size = Bid.numhy_z * sizeof(real_t);
    cuda_malloc(GPU, UM, (void **)&d_Sendbuf_hy_z, size);
    cuda_malloc(GPU, UM, (void **)&d_Recvbuf_hy_z, size);

    // block
    bufBlock = dim3(16, 16);

    // parameter
    //cudaMemcpyToSymbol(d_bid, &bid, sizeof(bid_t));
}


// copy from device to host
static void copy_to_host()
{
    if (NFeed) {
        cuda_memcpy(GPU, VFeed, d_VFeed, Feed_size, cudaMemcpyDeviceToHost);
        cuda_memcpy(GPU, IFeed, d_IFeed, Feed_size, cudaMemcpyDeviceToHost);
    }

    if (NPoint) {
        cuda_memcpy(GPU, VPoint, d_VPoint, Point_size, cudaMemcpyDeviceToHost);
    }
}
