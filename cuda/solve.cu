/*
solve.cu (CUDA)
*/

#include "ofd.h"
#include "ofd_cuda.h"
#include "ofd_prototype.h"
#include "finc.h"      // TPA 検証用の透過率測定 (ofd.h の後)

#include "hdf5.h"
#include "ofd_hdf5.h"
#define FILE_NAME "time_series_data.h5"

static void copy_to_host();

// ── HDM (GPU かつ UM でない) 実行での HDF5 スナップショット ──────────────
// Ex..Hz は cudaMalloc した device メモリなので host から直接読めない
// (読むと 0xC0000005 で落ちる — RTX 3060 実機で確認。UM / -cpu では出ない)。
// スナップショットを書くステップだけ host の作業バッファへコピーし、
// そちらを出力側に渡す。UM / CPU 実行では場の配列をそのまま返す。
static real_t *h_snap[6] = { NULL, NULL, NULL, NULL, NULL, NULL };
static int     h_snap_failed = 0;

// 戻り値: 1 = f[0..5] に host から読める Ex,Ey,Ez,Hx,Hy,Hz を入れた /
//         0 = 作業バッファを確保できない (出力を諦める)
static int snapshot_host_fields(const real_t *f[6])
{
	real_t *dev[6];
	dev[0] = Ex; dev[1] = Ey; dev[2] = Ez;
	dev[3] = Hx; dev[4] = Hy; dev[5] = Hz;
	if (!GPU || UM) {
		for (int n = 0; n < 6; n++) f[n] = dev[n];
		return 1;
	}
	if (h_snap_failed) return 0;
	const size_t size = (size_t)NN * sizeof(real_t);
	for (int n = 0; n < 6; n++) {
		if (h_snap[n] == NULL) h_snap[n] = (real_t *)malloc(size);
		if (h_snap[n] == NULL) {
			fprintf(stderr, "*** snapshot host buffer malloc error (%zu bytes) : HDF5 /timeseries is skipped\n", size);
			fflush(stderr);
			h_snap_failed = 1;
			return 0;
		}
		cudaMemcpy(h_snap[n], dev[n], size, cudaMemcpyDeviceToHost);
		f[n] = h_snap[n];
	}
	return 1;
}

static void snapshot_host_free(void)
{
	for (int n = 0; n < 6; n++) {
		free(h_snap[n]);
		h_snap[n] = NULL;
	}
}

void solve(int io, double *tdft, FILE *fp)
{
    // HDF5ファイルの作成
    // 関数から?(fp の入替え?)
    // local

    double fmax[] = {0, 0};
    char   str[BUFSIZ];
    int    converged = 0;

    // setup host memory
    setup_host();

    // setup (GPU)
    if (GPU) {
        setup_gpu();
    }

    // initial field
    initfield();

    // TPA (二光子吸収) : material id -> β テーブル作成 (cuda/updateTpa.cu)
    if (NTpa) {
        setupTpa();
    }

    // TPA 検証用の透過率測定 (CPU 版 sol/solve.c と同じ)
    const int tpaMon = (NTpa && IPlanewave && (WaveAmp > 0) && (NPoint > 0));
    double tpaEmax = 0;
    int tpaStart = INT_MAX;
    if (tpaMon) {
        const int nper = (int)(2 * PI / (WaveOmega * Dt)) + 1;  // 1 周期のステップ数
        tpaStart = Solver.maxiter - nper;
    }

    // HDF5ファイルの作成
    hdf5_open(1);

    // time step iteration
    int itime;
    double t = 0;

    int lx = (iABC == 0) ? 1 : (iABC == 1) ? cPML.l : 0;
    int ly = (iABC == 0) ? 1 : (iABC == 1) ? cPML.l : 0;
    int lz = (iABC == 0) ? 1 : (iABC == 1) ? cPML.l : 0;

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
            pbcx();
        }
        if (PBCy) {
            pbcy();
        }
        if (PBCz) {
            pbcz();
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
        if (NTpa) {
            updateTpa(t);
        }

        // TPA 検証用 : 最終 1 周期の全電界振幅を測定
        // (GPU 実行時は E が device 側にあるため、UM でない場合は測定できない)
        if (tpaMon && (itime >= tpaStart)) {
            if (GPU) cudaDeviceSynchronize();
            const int pi = Point[0].i;
            const int pj = Point[0].j;
            const int pk = Point[0].k;
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

        // point
        if (NPoint) {
            vpoint(itime);
        }

        // DFT
        if (GPU) cudaDeviceSynchronize();
        const double t0 = cputime();
        dftNear3d(itime);
        if (GPU) cudaDeviceSynchronize();
        *tdft += cputime() - t0;

        // average and convergence
        if ((itime % Solver.nout == 0) || (itime == Solver.maxiter)) {
            // average
            double fsum[2];
            average(fsum);

            // average (plot)
            // Niter は格納直後に加算する。以前は収束判定の後 (下の
            // 「Niterを増加」) で加算していたため、収束して break した
            // 最後の 1 点が収束履歴から落ちていた (CPU / MPI / CUDA+MPI は
            // どれも格納直後に加算しており、この実装だけがずれていた。
            // 実測: dipole の iter が CPU の 0..550 に対し CUDA は 0..500、
            // HDF5 の metadata/Niter も 12 に対し 11)。
            Eiter[Niter] = fsum[0];
            Hiter[Niter] = fsum[1];
            Niter++;

            // monitor
            if (io) {
                sprintf(str, "%7d %.6f %.6f", itime, fsum[0], fsum[1]);
                fprintf(fp,     "%s\n", str);
                fprintf(stdout, "%s\n", str);
                fflush(fp);
                fflush(stdout);

                // copy near3d from device to host
                memcopy3_gpu();

                // HDF5 : 瞬時値スナップショット (sol/outputHdf5.c)。
                // HDM では device → host コピーを挟む (snapshot_host_fields)。
                if (hdf5_snapshot_enabled(itime)) {
                    if (GPU) cudaDeviceSynchronize();
                    const real_t *f[6];
                    if (snapshot_host_fields(f)) {
                        hdf5_write_snapshot(itime, t, f[0], f[1], f[2], f[3], f[4], f[5]);
                    }
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
    if (io && tpaMon) {
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
        // DFT 配列 (cEx_r 等) は UM なので、ホストから読む前に同期する
        cudaDeviceSynchronize();
    }

    // HDF5 : 周波数領域の最終結果・収束履歴・メタデータ (sol/outputHdf5.c)
    hdf5_write_freqdomain(
        cEx_r, cEx_i, cEy_r, cEy_i, cEz_r, cEz_i,
        cHx_r, cHx_i, cHy_r, cHy_i, cHz_r, cHz_i);
    hdf5_write_convergence(Niter, Eiter, Hiter);
    hdf5_close();

    // free
    snapshot_host_free();
    memfree2_gpu();

    // copy near3d from device to host
    memcopy3_gpu();

    // free
    memfree3_gpu();
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

