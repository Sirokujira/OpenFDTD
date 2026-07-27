#include "ofd.h"
#include "complex.h"
#include "ofd_prototype.h"
#include "ev.h"
#include "finc.h"

#include "hdf5.h"

#define FILE_NAME "time_series_data.h5"

// 温度更新関数
//void updateTemperature(double *T, int Nx, int Ny, int Nz, double alpha, double Dt, double *P_loss) {
// 温度更新関数
void updateTemperature(double *T, int64_t NN, int NFreq2, double alpha, double Dt, double *P_loss, double Dx, double Dy, double Dz) {
    for (int ifreq = 0; ifreq < NFreq2; ifreq++) {
        const int64_t base_idx = (int64_t)ifreq * NN;
        // NA マクロで 3 次元格子の隣接セルを参照する (内点のみ更新、境界は固定境界)
        for (int i = 1; i < Nx - 1; i++) {
        for (int j = 1; j < Ny - 1; j++) {
        for (int k = 1; k < Nz - 1; k++) {
            const int64_t idx = base_idx + NA(i, j, k);
            const double d2Tdx2 = (T[base_idx + NA(i + 1, j, k)] - 2.0 * T[idx] + T[base_idx + NA(i - 1, j, k)]) / (Dx * Dx);
            const double d2Tdy2 = (T[base_idx + NA(i, j + 1, k)] - 2.0 * T[idx] + T[base_idx + NA(i, j - 1, k)]) / (Dy * Dy);
            const double d2Tdz2 = (T[base_idx + NA(i, j, k + 1)] - 2.0 * T[idx] + T[base_idx + NA(i, j, k - 1)]) / (Dz * Dz);

            // 温度の更新
            T[idx] += alpha * Dt * (d2Tdx2 + d2Tdy2 + d2Tdz2) + Dt * P_loss[idx];
        }
        }
        }
    }
}


// 発熱量の計算 (セル毎に材料を参照する)
//
// 各セル・各成分について
//   p = (1/2) σ_e |E|^2 + (1/2) σ_m |H|^2   [W/m^3]
// を足し合わせる。σ_e は導電率 [S/m] (material の esgm)、σ_m は磁気
// 導電率 [Ω/m] (同 msgm) で、FDTD の損失項そのものなので μ'' を別に
// 与える必要はない (σ_m = ω μ'' の関係にある)。
//
// Yee 格子では 6 成分がそれぞれ別の位置にあり、材料 ID も成分ごとに
// 独立の配列 (iEx..iHz) を持つ。以前はセル位置に依らず材料 0 (通常は
// 真空 = σ_e 0) の σ を全セルに使い、磁気損失は定数 1e-3 を仮置きして
// いたため、発熱密度は実質的に構造を反映していなかった。
//
// 注意: DFT 配列 (cEx_r 等) は float。double* で受けると読み越しになり
// Windows ではアクセス違反になる (glibc では偶然動作していた)。
// 材料 ID 配列 (iEx 等) は id_t で、ビルド構成により型幅が変わる。
static void calculatePowerLoss(double *P_loss, int64_t nn, int nfreq2,
                        const double *esgm, const double *msgm,
                        const id_t *iex, const id_t *iey, const id_t *iez,
                        const id_t *ihx, const id_t *ihy, const id_t *ihz,
                        const float *cEx_r, const float *cEx_i, const float *cEy_r, const float *cEy_i,
                        const float *cEz_r, const float *cEz_i,
                        const float *cHx_r, const float *cHx_i, const float *cHy_r, const float *cHy_i,
                        const float *cHz_r, const float *cHz_i) {
    for (int ifreq = 0; ifreq < nfreq2; ifreq++) {
        int64_t base_idx = (int64_t)ifreq * nn;
        int64_t i;
#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (i = 0; i < nn; i++) {
            const int64_t n = base_idx + i;

            // 電界による損失 : 成分ごとに自分の位置の材料の σ_e を使う
            const double ex2 = ((double)cEx_r[n] * cEx_r[n]) + ((double)cEx_i[n] * cEx_i[n]);
            const double ey2 = ((double)cEy_r[n] * cEy_r[n]) + ((double)cEy_i[n] * cEy_i[n]);
            const double ez2 = ((double)cEz_r[n] * cEz_r[n]) + ((double)cEz_i[n] * cEz_i[n]);
            const double pe = (esgm[iex[i]] * ex2)
                            + (esgm[iey[i]] * ey2)
                            + (esgm[iez[i]] * ez2);

            // 磁界による損失 : 同様に σ_m を使う
            const double hx2 = ((double)cHx_r[n] * cHx_r[n]) + ((double)cHx_i[n] * cHx_i[n]);
            const double hy2 = ((double)cHy_r[n] * cHy_r[n]) + ((double)cHy_i[n] * cHy_i[n]);
            const double hz2 = ((double)cHz_r[n] * cHz_r[n]) + ((double)cHz_i[n] * cHz_i[n]);
            const double pm = (msgm[ihx[i]] * hx2)
                            + (msgm[ihy[i]] * hy2)
                            + (msgm[ihz[i]] * hz2);

            P_loss[n] = 0.5 * (pe + pm);
        }
    }
}

// 材料 id -> 損失パラメータの表を作る (発熱密度の計算用)
//   esgm[m] : 導電率     σ_e [S/m]
//   msgm[m] : 磁気導電率 σ_m [Ω/m]
// 材料 ID 配列 (iEx 等) が指す先を毎セル Material[] から引くと構造体
// アクセスが入るので、平坦な double 配列にしておく。
// PEC (id = 1) は E も H も更新されず場が 0 なので、σ を入れても
// 発熱には寄与しない (Material[1] の値をそのまま使ってよい)。
static void setup_loss_table(double **esgm, double **msgm) {
    *esgm = (double *)malloc(NMaterial * sizeof(double));
    *msgm = (double *)malloc(NMaterial * sizeof(double));
    if ((*esgm == NULL) || (*msgm == NULL)) {
        // 既存の malloc エラー表示 (NN) と同じく size_t にキャストして %zu で出す
        fprintf(stderr, "*** loss table malloc error (NMaterial=%zu)\n", (size_t)NMaterial);
        exit(1);
    }
    for (int64_t m = 0; m < NMaterial; m++) {
        (*esgm)[m] = Material[m].esgm;
        (*msgm)[m] = Material[m].msgm;
    }
}

void solve(int io, double *tdft, FILE *fp) {
    // HDF5ファイルの作成
        // 関数から?(fp の入替え?)
    hid_t file_id;
        // local
        hid_t group_id, dataset_id, dataspace_id, memspace_id = -1;
    herr_t status;

    double fmax[] = {0, 0};
    char str[BUFSIZ];
    int converged = 0;

    // initial field
    initfield();

    // TPA (二光子吸収) : material id -> β テーブル作成 (sol/updateTpa.c)
    if (NTpa) {
        setupTpa();
    }

    // TPA 検証用の透過率測定 : CW 波源 (waveamp) + 平面波 + point がある場合、
    // 最終 1 周期の point #1 位置の全電界振幅 |E_tot| から T = (|E_t|/E0)^2 を求める
    // (point #1 の方向成分が入射偏波と一致していること。data/sample/tpa_slab.ofd 参照)
    const int tpaMon = (NTpa && IPlanewave && (WaveAmp > 0) && (NPoint > 0));
    double tpaEmax = 0;
    int tpaStart = INT_MAX;
    if (tpaMon) {
        const int nper = (int)(2 * PI / (WaveOmega * Dt)) + 1;  // 1 周期のステップ数
        tpaStart = Solver.maxiter - nper;
    }

    // 温度配列の初期化
    //int Nx = 100, Ny = 100, Nz = 100;
    double alpha = 0.01;  // 熱拡散係数
    //double *T = (double *)malloc(Nx * Ny * Nz * sizeof(double));
    //double *P_loss = (double *)malloc(Nx * Ny * Nz * sizeof(double));
    //memset(T, 0, Nx * Ny * Nz * sizeof(double));
    //memset(P_loss, 0, Nx * Ny * Nz * sizeof(double));
    //NN
    double *T = (double *)malloc(NFreq2 * NN * sizeof(double));
    double *P_losses = (double *)malloc(NFreq2 * NN * sizeof(double));
    if ((T == NULL) || (P_losses == NULL)) {
        fprintf(stderr, "*** temperature array malloc error (NFreq2=%d NN=%zu)\n", NFreq2, (size_t)NN);
        exit(1);
    }
    memset(T, 0, NFreq2 * NN * sizeof(double));
    memset(P_losses, 0, NFreq2 * NN * sizeof(double));

    // セルの幅（空間ステップ）を計算
    double Dx = (Xn[Nx] - Xn[0]) / Nx;
    double Dy = (Yn[Ny] - Yn[0]) / Ny;
    double Dz = (Zn[Nz] - Zn[0]) / Nz;
    sprintf(str, "%.6f %.6f %.6f", Dx, Dy, Dz);
    fprintf(stdout, "%s\n", str);

    // 初期温度設定
    for (int i = 0; i < NFreq2 * NN; i++) {
        T[i] = 20.0;  // 初期温度（20度）
    }

    // HDF5ファイルの作成
    file_id = H5Fcreate(FILE_NAME, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    // 発熱密度に使う材料毎の損失パラメータ (セル毎に材料 ID 配列から引く)
    double *loss_esgm = NULL, *loss_msgm = NULL;
    setup_loss_table(&loss_esgm, &loss_msgm);

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
            murH(numMurHx, fMurHx, Hx);
            murH(numMurHy, fMurHy, Hy);
            murH(numMurHz, fMurHz, Hz);
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

        // TPA (二光子吸収) 非線形減衰 (sol/updateTpa.c)
        if (NTpa) {
            updateTpa(t);
        }

        // point
        if (NPoint) {
            vpoint(itime);
        }

        // TPA 検証用 : 最終 1 周期の全電界振幅を測定
        if (tpaMon && (itime >= tpaStart)) {
            const int pi = Point[0].i;
            const int pj = Point[0].j;
            const int pk = Point[0].k;
            real_t fi = 0, dfi = 0;
            double e = 0;
            if      (Point[0].dir == 'X') {
                finc(Xc[pi], Yn[pj], Zn[pk], t, Planewave.r0, Planewave.ri, Planewave.ei[0], Planewave.ai, Dt, &fi, &dfi);
                e = EX(pi, pj, pk) + fi;
            }
            else if (Point[0].dir == 'Y') {
                finc(Xn[pi], Yc[pj], Zn[pk], t, Planewave.r0, Planewave.ri, Planewave.ei[1], Planewave.ai, Dt, &fi, &dfi);
                e = EY(pi, pj, pk) + fi;
            }
            else if (Point[0].dir == 'Z') {
                finc(Xn[pi], Yn[pj], Zc[pk], t, Planewave.r0, Planewave.ri, Planewave.ei[2], Planewave.ai, Dt, &fi, &dfi);
                e = EZ(pi, pj, pk) + fi;
            }
            tpaEmax = MAX(tpaEmax, fabs(e));
        }

        // DFT
        const double t0 = cputime();
        dftNear3d(itime);
        *tdft += cputime() - t0;

        // 発熱量の計算 (セル毎に材料を参照する)
        calculatePowerLoss(P_losses, NN, NFreq2, loss_esgm, loss_msgm,
                           iEx, iEy, iEz, iHx, iHy, iHz,
                           cEx_r, cEx_i, cEy_r, cEy_i, cEz_r, cEz_i,
                           cHx_r, cHx_i, cHy_r, cHy_i, cHz_r, cHz_i);

        // 温度の更新
        //updateTemperature(T, Nx, Ny, Nz, alpha, Dt, P_loss);
        // 温度の更新
        updateTemperature(T, NN, NFreq2, alpha, Dt, P_losses, Dx, Dy, Dz);

        // average and convergence
        if ((itime % Solver.nout == 0) || (itime == Solver.maxiter)) {
            // average
            double fsum[2];
            average(fsum);

            // average (post)
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

                // 各時間ステップごとにグループを作成
                char group_name[32];
                snprintf(group_name, sizeof(group_name), "/data%06d", itime);
                group_id = H5Gcreate(file_id, group_name, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

                // Eフィールドデータセットの作成と書き込み
                hsize_t e_dims[4] = {1, NFreq2, NN, 6};
                dataspace_id = H5Screate_simple(4, e_dims, NULL);
                dataset_id = H5Dcreate(group_id, "E", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

                // 書き込み用のメモリスペースを修正
                hsize_t mem_dims[1] = {6};
                memspace_id = H5Screate_simple(1, mem_dims, NULL);

                for (int ifreq = 0; ifreq < NFreq2; ifreq++) {
                    int64_t n0 = ifreq * NN;
                    for (int nn = 0; nn < NN; nn++) {
                        double e_value[6] = {
                            cEx_r[n0 + nn], cEy_r[n0 + nn], cEz_r[n0 + nn],
                            cEx_i[n0 + nn], cEy_i[n0 + nn], cEz_i[n0 + nn]
                        };

                        hsize_t e_offset[4] = {0, ifreq, nn, 0};
                        hsize_t e_count[4] = {1, 1, 1, 6};
                        H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, e_offset, NULL, e_count, NULL);

                        // 書き込み
                        status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, memspace_id, dataspace_id, H5P_DEFAULT, e_value);
                        if (status < 0) {
                            fprintf(stderr, "Error writing E data at itime=%d, ifreq=%d, nn=%d\n", itime, ifreq, nn);
                        }
                    }
                }
                H5Dclose(dataset_id);
                H5Sclose(dataspace_id);

                // Hフィールドデータセットの作成と書き込み
                hsize_t h_dims[4] = {1, NFreq2, NN, 6};
                dataspace_id = H5Screate_simple(4, h_dims, NULL);
                dataset_id = H5Dcreate(group_id, "H", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

                for (int ifreq = 0; ifreq < NFreq2; ifreq++) {
                    int64_t n0 = ifreq * NN;
                    for (int nn = 0; nn < NN; nn++) {
                        double h_value[6] = {
                            cHx_r[n0 + nn], cHy_r[n0 + nn], cHz_r[n0 + nn],
                            cHx_i[n0 + nn], cHy_i[n0 + nn], cHz_i[n0 + nn]
                        };

                        hsize_t h_offset[4] = {0, ifreq, nn, 0};
                        hsize_t h_count[4] = {1, 1, 1, 6};
                        H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, h_offset, NULL, h_count, NULL);
                        status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, memspace_id, dataspace_id, H5P_DEFAULT, h_value);
                        if (status < 0) {
                            fprintf(stderr, "Error writing H data at itime=%d, ifreq=%d, nn=%d\n", itime, ifreq, nn);
                        }
                    }
               }
                H5Dclose(dataset_id);
                H5Sclose(dataspace_id);

                // 複素数用のHDF5データ型を定義
                hid_t complex_datatype = H5Tcreate(H5T_COMPOUND, sizeof(d_complex_t));
                H5Tinsert(complex_datatype, "real", HOFFSET(d_complex_t, r), H5T_NATIVE_DOUBLE);
                H5Tinsert(complex_datatype, "imag", HOFFSET(d_complex_t, i), H5T_NATIVE_DOUBLE);

                // Surfaceフィールドデータセットの作成と書き込み
                hsize_t surf_dims[4] = {1, NFreq2, NN, 6};
                dataspace_id = H5Screate_simple(4, surf_dims, NULL);
                dataset_id = H5Dcreate(group_id, "Surface", complex_datatype, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

                for (int ifreq = 0; ifreq < NFreq2; ifreq++) {
                    //int64_t surf0 = ifreq * NSurface;
                    for (int surf = 0; surf < NSurface; surf++) {
                        d_complex_t surf_value[6] = {
                            SurfaceEx[ifreq][surf], SurfaceEy[ifreq][surf], SurfaceEz[ifreq][surf],
                            SurfaceHx[ifreq][surf], SurfaceHy[ifreq][surf], SurfaceHz[ifreq][surf]
                        };

                        hsize_t surf_offset[4] = {0, ifreq, surf, 0};
                        hsize_t surf_count[4] = {1, 1, 1, 6};
                        H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, surf_offset, NULL, surf_count, NULL);
                        status = H5Dwrite(dataset_id, complex_datatype, memspace_id, dataspace_id, H5P_DEFAULT, surf_value);
                        if (status < 0) {
                            fprintf(stderr, "Error writing H data at itime=%d, ifreq=%d, surf=%d\n", itime, ifreq, surf);
                        }
                    }
                }
                H5Dclose(dataset_id);
                H5Sclose(dataspace_id);

                // Pフィールドデータセットの作成と書き込み（仮の例）
                //hsize_t p_dims[4] = {1, NFreq2, NN, 3};
            	hsize_t p_dims[3] = {NFreq2, NN, 3};
                dataspace_id = H5Screate_simple(3, p_dims, NULL);
                dataset_id = H5Dcreate(group_id, "P", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

                // 書き込み用のメモリスペースを修正
                hsize_t mem_dims2[1] = {3};
                if (memspace_id >= 0) H5Sclose(memspace_id);
                memspace_id = H5Screate_simple(1, mem_dims2, NULL);

                for (int ifreq = 0; ifreq < NFreq2; ifreq++) {
                    int64_t n0 = ifreq * NN;
                    for (int nn = 0; nn < NN; nn++) {
                        double p_value[3] = {
                            cEx_r[n0 + nn] * cHy_r[n0 + nn] - cEy_r[n0 + nn] * cHx_r[n0 + nn],
                            cEy_r[n0 + nn] * cHz_r[n0 + nn] - cEz_r[n0 + nn] * cHy_r[n0 + nn],
                            cEz_r[n0 + nn] * cHx_r[n0 + nn] - cEx_r[n0 + nn] * cHz_r[n0 + nn]
                        };

                        hsize_t p_offset[3] = {ifreq, nn, 0};
                        hsize_t p_count[3] = {1, 1, 3};
                        H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, p_offset, NULL, p_count, NULL);
                        status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, memspace_id, dataspace_id, H5P_DEFAULT, p_value);
                        if (status < 0) {
                            fprintf(stderr, "Error writing P data at itime=%d, ifreq=%d, nn=%d\n", itime, ifreq, nn);
                        }
                    }
                }
                H5Dclose(dataset_id);
                H5Sclose(dataspace_id);

                // 発熱量の計算
                hsize_t p_dims2[4] = {1, NFreq2, NN, 1};
                dataspace_id = H5Screate_simple(4, p_dims2, NULL);
                dataset_id = H5Dcreate(group_id, "P_loss", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

                // 書き込み用のメモリスペースを修正
                hsize_t mem_dims3[1] = {1};
                if (memspace_id >= 0) H5Sclose(memspace_id);
                memspace_id = H5Screate_simple(1, mem_dims3, NULL);

                // 材料の特性設定
                for (int ifreq = 0; ifreq < NFreq2; ifreq++) {
                    int64_t n0 = ifreq * NN;
                    for (int nn = 0; nn < NN; nn++) {
                        // 発熱量密度の計算
                        double P_loss[1] = { P_losses[n0 + nn] };

                        hsize_t p_offset[4] = {0, ifreq, nn, 0};
                        hsize_t p_count[4] = {1, 1, 1, 1};
                        H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, p_offset, NULL, p_count, NULL);
                        status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, memspace_id, dataspace_id, H5P_DEFAULT, P_loss);
                        if (status < 0) {
                            fprintf(stderr, "Error writing P data at itime=%d, ifreq=%d, nn=%d\n", itime, ifreq, nn);
                        }
                    }
                }
                H5Dclose(dataset_id);
                H5Sclose(dataspace_id);

                // グループのクローズ
                H5Gclose(group_id);
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
    // メモリスペース、データセットとデータスペースのクローズ
    if (memspace_id >= 0) status = H5Sclose(memspace_id);

    // result
    if (io) {
        sprintf(str, "    --- %s ---", (converged ? "converged" : "max steps"));
        fprintf(fp,     "%s\n", str);
        fprintf(stdout, "%s\n", str);
        fflush(fp);
        fflush(stdout);
    }

    // TPA 検証用 : 透過率を出力 (CI が ofd.log のこの行を判定に使う)
    // I0 = (1/2) ε0 c E0^2 : 入射平面波 (真空中) の強度
    if (io && tpaMon) {
        const double i0 = 0.5 * EPS0 * C * WaveAmp * WaveAmp;
        const double trans = (tpaEmax / WaveAmp) * (tpaEmax / WaveAmp);
        sprintf(str, "TPA: transmission = %.6f (I0=%.6e W/m^2)", trans, i0);
        fprintf(fp,     "%s\n", str);
        fprintf(stdout, "%s\n", str);
        fflush(fp);
        fflush(stdout);
    }

    // 熱解析レイヤの診断 : 発熱密度をセル体積で重み付けして全体を積算する。
    // DFT (cEx_r 等) は入射スペクトルで正規化されていないので、この値は
    // 絶対的な W ではなく相対量。それでも
    //   - 損失材料が無ければ 0 になる
    //   - 同じ形状・同じ材料定数なら material id の付け方に依らず一致する
    // という性質は成り立つので、セル毎の材料参照が効いていることを
    // 外から確認できる (data/sample/thermal_material_check.sh)。
    if (io) {
        for (int ifreq = 0; ifreq < NFreq2; ifreq++) {
            const int64_t base_idx = (int64_t)ifreq * NN;
            double psum = 0;
            for (int i = 0; i < Nx; i++) {
            for (int j = 0; j < Ny; j++) {
            for (int k = 0; k < Nz; k++) {
                const double dv = (Xn[i + 1] - Xn[i])
                                * (Yn[j + 1] - Yn[j])
                                * (Zn[k + 1] - Zn[k]);
                psum += P_losses[base_idx + NA(i, j, k)] * dv;
            }
            }
            }
            sprintf(str, "Thermal: dissipated[%d] = %.6e (f=%.6e Hz)", ifreq, psum, Freq2[ifreq]);
            fprintf(fp,     "%s\n", str);
            fprintf(stdout, "%s\n", str);
        }
        fflush(fp);
        fflush(stdout);
    }

    // 熱解析レイヤのメモリ解放 (上の診断で P_losses を読み終えてから)
    free(T);
    free(P_losses);
    free(loss_esgm);
    free(loss_msgm);

    // time steps
    Ntime = itime + converged;

    // メタデータの作成
    hid_t metadata_group_id = H5Gcreate(file_id, "/metadata", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    // 時間に関するメタデータの書き込み(収束条件で終了時の対応)
    //int maxIter = MIN(Solver.maxiter, Ntime);
    //int maxNOut = MIN(Solver.nout, Niter);
    //double time_metadata[1] = {maxIter * Dt};
    double time_metadata[1] = {Solver.maxiter * Dt};
    //dataspace_id = H5Screate_simple(1, count, NULL);
    hsize_t time_count[1] = {1};
    dataspace_id = H5Screate_simple(1, time_count, NULL);
    dataset_id = H5Dcreate(metadata_group_id, "time", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, time_metadata);
    H5Dclose(dataset_id);
    H5Sclose(dataspace_id);

    // グリッドに関するメタデータの書き込み
    //double grid_metadata[3] = {Dx, Dy, Dz};
    //hsize_t grid_count[1] = {3};
    //dataspace_id = H5Screate_simple(1, grid_count, NULL);
    //dataset_id = H5Dcreate(metadata_group_id, "grid", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    //status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, grid_metadata);
    //H5Dclose(dataset_id);
    //H5Sclose(dataspace_id);

    // その他のメタデータの書き込み
/*
    // title, dt, source, fPlanewave, z0, Ni, Nj, Nk, N0, NN
    const char *title = Title;
    dataspace_id = H5Screate(H5S_SCALAR);
    dataset_id = H5Dcreate(metadata_group_id, "title", H5T_C_S1, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dataset_id, H5T_C_S1, H5S_ALL, H5S_ALL, H5P_DEFAULT, title);
    H5Dclose(dataset_id);
    H5Sclose(dataspace_id);

    // 平面波が伝わるときのインピーダンスの値
    //double metadata_values[8] = {Dt, Planewave.z0, Ni, Nj, Nk, N0, NN};
    double metadata_values[8] = {Dt, 0.0, Ni, Nj, Nk, N0, NN};
    hsize_t metadata_count[1] = {8};
    dataspace_id = H5Screate_simple(1, metadata_count, NULL);
    dataset_id = H5Dcreate(metadata_group_id, "metadata_values", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, metadata_values);
    H5Dclose(dataset_id);
    H5Sclose(dataspace_id);

    // 配列に関するメタデータの書き込み (Xn, Yn, Zn, Freq1, Freq2)
    hsize_t array_count[1];
    double *arrays[] = {Xn, Yn, Zn, Freq1, Freq2};
    const char *array_names[] = {"Xn", "Yn", "Zn", "Freq1", "Freq2"};
    size_t array_sizes[] = {Nx + 1, Ny + 1, Nz + 1, NFreq1, NFreq2};

    for (int i = 0; i < 5; i++) {
        array_count[0] = array_sizes[i];
        dataspace_id = H5Screate_simple(1, array_count, NULL);
        dataset_id = H5Dcreate(metadata_group_id, array_names[i], H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, arrays[i]);
        H5Dclose(dataset_id);
        H5Sclose(dataspace_id);
    }
*/
    // Title
    hsize_t title_dims[1] = {256};
    dataspace_id = H5Screate_simple(1, title_dims, NULL);
    dataset_id = H5Dcreate(metadata_group_id, "Title", H5T_NATIVE_CHAR, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dataset_id, H5T_NATIVE_CHAR, H5S_ALL, H5S_ALL, H5P_DEFAULT, Title);
    H5Dclose(dataset_id);
    H5Sclose(dataspace_id);

    // 各種整数型メタデータの書き込み
    struct {
        const char *name;
        void *value;
        hid_t type;
    } metadata[] = {
        {"Nx", &Nx, H5T_NATIVE_INT},
        {"Ny", &Ny, H5T_NATIVE_INT},
        {"Nz", &Nz, H5T_NATIVE_INT},
        {"Ni", &Ni, H5T_NATIVE_INT},
        {"Nj", &Nj, H5T_NATIVE_INT},
        {"Nk", &Nk, H5T_NATIVE_INT},
        {"N0", &N0, H5T_NATIVE_INT},
        {"NN", &NN, H5T_NATIVE_INT64},
        {"NFreq1", &NFreq1, H5T_NATIVE_INT},
        {"NFreq2", &NFreq2, H5T_NATIVE_INT},
        {"NFeed", &NFeed, H5T_NATIVE_INT},
        {"NPoint", &NPoint, H5T_NATIVE_INT},
        {"Niter", &Niter, H5T_NATIVE_INT},
        {"Ntime", &Ntime, H5T_NATIVE_INT},
        {"Solver_maxiter", &Solver.maxiter, H5T_NATIVE_INT},
        {"Solver_nout", &Solver.nout, H5T_NATIVE_INT},
        {"NGline", &NGline, H5T_NATIVE_INT},
        {"IPlanewave", &IPlanewave, H5T_NATIVE_INT}
    };

    for (int i = 0; i < sizeof(metadata) / sizeof(metadata[0]); i++) {
        dataspace_id = H5Screate(H5S_SCALAR);
        dataset_id = H5Dcreate(metadata_group_id, metadata[i].name, metadata[i].type, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Dwrite(dataset_id, metadata[i].type, H5S_ALL, H5S_ALL, H5P_DEFAULT, metadata[i].value);
        H5Dclose(dataset_id);
        H5Sclose(dataspace_id);
    }

    // Dtの書き込み
    dataspace_id = H5Screate(H5S_SCALAR);
    dataset_id = H5Dcreate(metadata_group_id, "Dt", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &Dt);
    H5Dclose(dataset_id);
    H5Sclose(dataspace_id);

    // Planewaveの書き込み
    dataspace_id = H5Screate(H5S_SCALAR);
    dataset_id = H5Dcreate(metadata_group_id, "Planewave", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &Planewave);
    H5Dclose(dataset_id);
    H5Sclose(dataspace_id);

    // 配列データの書き込み
    struct {
        const char *name;
        double *data;
        size_t size;
    } arrays[] = {
        {"Xn", Xn, Nx + 1},
        {"Yn", Yn, Ny + 1},
        {"Zn", Zn, Nz + 1},
        {"Xc", Xc, Nx},
        {"Yc", Yc, Ny},
        {"Zc", Zc, Nz},
        {"Eiter", Eiter, Niter},
        {"Hiter", Hiter, Niter},
        {"VFeed", VFeed, NFeed * (Solver.maxiter + 1)},
        {"IFeed", IFeed, NFeed * (Solver.maxiter + 1)},
        {"VPoint", VPoint, NPoint * (Solver.maxiter + 1)},
        {"Freq1", Freq1, NFreq1},
        {"Freq2", Freq2, NFreq2},
        {"Gline", Gline, NGline * 2 * 3}
    };

    for (int i = 0; i < sizeof(arrays) / sizeof(arrays[0]); i++) {
        hsize_t array_dims[1] = {arrays[i].size};
        dataspace_id = H5Screate_simple(1, array_dims, NULL);
        dataset_id = H5Dcreate(metadata_group_id, arrays[i].name, H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, arrays[i].data);
        H5Dclose(dataset_id);
        H5Sclose(dataspace_id);
    }
    
    // NSurfaceデータの書き込み
    dataspace_id = H5Screate(H5S_SCALAR);
    dataset_id = H5Dcreate(metadata_group_id, "NSurface", H5T_NATIVE_INT, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &NSurface);
    H5Dclose(dataset_id);
    H5Sclose(dataspace_id);

    // Surfaceデータの書き込み
    // surface_t構造体に対応する複合データ型を定義
    hid_t memtype = H5Tcreate(H5T_COMPOUND, sizeof(surface_t));
    H5Tinsert(memtype, "nx", HOFFSET(surface_t, nx), H5T_NATIVE_DOUBLE);
    H5Tinsert(memtype, "ny", HOFFSET(surface_t, ny), H5T_NATIVE_DOUBLE);
    H5Tinsert(memtype, "nz", HOFFSET(surface_t, nz), H5T_NATIVE_DOUBLE);
    H5Tinsert(memtype, "x", HOFFSET(surface_t, x), H5T_NATIVE_DOUBLE);
    H5Tinsert(memtype, "y", HOFFSET(surface_t, y), H5T_NATIVE_DOUBLE);
    H5Tinsert(memtype, "z", HOFFSET(surface_t, z), H5T_NATIVE_DOUBLE);
    H5Tinsert(memtype, "ds", HOFFSET(surface_t, ds), H5T_NATIVE_DOUBLE);

    hsize_t surface_dims[1] = {NSurface};
    dataspace_id = H5Screate_simple(1, surface_dims, NULL);
    dataset_id = H5Dcreate(metadata_group_id, "Surface", memtype, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dataset_id, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, Surface);
    H5Dclose(dataset_id);
    H5Sclose(dataspace_id);
    H5Tclose(memtype);

    // メタデータグループのクローズ
    H5Gclose(metadata_group_id);

    status = H5Fclose(file_id);

    // free
    memfree2();
}
