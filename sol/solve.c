#include "ofd.h"
#include "complex.h"
#include "ofd_prototype.h"
#include "ev.h"
#include "finc.h"

#include "hdf5.h"
#include "ofd_hdf5.h"

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
    // HDF5 出力は sol/outputHdf5.c に集約している
    // (従来は sol/mpi/cuda/cuda_mpi の 4 つの solve に同じコードが重複していた)

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

    // 発熱密度 (ログの Thermal: dissipated と HDF5 の /loss/P_loss に使う)
    //
    // 温度分布そのものは計算しない。以前は 3 次元熱拡散 (updateTemperature) を
    // 毎ステップ回していたが、結果はどこにも出力されず捨てられていた上に、
    // 数値的にも物理的にも意味を持っていなかったので計算ごと外した:
    //
    //  - 時間刻みが FDTD の Dt (光学解析では ~1e-16 s) のままだった。
    //    T += Dt * P_loss は 1 ステップあたり ~7e-18 しか足されず、初期温度 20 の
    //    倍精度分解能 (20 * 2.2e-16 = 4.4e-15) に完全に吸収される。実際
    //    thermal_slab.ofd でも最後まで全セルが初期値のまま (T = [20, 20]) だった。
    //  - 温度を周波数ごとに持っていた。温度は全周波数の発熱の合計で決まる
    //    単一のスカラー場なので、NFreq2 面に分けた時点で物理が合わない。
    //  - ラプラシアンが一様格子 (平均 Dx/Dy/Dz) 前提だった。OpenFDTD の
    //    メッシュは不等間隔を許すので、そのままでは使えない。
    //
    // 意味のある温度分布を出すには、熱拡散用の別の時間刻み (fs ではなく
    // us〜s のオーダー)、入射電力からの絶対値換算 (P_loss は相対値)、材料ごとの
    // 熱拡散係数、不等間隔格子のラプラシアンが要る。ここで残す P_losses が
    // その入力になる。
    double *P_losses = (double *)malloc(NFreq2 * NN * sizeof(double));
    if (P_losses == NULL) {
        fprintf(stderr, "*** power loss array malloc error (NFreq2=%d NN=%zu)\n", NFreq2, (size_t)NN);
        exit(1);
    }
    memset(P_losses, 0, NFreq2 * NN * sizeof(double));

    // HDF5ファイルの作成
    hdf5_open(1);

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

                // HDF5 : 瞬時値スナップショット (sol/outputHdf5.c)
                hdf5_write_snapshot(itime, t, Ex, Ey, Ez, Hx, Hy, Hz);
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

    // HDF5 : 損失電力密度 (P_losses を解放する前に書く)
    hdf5_write_loss(P_losses);

    // 熱解析レイヤのメモリ解放 (上の診断で P_losses を読み終えてから)
    free(P_losses);
    free(loss_esgm);
    free(loss_msgm);

    // time steps
    Ntime = itime + converged;

    // HDF5 : 周波数領域の最終結果・収束履歴・メタデータ (sol/outputHdf5.c)
    hdf5_write_freqdomain(
        cEx_r, cEx_i, cEy_r, cEy_i, cEz_r, cEz_i,
        cHx_r, cHx_i, cHy_r, cHy_i, cHz_r, cHz_i);
    hdf5_write_convergence(Niter, Eiter, Hiter);
    hdf5_close();

    // free
    memfree2();
}
