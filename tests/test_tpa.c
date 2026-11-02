/*
test_tpa.c

二光子吸収 (TPA) の非線形減衰 updateTpa() / setupTpa() (sol/updateTpa.c) のテスト。
物理の通し検証は data/sample/tpa_slab_check.sh (解析解 T = 1/(1+βI0L))。
ここは 1 点の減衰の式と、格子上の扱い (成分の補間・散乱界) を押さえる。

期待値の出所:
  - 減衰の大きさはエネルギー収支から決める (実装の調和平衡の導出は使わない)。
    TPA は dI/dz = -β I² なので、速さ v = c/n で進む CW 波のエネルギーは
    1 ステップに  ΔU/U = -v β I dt = -(1/2) β ε0 c² E0² dt  減る
    (I = (1/2) n ε0 c E0²)。減衰は電界だけに掛かり、進行波では電界と磁界の
    エネルギーが等しいので、電界のエネルギーの減少 ⟨E² - E'²⟩ の半分が
    全エネルギーの減少になる:
        ⟨E² - E'²⟩ / (2 ⟨E²⟩) = (1/2) β ε0 c² E0² dt
    1 点に E(t) = E0 cos(ωt) を毎ステップ与えて減衰させ、1 周期平均で比べる。
    この関係は εr に依らない (材料を変えても同じ) こと。
  - |E|² は各成分の位置での値: 自成分の 2 乗 + 他の 2 成分を周囲 4 点で
    平均した値の 2 乗。周囲 4 点 (例: Ex 位置の Ey は x 方向 i, i+1 と
    y 方向 j-1, j) だけに値を置き、他は大きな値にして参照位置を検査する。
  - 平面波 (散乱界) では全電界 = 散乱界 + 入射界を減衰させ、散乱界に戻す。
  - β = 0 の材料のセルは書き換えない。
  - `if (b <= 0) continue` を `b < 0` にする変更は、β = 0 なら γ = 0 で値が
    変わらない (スキップは計算の省略だけ) 等価な書き換えなので検出対象にしない。
*/

#include "ofd.h"
#include "ofd_prototype.h"
#include "finc.h"
#include "ofd_test.h"

#include <stdlib.h>
#include <string.h>

#define N 4

static void tp_setup(void)
{
	Nx = Ny = Nz = N;
	iMin = jMin = kMin = 0;
	iMax = jMax = kMax = N;
	Xn = (double *)malloc((N + 1) * sizeof(double));
	Yn = (double *)malloc((N + 1) * sizeof(double));
	Zn = (double *)malloc((N + 1) * sizeof(double));
	Xc = (double *)malloc(N * sizeof(double));
	Yc = (double *)malloc(N * sizeof(double));
	Zc = (double *)malloc(N * sizeof(double));
	for (int i = 0; i <= N; i++) {
		Xn[i] = 1e-7 * i;
		Yn[i] = 1.2e-7 * i;
		Zn[i] = 0.9e-7 * i;
	}
	for (int i = 0; i < N; i++) {
		Xc[i] = (Xn[i] + Xn[i + 1]) / 2;
		Yc[i] = (Yn[i] + Yn[i + 1]) / 2;
		Zc[i] = (Zn[i] + Zn[i + 1]) / 2;
	}
	Nk = 1;
	Nj = (kMax - kMin) + 2 + 1;
	Ni = ((jMax - jMin) + 2 + 1) * Nj;
	N0 = -((iMin - 1) * Ni + (jMin - 1) * Nj + (kMin - 1) * Nk);
	NN = NA(iMax + 1, jMax + 1, kMax + 1) + 1;
	Ex = (real_t *)calloc(NN, sizeof(real_t));
	Ey = (real_t *)calloc(NN, sizeof(real_t));
	Ez = (real_t *)calloc(NN, sizeof(real_t));
	iEx = (id_t *)calloc(NN, sizeof(id_t));
	iEy = (id_t *)calloc(NN, sizeof(id_t));
	iEz = (id_t *)calloc(NN, sizeof(id_t));

	/* 0=真空, 1=PEC, 2=TPA 材料 (εr=12.1), 3=TPA 材料 (εr=2, 別の β), 4=TPA なしの誘電体 */
	NMaterial = 5;
	Material = (material_t *)calloc(NMaterial, sizeof(material_t));
	for (int m = 0; m < NMaterial; m++) { Material[m].type = 1; Material[m].epsr = Material[m].amur = 1; }
	Material[2].epsr = 12.1;
	Material[3].epsr = 2;
	Material[4].epsr = 4;
	NTpa = 2;
	Tpa = (tpa_t *)calloc(NTpa, sizeof(tpa_t));
	Tpa[0].m = 2; Tpa[0].beta = 2000e-11;   /* 2000 cm/GW */
	Tpa[1].m = 3; Tpa[1].beta = 500e-11;
	TpaBeta = NULL;
	setupTpa();

	Dt = 1e-17;
	IPlanewave = 0;
	WaveAmp = 0;
	memset(&Planewave, 0, sizeof(Planewave));
}

static void tp_teardown(void)
{
	free(Xn); free(Yn); free(Zn); free(Xc); free(Yc); free(Zc);
	free(Ex); free(Ey); free(Ez); free(iEx); free(iEy); free(iEz);
	free(Material); free(Tpa); free(TpaBeta);
	Xn = Yn = Zn = Xc = Yc = Zc = NULL;
	Ex = Ey = Ez = NULL;
	iEx = iEy = iEz = NULL;
	Material = NULL; Tpa = NULL; TpaBeta = NULL;
	NMaterial = 0; NTpa = 0;
	IPlanewave = 0;
	memset(&Planewave, 0, sizeof(Planewave));
}

static void test_table(void)
{
	CHECK(TpaBeta != NULL);
	CHECK(TpaBeta[0] == 0 && TpaBeta[1] == 0 && TpaBeta[4] == 0);
	CHECK(TpaBeta[2] == 2000e-11);
	CHECK(TpaBeta[3] == 500e-11);
}

/* 1 点 (Ex) に E0 cos(ωt) を与え続け、1 周期平均のエネルギー減少を比べる */
static void test_energy_balance(void)
{
	const int i = 1, j = 2, k = 2;
	for (id_t m = 2; m <= 3; m++) {
		memset(Ex, 0, NN * sizeof(real_t));
		memset(Ey, 0, NN * sizeof(real_t));
		memset(Ez, 0, NN * sizeof(real_t));
		memset(iEx, 0, NN * sizeof(id_t));
		iEx[NA(i, j, k)] = m;
		const double beta = TpaBeta[m];
		const double E0 = 3e7;              /* γ ≈ 1e-4 程度 (一次の展開が効く) */
		const int nper = 200;               /* 1 周期 200 ステップ */
		double sum_e2 = 0, sum_loss = 0;
		for (int s = 0; s < 4 * nper; s++) {
			const double e = E0 * cos(2 * PI * s / nper + 0.3);
			EX(i, j, k) = (real_t)e;
			updateTpa((s + 1) * Dt);
			const double e1 = EX(i, j, k);
			sum_e2 += e * e;
			sum_loss += e * e - e1 * e1;
		}
		const double measured = sum_loss / (2 * sum_e2);
		const double expect = 0.5 * beta * EPS0 * C * C * E0 * E0 * Dt;
		CHECK(expect > 1e-6 && expect < 1e-3);
		/* real_t (float) で E を持つので、相対 1e-2 で比べる */
		CHECK_NEAR(measured, expect, 1e-2 * expect);
		if (fabs(measured - expect) > 1e-2 * expect) {
			printf("    material %d : measured %.4e, expected %.4e\n", (int)m, measured, expect);
		}
	}
}

/* 他成分の 4 点平均 : 正しい 4 点にだけ値を置き、他は大きな値 */
static void test_colocation(void)
{
	real_t *E[3] = {Ex, Ey, Ez};
	id_t *ID[3] = {iEx, iEy, iEz};
	const int idx[3] = {1, 2, 2};
	const double own = 1e8, val[3] = {2e7, 4e7, 3e7};
	const double cf = (2.0 / 3.0) * EPS0 * C * C * Dt;   /* 上のエネルギー収支で確かめた係数 */

	for (int a = 0; a < 3; a++) {
		for (int c = 0; c < 3; c++) {
			for (int64_t q = 0; q < NN; q++) E[c][q] = (real_t)9e8;
			memset(ID[c], 0, NN * sizeof(id_t));
		}
		ID[a][NA(idx[0], idx[1], idx[2])] = 2;
		E[a][NA(idx[0], idx[1], idx[2])] = (real_t)own;
		/* 他の成分 c の 4 点 : a 方向に idx_a, idx_a + 1 (E_a はセル中心、E_c は節点)、
		   c 方向に idx_c - 1, idx_c (E_c はセル中心、E_a 位置は節点)。
		   4 点は互いに違う値 (平均 val[c]) にする。同じ値だと 4 点のうち
		   2 点だけの平均でも同じになり、取りこぼしを検出できない
		   (重み 0.1/0.5/1.2/2.2 は、どの 2 点・3 点の組の平均も 1 にならない) */
		const double w[4] = {0.1, 0.5, 1.2, 2.2};
		for (int c = 0; c < 3; c++) {
			if (c == a) continue;
			int nw = 0;
			for (int da = 0; da <= 1; da++) {
				for (int dc = -1; dc <= 0; dc++) {
					int q[3] = {idx[0], idx[1], idx[2]};
					q[a] += da;
					q[c] += dc;
					E[c][NA(q[0], q[1], q[2])] = (real_t)(val[c] * w[nw++]);
				}
			}
		}
		real_t *copy[3];
		for (int c = 0; c < 3; c++) {
			copy[c] = (real_t *)malloc(NN * sizeof(real_t));
			memcpy(copy[c], E[c], NN * sizeof(real_t));
		}
		updateTpa(Dt);

		double e2 = own * own;
		for (int c = 0; c < 3; c++) if (c != a) e2 += val[c] * val[c];
		const double gamma = cf * TpaBeta[2] * e2;
		CHECK(gamma > 1e-3);
		const int64_t q0 = NA(idx[0], idx[1], idx[2]);
		CHECK_NEAR(E[a][q0], own / (1 + gamma), 1e-6 * own);
		/* β = 0 の点 (対象の 1 点以外すべて) は書き換えない */
		int changed = 0;
		for (int c = 0; c < 3; c++) {
			for (int64_t q = 0; q < NN; q++) {
				if ((c == a) && (q == q0)) continue;
				changed += (E[c][q] != copy[c][q]);
			}
			free(copy[c]);
		}
		CHECK(changed == 0);
	}
}

/* 平面波 : 全電界を減衰させ、散乱界に戻す */
static void test_planewave(void)
{
	const int i = 2, j = 1, k = 2;
	memset(Ex, 0, NN * sizeof(real_t));
	memset(Ey, 0, NN * sizeof(real_t));
	memset(Ez, 0, NN * sizeof(real_t));
	memset(iEx, 0, NN * sizeof(id_t));
	memset(iEy, 0, NN * sizeof(id_t));
	memset(iEz, 0, NN * sizeof(id_t));
	iEz[NA(i, j, k)] = 2;

	/* CW 波源 (waveamp) にして入射界を散乱界と同程度 (1e8 V/m) にする。
	   振幅 1 のパルスだと float の丸めに埋もれて入射界の扱いを検査できない */
	IPlanewave = 1;
	Planewave.ri[0] = 0; Planewave.ri[1] = 0; Planewave.ri[2] = -1;
	Planewave.ei[0] = 0.6; Planewave.ei[1] = 0; Planewave.ei[2] = 0.8;
	Planewave.r0[0] = Planewave.r0[1] = 0; Planewave.r0[2] = 2e-7;
	WaveAmp = 1e8;
	WaveOmega = 2 * PI / (40 * Dt);

	const double es = 3e7;
	EZ(i, j, k) = (real_t)es;
	const double t = 307 * Dt;           /* ランプ (5 周期 = 200 dt) の後 */
	updateTpa(t);

	real_t fi, dfi;
	finc(Xn[i], Yn[j], Zc[k], t, Planewave.r0, Planewave.ri, 1, Planewave.ai, Dt, &fi, &dfi);
	CHECK(fabs(fi) > 3e7);
	const double ex = 0 + fi * Planewave.ei[0];
	const double ez = es + fi * Planewave.ei[2];
	const double cf = (2.0 / 3.0) * EPS0 * C * C * Dt;
	const double gamma = cf * TpaBeta[2] * (ex * ex + ez * ez);
	const double expect = ez / (1 + gamma) - fi * Planewave.ei[2];
	CHECK(gamma * fabs(ez) > 100 * 1e-6 * es);   /* 減衰が許容値より十分大きい */
	CHECK_NEAR(EZ(i, j, k), expect, 1e-6 * es);
	WaveAmp = 0;
	WaveOmega = 0;
	IPlanewave = 0;
}

void test_tpa(void)
{
	const double dt0 = Dt;
	tp_setup();
	test_table();
	test_energy_balance();
	test_colocation();
	test_planewave();
	tp_teardown();
	Dt = dt0;
}
