/*
test_update.c

電磁界の更新 updateE[xyz]() / updateH[xyz]() (sol/update??.c) のテスト。
給電モード (全電界, _f_) と平面波モード (散乱界, _p_) のそれぞれについて、
vector / no-vector の両経路を調べる。

期待値の出所 (Maxwell の積分形と散乱界の定義):
  - Ampère: 辺 (E_a) の 1 ステップの変化は、その辺を右回りに囲む H の
    周回積分を面積で割ったもの
        ΔE_a = c dt ∮H·dl / A        (H は η0 倍で格納)
    E は双対格子の面 (双対セル幅 × 双対セル幅) を囲む。
  - Faraday: ΔH_a = -c dt ∮E·dl / A。H は主格子の面 (セル幅 × セル幅)。
    周回は辺の向き â に対して右ねじ。テストは周回の各辺 (H_c / E_c の位置と
    向き・長さ) を Yee 格子の配置から組み立て、実装の差分の書き方を写さない。
  - 材料: E_new = C1 E_old + C2 ΔE、H_new = D1 H_old + D2 ΔH (係数の値
    自体は sol/setup.c の担当。ここでは「どちらに掛かるか」と材料 id の参照)。
    PEC は C1 = C2 = 0 なので給電モードでは 0。
  - 散乱界 (平面波): 全電界 = 散乱界 + 入射界。全電界は材料の式に、入射界は
    真空の式 (ΔE_i = E_i(t) - E_i(t-dt)) に従うので、差をとると
        E_s(t) = C1 E_s + C2 ΔE_s + C1 E_i(t-dt) + C2 (E_i(t) - E_i(t-dt)) - E_i(t)
    真空 (C1 = C2 = 1) では入射の項が消え、PEC (0, 0) では E_s = -E_i(t)
    (全電界 0)。実装は E_i(t) - E_i(t-dt) を dt dE_i/dt で近似しているので、
    損失のある材料では O(dt^2) の差が出る (パルスを dt に比べて十分
    滑らかにして許容値に収める)。H も同様 (hi 成分、半整数時刻)。
  - 更新範囲: E_a は a 方向にセル [0, N_a)、他は節点 [0, N]。H_a は a 方向に
    節点 [0, N_a]、他はセル [0, N)。範囲外 (ゴースト) は書き換えない。
*/

#include "ofd.h"
#include "ofd_prototype.h"
#include "finc.h"
#include "ofd_test.h"

#include <stdint.h>
#include <stdlib.h>
#include <string.h>

static const double DX[] = {1.0e-3, 1.3e-3, 0.8e-3, 1.6e-3};
static const double DY[] = {0.7e-3, 0.9e-3, 1.2e-3, 1.0e-3, 1.25e-3};
static const double DZ[] = {1.5e-3, 1.0e-3, 0.9e-3, 1.1e-3, 1.3e-3, 0.6e-3};

static int     n_[3];
static double *node_[3], *cell_[3];
static real_t **ef[3], **hf[3];
static id_t   **eid_[3], **hid_[3];
static real_t **k1e[3], **k2e[3], **k1h[3], **k2h[3];
static double  cdt;

static uint32_t rng = 4242;
static double urand(void)
{
	rng = rng * 1664525u + 1013904223u;
	return ((rng >> 8) & 0xffff) / 65536.0 - 0.5;
}

static void make_axis(const double *d, int n, double **nodes, double **cells)
{
	double *xn = (double *)malloc((n + 1) * sizeof(double));
	double *xc = (double *)malloc(n * sizeof(double));
	xn[0] = -0.3e-3;
	for (int i = 0; i < n; i++) {
		xn[i + 1] = xn[i] + d[i];
		xc[i] = (xn[i] + xn[i + 1]) / 2;
	}
	*nodes = xn;
	*cells = xc;
}

static double dual_width(int a, int i)
{
	const int n = n_[a];
	const double *x = node_[a];
	if (i <= 0) return x[1] - x[0];
	if (i >= n) return x[n] - x[n - 1];
	return (x[i + 1] - x[i - 1]) / 2;
}

static double cell_width(int a, int i)
{
	return node_[a][i + 1] - node_[a][i];
}

static void up_setup(void)
{
	n_[0] = Nx = (int)(sizeof(DX) / sizeof(DX[0]));
	n_[1] = Ny = (int)(sizeof(DY) / sizeof(DY[0]));
	n_[2] = Nz = (int)(sizeof(DZ) / sizeof(DZ[0]));
	iMin = jMin = kMin = 0;
	iMax = Nx; jMax = Ny; kMax = Nz;
	make_axis(DX, Nx, &Xn, &Xc);
	make_axis(DY, Ny, &Yn, &Yc);
	make_axis(DZ, Nz, &Zn, &Zc);
	node_[0] = Xn; node_[1] = Yn; node_[2] = Zn;
	cell_[0] = Xc; cell_[1] = Yc; cell_[2] = Zc;

	Dt = 1.0e-12;
	cdt = C * Dt;

	Nk = 1;
	Nj = (kMax - kMin) + 2 + 1;
	Ni = ((jMax - jMin) + 2 + 1) * Nj;
	N0 = -((iMin - 1) * Ni + (jMin - 1) * Nj + (kMin - 1) * Nk);
	NN = NA(iMax + 1, jMax + 1, kMax + 1) + 1;

	real_t **E[3] = {&Ex, &Ey, &Ez}, **H[3] = {&Hx, &Hy, &Hz};
	id_t **IE[3] = {&iEx, &iEy, &iEz}, **IH[3] = {&iHx, &iHy, &iHz};
	real_t **K1E[3] = {&K1Ex, &K1Ey, &K1Ez}, **K2E[3] = {&K2Ex, &K2Ey, &K2Ez};
	real_t **K1H[3] = {&K1Hx, &K1Hy, &K1Hz}, **K2H[3] = {&K2Hx, &K2Hy, &K2Hz};
	for (int a = 0; a < 3; a++) {
		*E[a] = (real_t *)calloc(NN, sizeof(real_t));
		*H[a] = (real_t *)calloc(NN, sizeof(real_t));
		*IE[a] = (id_t *)calloc(NN, sizeof(id_t));
		*IH[a] = (id_t *)calloc(NN, sizeof(id_t));
		*K1E[a] = (real_t *)calloc(NN, sizeof(real_t));
		*K2E[a] = (real_t *)calloc(NN, sizeof(real_t));
		*K1H[a] = (real_t *)calloc(NN, sizeof(real_t));
		*K2H[a] = (real_t *)calloc(NN, sizeof(real_t));
		ef[a] = E[a]; hf[a] = H[a]; eid_[a] = IE[a]; hid_[a] = IH[a];
		k1e[a] = K1E[a]; k2e[a] = K2E[a]; k1h[a] = K1H[a]; k2h[a] = K2H[a];
	}

	/* 節点は双対セル幅、セルはセル幅 (端の節点は境界セル幅) */
	real_t **rn[3] = {&RXn, &RYn, &RZn}, **rc[3] = {&RXc, &RYc, &RZc};
	for (int a = 0; a < 3; a++) {
		*rn[a] = (real_t *)malloc((n_[a] + 1) * sizeof(real_t));
		*rc[a] = (real_t *)malloc(n_[a] * sizeof(real_t));
		for (int i = 0; i <= n_[a]; i++) (*rn[a])[i] = (real_t)(cdt / dual_width(a, i));
		for (int i = 0; i <  n_[a]; i++) (*rc[a])[i] = (real_t)(cdt / cell_width(a, i));
	}

	/* 0=真空, 1=PEC, 2,3=損失のある材料 (C1≠C2, D1≠D2) */
	NMaterial = 4;
	C1 = (real_t *)malloc(NMaterial * sizeof(real_t));
	C2 = (real_t *)malloc(NMaterial * sizeof(real_t));
	D1 = (real_t *)malloc(NMaterial * sizeof(real_t));
	D2 = (real_t *)malloc(NMaterial * sizeof(real_t));
	const real_t c1[4] = {1, 0, 0.8f, 0.6f}, c2[4] = {1, 0, 0.35f, 0.2f};
	const real_t d1[4] = {1, 0, 0.9f, 0.7f}, d2[4] = {1, 0, 0.5f, 0.45f};
	for (int m = 0; m < 4; m++) { C1[m] = c1[m]; C2[m] = c2[m]; D1[m] = d1[m]; D2[m] = d2[m]; }

	/* 材料 id はセルごとにばらばら。vector 経路の係数配列も同じ材料から作る
	   (sol/setup_vector.c と同じく K1 = C1[m], K2 = C2[m]。PEC は 0) */
	for (int a = 0; a < 3; a++) {
		for (int64_t q = 0; q < NN; q++) {
			const id_t me = (id_t)((uint32_t)(urand() * 1e6 + 1e6) % 4);
			const id_t mh = (id_t)((uint32_t)(urand() * 1e6 + 1e6) % 4);
			(*eid_[a])[q] = me;
			(*hid_[a])[q] = mh;
			(*k1e[a])[q] = C1[me]; (*k2e[a])[q] = C2[me];
			(*k1h[a])[q] = D1[mh]; (*k2h[a])[q] = D2[mh];
		}
	}

	/* 平面波 : 斜め入射、パルス幅は dt の約 30 倍 (ai dt = 0.03) */
	memset(&Planewave, 0, sizeof(Planewave));
	const double th = 50 * DTOR, ph = 35 * DTOR;
	Planewave.ri[0] = -sin(th) * cos(ph);
	Planewave.ri[1] = -sin(th) * sin(ph);
	Planewave.ri[2] = -cos(th);
	Planewave.ei[0] = 0.6;  Planewave.ei[1] = -0.3; Planewave.ei[2] = 0.74;
	Planewave.hi[0] = -0.2; Planewave.hi[1] = 0.85; Planewave.hi[2] = 0.49;
	Planewave.ai = 0.03 / Dt;
	/* 基準点を格子の中ほどに置き、格子内の遅延を ±7 dt 程度に収める */
	Planewave.r0[0] = 2.0e-3; Planewave.r0[1] = 2.5e-3; Planewave.r0[2] = 3.0e-3;
	WaveAmp = 0;
}

static void up_teardown(void)
{
	for (int a = 0; a < 3; a++) {
		free(*ef[a]); free(*hf[a]); free(*eid_[a]); free(*hid_[a]);
		free(*k1e[a]); free(*k2e[a]); free(*k1h[a]); free(*k2h[a]);
		*ef[a] = *hf[a] = NULL; *eid_[a] = *hid_[a] = NULL;
		*k1e[a] = *k2e[a] = *k1h[a] = *k2h[a] = NULL;
	}
	free(RXn); free(RYn); free(RZn); free(RXc); free(RYc); free(RZc);
	RXn = RYn = RZn = RXc = RYc = RZc = NULL;
	free(C1); free(C2); free(D1); free(D2);
	C1 = C2 = D1 = D2 = NULL;
	free(Xn); free(Yn); free(Zn); free(Xc); free(Yc); free(Zc);
	Xn = Yn = Zn = Xc = Yc = Zc = NULL;
	NMaterial = 0;
	memset(&Planewave, 0, sizeof(Planewave));
	IPlanewave = 0;
	NFeed = 0;
	VECTOR = 0;
}

static void fill_random(void)
{
	for (int a = 0; a < 3; a++) {
		for (int64_t q = 0; q < NN; q++) {
			(*ef[a])[q] = (real_t)urand();
			(*hf[a])[q] = (real_t)urand();
		}
	}
}

static real_t val(real_t **f, const int idx[3])
{
	return (*f)[NA(idx[0], idx[1], idx[2])];
}

/* 成分 (isE, a) の格子点 idx の座標 */
static void position(int isE, int a, const int idx[3], double r[3])
{
	for (int b = 0; b < 3; b++) {
		const int on_cell = isE ? (b == a) : (b != a);
		r[b] = on_cell ? cell_[b][idx[b]] : node_[b][idx[b]];
	}
}

/* 更新範囲 : E_a は a 方向セル [0,N_a)、他は節点 [0,N]。H_a は逆 */
static int in_range(int isE, int a, const int idx[3])
{
	for (int b = 0; b < 3; b++) {
		const int on_cell = isE ? (b == a) : (b != a);
		const int hi = on_cell ? (n_[b] - 1) : n_[b];
		if ((idx[b] < 0) || (idx[b] > hi)) return 0;
	}
	return 1;
}

/* Maxwell の積分形による増分 (E: +c dt ∮H·dl/A、H: -c dt ∮E·dl/A)。
   src は更新前の相手の場 (E の更新なら H) */
static double maxwell_delta(int isE, int a, const int idx[3], real_t **src[3])
{
	double loop = 0, area = 1;
	for (int s = 1; s <= 2; s++) {
		const int b = (a + s) % 3;           /* 辺に垂直な軸 (周回の一方の辺の向き) */
		const int c = 3 - a - b;             /* そこにある相手の場の成分 */
		/* â × b̂ = +ĉ (b = a+1) / -ĉ (b = a+2) : +b 側の辺は ±ĉ 向きに回る */
		const double sgn = (b == (a + 1) % 3) ? +1.0 : -1.0;
		int plus[3] = {idx[0], idx[1], idx[2]};
		int minus[3] = {idx[0], idx[1], idx[2]};
		double len;
		if (isE) {
			/* H_c はセル中心 : +b 側は添字 idx_b、-b 側は idx_b - 1。
			   周回の辺の長さは c 方向の双対セル幅 */
			minus[b] -= 1;
			len = dual_width(c, idx[c]);
			area *= dual_width(b, idx[b]);
		}
		else {
			/* E_c は節点 : +b 側は idx_b + 1、-b 側は idx_b。長さは c 方向のセル幅 */
			plus[b] += 1;
			len = cell_width(c, idx[c]);
			area *= cell_width(b, idx[b]);
		}
		loop += sgn * (val(src[c], plus) - val(src[c], minus)) * len;
	}
	return (isE ? +1 : -1) * cdt * loop / area;
}

/* 1 成分を更新して全点を期待値と比べる。戻り値 : 最大誤差 / 範囲外の書き換え数 */
static void check_component(int isE, int a, int planewave, double t, double *maxerr, int *touched)
{
	real_t **f = isE ? ef[a] : hf[a];
	real_t **src[3];
	for (int c = 0; c < 3; c++) src[c] = isE ? hf[c] : ef[c];
	id_t **ids = isE ? eid_[a] : hid_[a];
	const real_t *A1 = isE ? C1 : D1;
	const real_t *A2 = isE ? C2 : D2;
	const double comp = isE ? Planewave.ei[a] : Planewave.hi[a];

	real_t *before = (real_t *)malloc(NN * sizeof(real_t));
	memcpy(before, *f, NN * sizeof(real_t));

	void (*upE[3])(double) = {updateEx, updateEy, updateEz};
	void (*upH[3])(double) = {updateHx, updateHy, updateHz};
	(isE ? upE[a] : upH[a])(t);

	for (int i = -1; i <= Nx + 1; i++) {
	for (int j = -1; j <= Ny + 1; j++) {
	for (int k = -1; k <= Nz + 1; k++) {
		const int idx[3] = {i, j, k};
		const int64_t q = NA(i, j, k);
		if ((q < 0) || (q >= NN)) continue;
		if (!in_range(isE, a, idx)) {
			*touched += ((*f)[q] != before[q]);
			continue;
		}
		const id_t m = (*ids)[q];
		const double d = maxwell_delta(isE, a, idx, src);
		double expect = A1[m] * before[q] + A2[m] * d;
		if (planewave) {
			/* 入射界 (finc は test_finc で検証済み) を格子点の位置で */
			double r[3];
			real_t fnew, fold, dfi;
			position(isE, a, idx, r);
			finc(r[0], r[1], r[2], t,      Planewave.r0, Planewave.ri, comp, Planewave.ai, Dt, &fnew, &dfi);
			finc(r[0], r[1], r[2], t - Dt, Planewave.r0, Planewave.ri, comp, Planewave.ai, Dt, &fold, &dfi);
			expect += A1[m] * fold + A2[m] * (fnew - fold) - fnew;
		}
		*maxerr = fmax(*maxerr, fabs((*f)[q] - expect));
	}
	}
	}

	free(before);
}

static void run_mode(int planewave, int vector)
{
	NFeed = planewave ? 0 : 1;
	IPlanewave = planewave;
	VECTOR = vector;

	/* 入射パルスが格子を通過中の時刻 (格子内で at ≈ 0.5〜0.9、f ≈ 0.4) */
	const double tE = 24 * Dt, tH = 23.5 * Dt;
	const char *name[2] = {"E", "H"};
	for (int isE = 1; isE >= 0; isE--) {
		for (int a = 0; a < 3; a++) {
			fill_random();
			double err = 0;
			int touched = 0;
			check_component(isE, a, planewave, isE ? tE : tH, &err, &touched);
			/* 給電モードは厳密 (float の丸め)。平面波モードは dt dE_i/dt 近似の O(dt^2) */
			const double tol = planewave ? 2e-3 : 1e-5;
			CHECK(err < tol);
			CHECK(touched == 0);
			if ((err >= tol) || touched) {
				printf("    %s%c (%s, %s) : max error %.3e, touched %d\n", name[1 - isE], 'x' + a,
					planewave ? "planewave" : "feed", vector ? "vector" : "no-vector", err, touched);
			}
		}
	}
}

/* 入射界が本当に効いていること (許容値に対して十分大きい) を確かめる */
static void check_incident_scale(void)
{
	real_t fi, dfi;
	double fmax_ = 0;
	for (int i = 0; i < Nx; i++) {
		finc(Xc[i], Yn[2], Zn[3], 24 * Dt, Planewave.r0, Planewave.ri, 1.0, Planewave.ai, Dt, &fi, &dfi);
		fmax_ = fmax(fmax_, fabs(fi));
	}
	CHECK(fmax_ > 0.3);
}

void test_update(void)
{
	const double dt0 = Dt;

	up_setup();
	check_incident_scale();
	run_mode(0, 0);
	run_mode(0, 1);
	run_mode(1, 0);
	run_mode(1, 1);
	up_teardown();

	Dt = dt0;
}
