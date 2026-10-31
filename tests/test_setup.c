/*
test_setup.c

セットアップ setup() (sol/setup.c) のテスト。小さな問題 (Mur、物体なし) で
setup() を通しで呼び、決まる量を物理の関係と比べる。

期待値の出所:
  - 時間刻み: 3 次元 Yee 格子の安定条件 (Courant)
        c dt <= 1 / √(1/dx² + 1/dy² + 1/dz²)
    の等号を、各軸の最小セル幅で。最小セルを各軸の末尾に置き、最小値の
    探索範囲 (端のセルの取りこぼし) を検出する。
  - パルス幅: Tw = 1.27 / f0、f0 は frequency1 の中央 (無ければ frequency2、
    それも無ければ 1/(20 dt))。給電波形のスペクトルの山を f0 付近に置く値。
  - 材料係数: 損失を陰的に入れた Ampère の差分式
        εr (E(n+1) - E(n)) + (σ dt / ε0) E(n+1) = X       (X = c dt ∇×H)
    を E(n+1) = C1 E(n) + C2 X が任意の E(n), X について満たすこと
    (恒等式として確かめる)。H は μr と σm dt / μ0。真空 C1 = C2 = 1、
    PEC 0。回転 0 の導体中の減衰率 C1 ≈ exp(-σ dt / (ε0 εr)) (O((σdt/ε)²))。
    分散性材料は D の保存 (tests/test_dispersion.c) と同じ C1 = ε∞/(ε∞+ξ0)。
  - メッシュ係数: RXn は節点の双対セル幅 (端は境界セル幅)、RXc はセル幅で
    c dt を割ったもの。
  - 既定値の上書き: 入力で dt / Tw を与えたら setup() はそれを変えない。
*/

#include "ofd.h"
#include "ofd_prototype.h"
#include "ofd_test.h"

#include <stdlib.h>
#include <string.h>

/* 最小セルは各軸の末尾 (0.5 / 0.4 / 0.3 mm) */
static const double DX[] = {1.0e-3, 1.3e-3, 0.8e-3, 0.5e-3};
static const double DY[] = {0.7e-3, 0.9e-3, 1.2e-3, 1.0e-3, 0.4e-3};
static const double DZ[] = {1.5e-3, 1.0e-3, 0.9e-3, 1.1e-3, 1.3e-3, 0.3e-3};

static int     n_[3];
static double *node_[3];

static void make_axis(const double *d, int n, double **nodes, double **cells)
{
	double *xn = (double *)malloc((n + 1) * sizeof(double));
	double *xc = (double *)malloc(n * sizeof(double));
	xn[0] = 0;
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

static void st_setup(void)
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

	iABC = 0;
	VECTOR = 0;
	NGeometry = 0;
	Geometry = NULL;
	NFeed = 0;
	IPlanewave = 0;
	NFreq2 = 0;

	Nk = 1;
	Nj = (kMax - kMin) + 2 + 1;
	Ni = ((jMax - jMin) + 2 + 1) * Nj;
	N0 = -((iMin - 1) * Ni + (jMin - 1) * Nj + (kMin - 1) * Nk);
	NN = NA(iMax + 1, jMax + 1, kMax + 1) + 1;

	iEx = (id_t *)calloc(NN, sizeof(id_t));
	iEy = (id_t *)calloc(NN, sizeof(id_t));
	iEz = (id_t *)calloc(NN, sizeof(id_t));
	iHx = (id_t *)calloc(NN, sizeof(id_t));
	iHy = (id_t *)calloc(NN, sizeof(id_t));
	iHz = (id_t *)calloc(NN, sizeof(id_t));

	/* 0=真空, 1=PEC, 2=損失誘電体, 3=磁気損失磁性体, 4=分散性 (Debye) */
	NMaterial = 5;
	Material = (material_t *)calloc(NMaterial, sizeof(material_t));
	for (int m = 0; m < NMaterial; m++) {
		Material[m].type = 1;
		Material[m].epsr = Material[m].amur = 1;
	}
	Material[2].epsr = 4.0;  Material[2].esgm = 0.5;
	Material[3].amur = 2.5;  Material[3].msgm = 30000;
	Material[4].type = 2;    Material[4].einf = 2; Material[4].be = 4e10; Material[4].ce = 1e10;

	C1 = (real_t *)calloc(NMaterial, sizeof(real_t));
	C2 = (real_t *)calloc(NMaterial, sizeof(real_t));
	D1 = (real_t *)calloc(NMaterial, sizeof(real_t));
	D2 = (real_t *)calloc(NMaterial, sizeof(real_t));

	RXn = (real_t *)calloc(Nx + 1, sizeof(real_t));
	RYn = (real_t *)calloc(Ny + 1, sizeof(real_t));
	RZn = (real_t *)calloc(Nz + 1, sizeof(real_t));
	RXc = (real_t *)calloc(Nx, sizeof(real_t));
	RYc = (real_t *)calloc(Ny, sizeof(real_t));
	RZc = (real_t *)calloc(Nz, sizeof(real_t));
}

static void st_alloc_mur(void)
{
	numMurHx = numMurHy = numMurHz = 0;
	setupMurHx(0);
	setupMurHy(0);
	setupMurHz(0);
	fMurHx = (mur_t *)calloc(numMurHx, sizeof(mur_t));
	fMurHy = (mur_t *)calloc(numMurHy, sizeof(mur_t));
	fMurHz = (mur_t *)calloc(numMurHz, sizeof(mur_t));
}

static void st_free_run(void)
{
	free(fMurHx); free(fMurHy); free(fMurHz);
	fMurHx = fMurHy = fMurHz = NULL;
	free(mDispersionEx); free(mDispersionEy); free(mDispersionEz);
	free(DispersionEx); free(DispersionEy); free(DispersionEz);
	mDispersionEx = mDispersionEy = mDispersionEz = NULL;
	DispersionEx = DispersionEy = DispersionEz = NULL;
	numDispersionEx = numDispersionEy = numDispersionEz = 0;
}

static void st_teardown(void)
{
	st_free_run();
	free(Xn); free(Yn); free(Zn); free(Xc); free(Yc); free(Zc);
	free(iEx); free(iEy); free(iEz); free(iHx); free(iHy); free(iHz);
	free(Material); free(C1); free(C2); free(D1); free(D2);
	free(RXn); free(RYn); free(RZn); free(RXc); free(RYc); free(RZc);
	Xn = Yn = Zn = Xc = Yc = Zc = NULL;
	iEx = iEy = iEz = iHx = iHy = iHz = NULL;
	Material = NULL; C1 = C2 = D1 = D2 = NULL;
	RXn = RYn = RZn = RXc = RYc = RZc = NULL;
	NMaterial = 0;
	NFreq1 = 0;
	free(Freq1);
	Freq1 = NULL;
}

static void run_setup(double dt, double tw)
{
	Dt = dt;
	Tw = tw;
	st_free_run();
	st_alloc_mur();
	setup();
}

void test_setup(void)
{
	const double dt0 = Dt, tw0 = Tw;

	st_setup();
	NFreq1 = 3;
	Freq1 = (double *)malloc(NFreq1 * sizeof(double));
	Freq1[0] = 2e9; Freq1[1] = 3e9; Freq1[2] = 5e9;

	/* 既定の dt / Tw */
	run_setup(0, 0);

	/* Courant (各軸の最小セル : 末尾) */
	const double dtc = 1 / (C * sqrt(1 / (0.5e-3 * 0.5e-3) + 1 / (0.4e-3 * 0.4e-3) + 1 / (0.3e-3 * 0.3e-3)));
	CHECK_NEAR(Dt, dtc, 1e-12 * dtc);

	/* パルス幅 : f0 = (2 + 5)/2 GHz */
	CHECK_NEAR(Tw, 1.27 / 3.5e9, 1e-20);

	/* 材料係数 : 陰的損失の Ampère / Faraday を恒等式として満たす */
	const double cdt = C * Dt;
	const double en = 0.37, xx = -0.81;   /* 任意の E(n) と回転項 */
	for (int m = 0; m < NMaterial; m++) {
		const material_t *p = &Material[m];
		if (m == PEC) {
			CHECK(C1[m] == 0 && C2[m] == 0 && D1[m] == 0 && D2[m] == 0);
			continue;
		}
		const double hn = 0.37, yy = -0.81;
		const double e1 = C1[m] * en + C2[m] * xx;
		const double h1 = D1[m] * hn + D2[m] * yy;
		if (p->type == 1) {
			const double le = p->esgm * Dt / EPS0;     /* σ dt / ε0 */
			CHECK_NEAR(p->epsr * (e1 - en) + le * e1, xx, 1e-6);
		}
		else {
			/* 分散性 : ε∞ (E1 - E0) + ξ0 E1 = X (D の保存と同じ形) */
			double xi0, dxi0, ke;
			dispersionFactor(p->ae, p->be, p->ce, Dt, &xi0, &dxi0, &ke);
			CHECK_NEAR(p->einf * (e1 - en) + xi0 * e1, xx, 1e-6);
		}
		const double lm = p->msgm * Dt / MU0;          /* σm dt / μ0 */
		CHECK_NEAR(p->amur * (h1 - hn) + lm * h1, yy, 1e-6);
	}
	CHECK(C1[0] == 1 && C2[0] == 1 && D1[0] == 1 && D2[0] == 1);
	/* 導体中の減衰率 ≈ exp(-σ dt / ε) */
	{
		const double r = Material[2].esgm * Dt / (EPS0 * Material[2].epsr);
		CHECK(r > 1e-3 && r < 0.1);
		CHECK_NEAR(C1[2], exp(-r), r * r);
		const double rm = Material[3].msgm * Dt / (MU0 * Material[3].amur);
		CHECK(rm > 1e-3 && rm < 0.1);
		CHECK_NEAR(D1[3], exp(-rm), rm * rm);
	}

	/* メッシュ係数 */
	{
		int ok = 1;
		real_t *rn[3] = {RXn, RYn, RZn}, *rc[3] = {RXc, RYc, RZc};
		for (int a = 0; a < 3; a++) {
			for (int i = 0; i <= n_[a]; i++) ok &= (fabs(rn[a][i] - cdt / dual_width(a, i)) < 1e-6 * rn[a][i]);
			for (int i = 0; i <  n_[a]; i++) ok &= (fabs(rc[a][i] - cdt / (node_[a][i + 1] - node_[a][i])) < 1e-6 * rc[a][i]);
		}
		CHECK(ok);
	}

	/* 入力で与えた dt / Tw は変えない */
	run_setup(0.5e-12, 7e-10);
	CHECK(Dt == 0.5e-12);
	CHECK(Tw == 7e-10);

	/* f0 の決め方 : frequency1 が無ければ frequency2、どちらも無ければ 1/(20 dt) */
	NFreq1 = 0;
	NFreq2 = 2;
	Freq2 = (double *)malloc(NFreq2 * sizeof(double));
	Freq2[0] = 1e9; Freq2[1] = 2e9;
	Solver.maxiter = 0;                 /* setupDft は何もしない */
	run_setup(1e-12, 0);
	CHECK_NEAR(Tw, 1.27 / 1.5e9, 1e-20);
	free(Freq2);
	Freq2 = NULL;
	NFreq2 = 0;
	run_setup(1e-12, 0);
	CHECK_NEAR(Tw, 1.27 * 20 * 1e-12, 1e-24);

	st_teardown();
	Dt = dt0;
	Tw = tw0;
}
