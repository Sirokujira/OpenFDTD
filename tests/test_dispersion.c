/*
test_dispersion.c

分散性材料 (material の type = 2) のテスト。

    eps(w) = einf + ae / (jw) + be / (ce + jw)
    chi(t) = ae + be exp(-ce t)          (t >= 0)

対象:
  dispersionFactor()  sol/setupDispersion.c  RC (再帰的畳み込み) の係数
  setupDispersion()   sol/setupDispersion.c  係数と格子点リストの組み立て
  dispersionE[xyz]()  sol/dispersionE[xyz].c 電界の補正と畳み込み和の更新

期待値の出所:
  - 係数は chi(t) の区間積分を数値積分 (Simpson) で求めて比較する
    (xi0 = ∫[0,dt] chi, dxi0 = ∫[0,dt] - ∫[dt,2dt] の be 項, ke = 隣り合う区間の比)。
  - 時間発展は「回転 (curl H) が 0 の一様媒質」で解析解と比較する。
    このとき dD/dt = 0 なので D は保存し、E は
        E(s) = D0 / (s eps(s))
    の逆 Laplace 変換で与えられる。
      Debye (ae = 0) :
        E(t) = E0 [einf/es + (1 - einf/es) exp(-t/te)],
        es = einf + be/ce, te = einf / (ce es)
        十分時間が経つと E -> E0 einf / es (離散でも D 保存から厳密に成り立つ)
      導電率のみ (be = ce = 0) :
        E(t) = E0 exp(-ae t / einf)
  - 平面波入射では畳み込みに全電界 (散乱界 + 入射界) を使う。入射界を
    各成分の Yee 位置 (Ex は (Xc, Yn, Zn) など) と偏波成分 ei[0..2] で
    finc() から求めて比較する。斜め入射にして 3 座標すべてが効くようにする。

curl H = 0 のとき、電界の通常更新 (updateEx 等) は E <- C1 E になる。
C1 = einf / (einf + xi0) (sol/setup.c) はここでは数値積分した xi0 で作る。
コード側の xi0 が違えば D が保存せず、最終値がずれて検出される。
*/

#include "ofd.h"
#include "ofd_prototype.h"
#include "finc.h"
#include "ofd_test.h"

#include <stdlib.h>
#include <string.h>

#define N      4
#define DELTA  1e-3
#define DT     1e-12

/* chi(t) の be 項 / ae 項を [t0, t1] で Simpson 積分する */
static double chi_integral(double ae, double be, double ce, double t0, double t1)
{
	const int n = 2000;
	const double h = (t1 - t0) / n;
	double s = 0;
	for (int i = 0; i <= n; i++) {
		const double t = t0 + i * h;
		const double f = ae + be * exp(-ce * t);
		const double w = ((i == 0) || (i == n)) ? 1 : ((i % 2) ? 4 : 2);
		s += w * f;
	}
	return s * h / 3;
}

static void test_factor(void)
{
	const double ae = 3e9, be = 4e10, ce = 1e10;
	double xi0, dxi0, ke;

	/* Debye + 導電率 */
	dispersionFactor(ae, be, ce, DT, &xi0, &dxi0, &ke);
	const double i0 = chi_integral(ae, be, ce, 0, DT);
	const double i1b = chi_integral(0, be, ce, DT, 2 * DT);
	const double i0b = chi_integral(0, be, ce, 0, DT);
	const double i2b = chi_integral(0, be, ce, 2 * DT, 3 * DT);
	CHECK_NEAR(xi0, i0, 1e-10 * i0);
	CHECK_NEAR(dxi0, i0b - i1b, 1e-8 * (i0b - i1b));
	CHECK_NEAR(ke, i1b / i0b, 1e-10);
	CHECK_NEAR(ke, i2b / i1b, 1e-10);

	/* 総和 sum_m chi_m (be 項) = be/ce : xi0 の be 分 / (1 - ke) */
	CHECK_NEAR((xi0 - ae * DT) / (1 - ke), be / ce, 1e-9 * be / ce);

	/* ce = 0 (導電率のみ / be を定数とみなす極限) : 0/0 にならない */
	dispersionFactor(ae, 0, 0, DT, &xi0, &dxi0, &ke);
	CHECK(isfinite(xi0) && isfinite(dxi0) && isfinite(ke));
	CHECK_NEAR(xi0, ae * DT, 1e-12 * ae * DT);
	CHECK_NEAR(dxi0, 0, 1e-30);
	CHECK_NEAR(ke, 1, 1e-15);

	dispersionFactor(ae, be, 0, DT, &xi0, &dxi0, &ke);
	CHECK_NEAR(xi0, (ae + be) * DT, 1e-12 * (ae + be) * DT);
	CHECK_NEAR(dxi0, 0, 1e-30);

	/* ce -> 0 の連続性 */
	double xs, ds, ks;
	dispersionFactor(ae, be, 1e-6 / DT, DT, &xs, &ds, &ks);
	CHECK_NEAR(xs, xi0, 1e-5 * xi0);
	CHECK_NEAR(ds, dxi0, 1e-5 * be * DT);
}

/* ---- 格子と材料 ---- */

static void setup_grid(void)
{
	Nx = Ny = Nz = N;
	iMin = jMin = kMin = 0;
	iMax = jMax = kMax = N;
	Dt = DT;

	Xn = (double *)malloc((N + 1) * sizeof(double));
	Yn = (double *)malloc((N + 1) * sizeof(double));
	Zn = (double *)malloc((N + 1) * sizeof(double));
	Xc = (double *)malloc(N * sizeof(double));
	Yc = (double *)malloc(N * sizeof(double));
	Zc = (double *)malloc(N * sizeof(double));
	/* 軸ごとに座標をずらし、x/y/z の取り違えを位相差で見えるようにする */
	for (int i = 0; i <= N; i++) {
		Xn[i] = DELTA * i;
		Yn[i] = DELTA * i + 0.3 * DELTA;
		Zn[i] = DELTA * i - 0.7 * DELTA;
	}
	for (int i = 0; i < N; i++) {
		Xc[i] = (Xn[i] + Xn[i + 1]) / 2;
		Yc[i] = (Yn[i] + Yn[i + 1]) / 2;
		Zc[i] = (Zn[i] + Zn[i + 1]) / 2;
	}

	/* setupSize() と同じ配列レイアウト (tests/test_setupid.c 参照) */
	Nk = 1;
	Nj = (kMax - kMin) + 2 + 1;
	Ni = ((jMax - jMin) + 2 + 1) * Nj;
	N0 = -((iMin - 1) * Ni + (jMin - 1) * Nj + (kMin - 1) * Nk);
	NN = NA(iMax + 1, jMax + 1, kMax + 1) + 1;

	iEx = (id_t *)calloc(NN, sizeof(id_t));
	iEy = (id_t *)calloc(NN, sizeof(id_t));
	iEz = (id_t *)calloc(NN, sizeof(id_t));
	Ex = (real_t *)calloc(NN, sizeof(real_t));
	Ey = (real_t *)calloc(NN, sizeof(real_t));
	Ez = (real_t *)calloc(NN, sizeof(real_t));

	/* 0=真空, 1=PEC, 2=Debye, 3=導電率のみ (ce = 0), 4=誘電体 */
	NMaterial = 5;
	Material = (material_t *)calloc(NMaterial, sizeof(material_t));
	for (int m = 0; m < NMaterial; m++) {
		Material[m].type = 1;
		Material[m].epsr = Material[m].amur = 1;
	}
	Material[2].type = 2;
	Material[2].einf = 2;
	Material[2].ae   = 0;
	Material[2].be   = 4e10;   /* es = 2 + 4 = 6 */
	Material[2].ce   = 1e10;   /* ce dt = 0.01 */
	Material[3].type = 2;
	Material[3].einf = 1.5;
	Material[3].ae   = 1e10;   /* tau = einf / ae = 150 dt */
	Material[3].be   = 0;
	Material[3].ce   = 0;
	Material[4].epsr = 4;

	IPlanewave = 0;
	WaveAmp = 0;
}

static void free_dispersion(void)
{
	free(mDispersionEx); free(mDispersionEy); free(mDispersionEz);
	free(DispersionEx);  free(DispersionEy);  free(DispersionEz);
	mDispersionEx = mDispersionEy = mDispersionEz = NULL;
	DispersionEx = DispersionEy = DispersionEz = NULL;
	numDispersionEx = numDispersionEy = numDispersionEz = 0;
}

static void free_grid(void)
{
	free_dispersion();
	free(Xn); free(Yn); free(Zn); free(Xc); free(Yc); free(Zc);
	free(iEx); free(iEy); free(iEz);
	free(Ex); free(Ey); free(Ez);
	free(Material);
	Xn = Yn = Zn = Xc = Yc = Zc = NULL;
	iEx = iEy = iEz = NULL;
	Ex = Ey = Ez = NULL;
	Material = NULL;
	NMaterial = 0;
	IPlanewave = 0;
}

static void clear_fields(void)
{
	memset(iEx, 0, NN * sizeof(id_t));
	memset(iEy, 0, NN * sizeof(id_t));
	memset(iEz, 0, NN * sizeof(id_t));
	memset(Ex, 0, NN * sizeof(real_t));
	memset(Ey, 0, NN * sizeof(real_t));
	memset(Ez, 0, NN * sizeof(real_t));
}

static int find(const dispersion_t *d, int64_t num, int i, int j, int k)
{
	for (int64_t n = 0; n < num; n++) {
		if ((d[n].i == i) && (d[n].j == j) && (d[n].k == k)) return (int)n;
	}
	return -1;
}

static int same_factor(const dispersion_t *d, id_t m)
{
	return (d->f1 == Material[m].edisp[0])
	    && (d->f2 == Material[m].edisp[1])
	    && (d->f3 == Material[m].edisp[2]);
}

/* 格子点リスト : 各成分の範囲 (自成分方向だけ < Max) と係数の対応 */
static void test_point_list(void)
{
	clear_fields();

	/* 範囲内 */
	IEX(1, 2, 3) = 2;  IEX(0, jMax, kMax) = 3;  IEX(2, 1, 1) = 4;
	IEY(2, 3, 1) = 3;  IEY(iMax, 0, kMax) = 2;
	IEZ(3, 1, 2) = 2;  IEZ(iMax, jMax, 0) = 3;
	/* 自成分方向の上端 (その成分は存在しない) */
	IEX(iMax, 1, 1) = 2;
	IEY(1, jMax, 1) = 2;
	IEZ(1, 1, kMax) = 2;

	setupDispersion();

	CHECK(numDispersionEx == 2);
	CHECK(numDispersionEy == 2);
	CHECK(numDispersionEz == 2);

	int n;
	n = find(mDispersionEx, numDispersionEx, 1, 2, 3);
	CHECK((n >= 0) && same_factor(&mDispersionEx[n], 2));
	n = find(mDispersionEx, numDispersionEx, 0, jMax, kMax);
	CHECK((n >= 0) && same_factor(&mDispersionEx[n], 3));
	n = find(mDispersionEy, numDispersionEy, 2, 3, 1);
	CHECK((n >= 0) && same_factor(&mDispersionEy[n], 3));
	n = find(mDispersionEy, numDispersionEy, iMax, 0, kMax);
	CHECK((n >= 0) && same_factor(&mDispersionEy[n], 2));
	n = find(mDispersionEz, numDispersionEz, 3, 1, 2);
	CHECK((n >= 0) && same_factor(&mDispersionEz[n], 2));
	n = find(mDispersionEz, numDispersionEz, iMax, jMax, 0);
	CHECK((n >= 0) && same_factor(&mDispersionEz[n], 3));

	/* 畳み込み和は 0 から始まる */
	CHECK(DispersionEx[0] == 0 && DispersionEx[1] == 0);
	CHECK(DispersionEy[0] == 0 && DispersionEy[1] == 0);
	CHECK(DispersionEz[0] == 0 && DispersionEz[1] == 0);

	/* 係数そのもの (f1 = 1/(einf + xi0), f2 = dxi0, f3 = ke) */
	const material_t *d = &Material[2];
	const double xi0 = chi_integral(d->ae, d->be, d->ce, 0, DT);
	const double dxi0 = chi_integral(0, d->be, d->ce, 0, DT) - chi_integral(0, d->be, d->ce, DT, 2 * DT);
	CHECK_NEAR(d->edisp[0], 1 / (d->einf + xi0), 1e-6);
	CHECK_NEAR(d->edisp[1], dxi0, 1e-6 * dxi0);
	CHECK_NEAR(d->edisp[2], exp(-d->ce * DT), 1e-6);

	free_dispersion();
}

/*
curl H = 0 の緩和。Ex/Ey/Ez の 1 点ずつに Debye と導電率のみの材料を置く
(成分ごとに位置と初期値を変え、添字・成分の取り違えを検出する)。
*/
typedef struct {
	int comp;        /* 0/1/2 : Ex/Ey/Ez */
	int i, j, k;
	id_t m;
	double e0;
} cell_t;

static real_t *field(int comp, int i, int j, int k)
{
	if (comp == 0) return &EX(i, j, k);
	if (comp == 1) return &EY(i, j, k);
	return &EZ(i, j, k);
}

static void test_relaxation(void)
{
	const cell_t cells[] = {
		{0, 1, 2, 3, 2, 1.0},
		{1, 2, 3, 1, 2, -0.7},
		{2, 3, 1, 2, 2, 0.4},
		{0, 2, 1, 1, 3, 0.9},
		{1, 1, 1, 2, 3, -0.5},
		{2, 2, 2, 1, 3, 0.6},
	};
	const int ncell = (int)(sizeof(cells) / sizeof(cells[0]));
	const real_t sentinel = (real_t)0.123;

	clear_fields();
	for (int c = 0; c < ncell; c++) {
		const cell_t *p = &cells[c];
		id_t *id = (p->comp == 0) ? iEx : ((p->comp == 1) ? iEy : iEz);
		id[NA(p->i, p->j, p->k)] = p->m;
	}
	setupDispersion();

	/* 分散性でない点は触られないこと */
	EX(3, 3, 3) = EY(3, 3, 3) = EZ(3, 3, 3) = sentinel;

	double c1[5];
	for (id_t m = 2; m <= 3; m++) {
		const material_t *d = &Material[m];
		c1[m] = d->einf / (d->einf + chi_integral(d->ae, d->be, d->ce, 0, DT));
	}

	for (int c = 0; c < ncell; c++) {
		*field(cells[c].comp, cells[c].i, cells[c].j, cells[c].k) = (real_t)cells[c].e0;
	}

	const int nstep = 3000;
	double maxerr_debye = 0, maxerr_cond = 0;
	int finite = 1;
	for (int n = 1; n <= nstep; n++) {
		/* 通常更新 (curl H = 0) */
		for (int c = 0; c < ncell; c++) {
			real_t *e = field(cells[c].comp, cells[c].i, cells[c].j, cells[c].k);
			*e = (real_t)(c1[cells[c].m] * *e);
		}
		const double t = n * DT;
		dispersionEx(t);
		dispersionEy(t);
		dispersionEz(t);

		for (int c = 0; c < ncell; c++) {
			const cell_t *p = &cells[c];
			const material_t *d = &Material[p->m];
			const double e = *field(p->comp, p->i, p->j, p->k);
			finite &= isfinite(e);
			if (p->m == 2) {
				const double es = d->einf + d->be / d->ce;
				const double te = d->einf / (d->ce * es);
				const double r = d->einf / es;
				const double ref = p->e0 * (r + (1 - r) * exp(-t / te));
				const double err = fabs(e - ref) / fabs(p->e0);
				if (err > maxerr_debye) maxerr_debye = err;
			}
			else {
				const double ref = p->e0 * exp(-d->ae * t / d->einf);
				const double err = fabs(e - ref) / fabs(p->e0);
				if (err > maxerr_cond) maxerr_cond = err;
			}
		}
	}

	CHECK(finite);
	/* 一次精度の離散化誤差 : Debye の緩和時定数 33 dt, 導電率 150 dt */
	CHECK(maxerr_debye < 1e-2);
	CHECK(maxerr_cond < 5e-3);
	printf("    relaxation max error : Debye %.2e, conductivity %.2e\n", maxerr_debye, maxerr_cond);

	/* 最終値 : Debye は E0 einf/es (D 保存から離散でも厳密), 導電率は 0 へ */
	for (int c = 0; c < ncell; c++) {
		const cell_t *p = &cells[c];
		const material_t *d = &Material[p->m];
		const double e = *field(p->comp, p->i, p->j, p->k);
		if (p->m == 2) {
			const double es = d->einf + d->be / d->ce;
			CHECK_NEAR(e, p->e0 * d->einf / es, 1e-4 * fabs(p->e0));
		}
		else {
			CHECK_NEAR(e, p->e0 * exp(-d->ae * nstep * DT / d->einf), 5e-3 * fabs(p->e0));
		}
	}

	CHECK(EX(3, 3, 3) == sentinel);
	CHECK(EY(3, 3, 3) == sentinel);
	CHECK(EZ(3, 3, 3) == sentinel);

	free_dispersion();
}

/* 平面波 : 畳み込みに入射界 (各成分の Yee 位置・偏波成分) を加えること */
static void test_planewave(void)
{
	clear_fields();
	IEX(1, 2, 3) = 2;
	IEY(2, 3, 1) = 2;
	IEZ(3, 1, 2) = 2;
	setupDispersion();

	/* 斜め入射 (theta = 60, phi = 30 相当) */
	const double th = 60 * DTOR, ph = 30 * DTOR;
	IPlanewave = 1;
	memset(&Planewave, 0, sizeof(Planewave));
	Planewave.ri[0] = -sin(th) * cos(ph);
	Planewave.ri[1] = -sin(th) * sin(ph);
	Planewave.ri[2] = -cos(th);
	Planewave.ei[0] = 0.3;
	Planewave.ei[1] = -0.5;
	Planewave.ei[2] = 0.81;
	/* 基準点を格子の中ほどに置き、各点の遅延をパルス幅程度に収める */
	Planewave.r0[0] = Planewave.r0[1] = Planewave.r0[2] = 2 * DELTA;
	/* パルス幅をセル通過時間 (3.3 ps) 程度にして位置の差を大きく出す */
	Planewave.ai = 1 / 4e-12;

	const double t = 3e-12;
	dispersionEx(t);
	dispersionEy(t);
	dispersionEz(t);

	const double f2 = Material[2].edisp[1];
	real_t fi, dfi;

	finc(Xc[1], Yn[2], Zn[3], t, Planewave.r0, Planewave.ri, 0.3, Planewave.ai, DT, &fi, &dfi);
	CHECK(fabs(fi) > 0.05);
	CHECK_NEAR(DispersionEx[0], f2 * fi, 1e-6 * f2);

	finc(Xn[2], Yc[3], Zn[1], t, Planewave.r0, Planewave.ri, -0.5, Planewave.ai, DT, &fi, &dfi);
	CHECK(fabs(fi) > 0.05);
	CHECK_NEAR(DispersionEy[0], f2 * fi, 1e-6 * f2);

	finc(Xn[3], Yn[1], Zc[2], t, Planewave.r0, Planewave.ri, 0.81, Planewave.ai, DT, &fi, &dfi);
	CHECK(fabs(fi) > 0.05);
	CHECK_NEAR(DispersionEz[0], f2 * fi, 1e-6 * f2);

	/* 散乱界は畳み込み前 0 なので、補正項 f1 * psi も 0 のまま */
	CHECK(EX(1, 2, 3) == 0 && EY(2, 3, 1) == 0 && EZ(3, 1, 2) == 0);

	IPlanewave = 0;
	memset(&Planewave, 0, sizeof(Planewave));
	free_dispersion();
}

void test_dispersion(void)
{
	const double dt0 = Dt;

	test_factor();

	setup_grid();
	test_point_list();
	test_relaxation();
	test_planewave();
	free_grid();

	Dt = dt0;
}
