/*
test_mur.c

Mur 1 次吸収境界 (abc = 0, 既定) のテスト。

対象:
  factorMur()     sol/setupMurHx.c   係数 g = (v dt - d) / (v dt + d)
  setupMurH[xyz]  sol/setupMurH?.c   境界外側の接線 H の点・内側の点・係数
  murH()          sol/murH.c         更新 h_out(n+1) = h_in(n) + g (h_in(n+1) - h_out(n))

期待値の出所:
  - Mur の吸収条件は外向きの 1 次元波動 (∂t - v ∂s) u = 0 (s は内向きの距離)。
    その厳密解 u(s, t) = F(t + s/v) を内側の点 (s = +d/2) に毎ステップ与え、
    外側の点 (s = -d/2) が厳密解に追従することを確かめる。式を写さずに、
    「外へ出ていく波を正しく延長する」という性質そのものを見る。
    差分は (t, s) = ((n+1/2)dt, 0) 中心の箱型差分 (局所打ち切り誤差 2 次)。
    内側には毎ステップ厳密値を与え、外側の誤差は 1 ステップごとに
    |g| < 1 倍されて減衰し蓄積しないので、観測される誤差は 1 ステップ分の
    局所誤差 O(h^3) になる: ω を半分にすると誤差は約 1/8。収束次数で
    判定し、誤差の絶対値の見積りには頼らない (ここでは d/v が dt の 5 倍
    程度あり、空間方向の分解能 ω d / v が誤差を決める)。符号・係数・
    保存値の取り違えは O(1) の誤差か低い次数の収束になる。
  - v dt = d (g = 0) では 1 ステップでちょうど 1 セル進むので、外側の値は
    1 ステップ前の内側の値に厳密に一致する (特性線)。
  - 位相速度 v = c / √(εr μr) の材料は、境界面上の接線 E のうち H に
    直交する成分 (Hx の y 面なら Ez、z 面なら Ey) のもの。面ごとに 1 回ずつ、
    その面の該当 E だけを誘電体 (εr = 4) にして実行し、別の成分・別の面・
    1 つずれた位置を参照すると係数が真空のものになって検出される。
  - d はその面の境界セル幅。不等間隔で、各面で異なる幅にする。
  - PEC に接する境界は g = -1 (外側の値 = 前の内側の値 - (新しい内側 - 前の外側))。
*/

#include "ofd.h"
#include "ofd_prototype.h"
#include "ofd_test.h"

#include <stdlib.h>
#include <string.h>

/* 各軸のセル幅 [m] (不等間隔。両端の境界セルはすべて別の幅) */
static const double DX[] = {1.0e-3, 1.3e-3, 0.8e-3, 1.6e-3};
static const double DY[] = {0.7e-3, 0.9e-3, 1.2e-3, 1.0e-3, 1.25e-3};
static const double DZ[] = {1.5e-3, 1.0e-3, 0.9e-3, 1.1e-3, 1.3e-3, 0.6e-3};

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

static void mur_setup(void)
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
	Dt = 2.0e-12;

	/* setupSize() と同じ配列レイアウト (Mur の余白 1 層) */
	Nk = 1;
	Nj = (kMax - kMin) + 2 + 1;
	Ni = ((jMax - jMin) + 2 + 1) * Nj;
	N0 = -((iMin - 1) * Ni + (jMin - 1) * Nj + (kMin - 1) * Nk);
	NN = NA(iMax + 1, jMax + 1, kMax + 1) + 1;

	iEx = (id_t *)calloc(NN, sizeof(id_t));
	iEy = (id_t *)calloc(NN, sizeof(id_t));
	iEz = (id_t *)calloc(NN, sizeof(id_t));
	Hx = (real_t *)calloc(NN, sizeof(real_t));
	Hy = (real_t *)calloc(NN, sizeof(real_t));
	Hz = (real_t *)calloc(NN, sizeof(real_t));

	/* 0=真空, 1=PEC, 2=誘電体 (εr = 4 → v = c/2), 3=磁性体 (μr = 2.25 → v = c/1.5) */
	NMaterial = 4;
	Material = (material_t *)calloc(NMaterial, sizeof(material_t));
	for (int m = 0; m < NMaterial; m++) {
		Material[m].type = 1;
		Material[m].epsr = Material[m].amur = 1;
	}
	Material[2].epsr = 4;
	Material[3].amur = 2.25;
}

static void mur_teardown(void)
{
	free(Xn); free(Yn); free(Zn); free(Xc); free(Yc); free(Zc);
	free(iEx); free(iEy); free(iEz);
	free(Hx); free(Hy); free(Hz);
	free(Material);
	Xn = Yn = Zn = Xc = Yc = Zc = NULL;
	iEx = iEy = iEz = NULL;
	Hx = Hy = Hz = NULL;
	Material = NULL;
	NMaterial = 0;
	numMurHx = numMurHy = numMurHz = 0;
}

static void test_factor(void)
{
	const double cdt = C * Dt;

	/* g = 0 ⇔ v dt = d : 真空で d = c dt、εr = 4 で d = c dt / 2 */
	CHECK_NEAR(factorMur(cdt, 0), 0, 1e-12);
	CHECK_NEAR(factorMur(cdt / 2, 2), 0, 1e-12);
	CHECK_NEAR(factorMur(cdt / 1.5, 3), 0, 1e-12);

	/* セルが粗いほど (v dt << d) g → -1、PEC は -1 */
	CHECK(factorMur(1e3 * cdt, 0) < -0.99);
	CHECK(factorMur(1e3 * cdt, 0) > -1);
	CHECK(factorMur(cdt, PEC) == -1);

	/* 一般の値 : (v dt - d)/(v dt + d)。v dt = 0.6 d なら -0.25 */
	CHECK_NEAR(factorMur(cdt / 2 / 0.6, 2), -0.25, 1e-12);
}

/* ---- 成分ごとの点リストと更新 ---- */

typedef struct {
	const char *name;
	int a;                  /* H の向き */
	void (*setup)(int);
	mur_t **list;
	int64_t *num;
	real_t **h;
} comp_t;

/* 点 p の軸 b の添字 */
static int coord(const mur_t *p, int b)
{
	return (b == 0) ? p->i : ((b == 1) ? p->j : p->k);
}

static id_t **eid(int c)
{
	return (c == 0) ? &iEx : ((c == 1) ? &iEy : &iEz);
}

static double face_width(int b, int side)
{
	const int n = n_[b];
	const double *x = node_[b];
	return side ? (x[n] - x[n - 1]) : (x[1] - x[0]);
}

/* 成分 c の面 (法線 b, side) での検査。その面の接線 E (向き e = 3 - a - b) を
   材料 m にして setup → 点・係数・更新を確かめる */
static double omega_dt;   /* 外向き波の角周波数 × dt */

static void check_face(const comp_t *c, int b, int side, id_t m, double *maxerr)
{
	const int a = c->a;
	const int e = 3 - a - b;
	const int nb = n_[b];

	/* E の材料 : 境界面 (b 方向の節点 0 か nb) 上の E_e だけを m にする */
	for (int t = 0; t < 3; t++) memset(*eid(t), 0, NN * sizeof(id_t));
	for (int i = 0; i <= Nx; i++) {
	for (int j = 0; j <= Ny; j++) {
	for (int k = 0; k <= Nz; k++) {
		const int idx[3] = {i, j, k};
		if (idx[b] == (side ? nb : 0)) (*eid(e))[NA(i, j, k)] = m;
	}
	}
	}

	c->setup(0);
	const int64_t num = *c->num;
	/* H_a は a 方向に節点 (N_a + 1 個)、他の 2 軸はセル。外側の点は、
	   法線方向に -1 か N、もう一方の軸は内部のセル (辺は含まない) */
	const int b2 = 3 - a - b;
	const int64_t expect = (int64_t)(n_[a] + 1) * 2 * (n_[b] + n_[b2]);
	CHECK(num == expect);
	*c->list = (mur_t *)malloc(num * sizeof(mur_t));
	c->setup(1);

	const double d = face_width(b, side);
	const double v = C / sqrt(Material[m].epsr * Material[m].amur);
	const double g = (m == PEC) ? -1 : (v * Dt - d) / (v * Dt + d);

	/* この面の点 : 法線方向の添字が -1 / nb */
	int nface = 0, ok_inner = 1, ok_g = 1;
	for (int64_t n = 0; n < num; n++) {
		const mur_t *p = &(*c->list)[n];
		const int idx[3] = {p->i, p->j, p->k};
		const int in1[3] = {p->i1, p->j1, p->k1};
		if (idx[b] != (side ? nb : -1)) continue;
		nface++;
		for (int t = 0; t < 3; t++) {
			const int want = (t == b) ? (side ? nb - 1 : 0) : idx[t];
			ok_inner &= (in1[t] == want);
		}
		ok_g &= (fabs(p->g - g) < 1e-6);
	}
	CHECK(nface == (n_[a] + 1) * n_[b2]);
	CHECK(ok_inner);
	CHECK(ok_g);

	/* 外向きの波 u(s, t) = F(t + s/v), F = sin(ω t)。s は内向きの距離。
	   内側 (s = +d/2) に厳密解を与え続け、外側 (s = -d/2) を比べる */
	const double omega = omega_dt / Dt;
	real_t *h = *c->h;
	for (int64_t n = 0; n < num; n++) {
		mur_t *p = &(*c->list)[n];
		h[NA(p->i, p->j, p->k)] = (real_t)sin(omega * (0 - d / 2 / v));
		h[NA(p->i1, p->j1, p->k1)] = (real_t)sin(omega * (0 + d / 2 / v));
		p->f = h[NA(p->i1, p->j1, p->k1)];
	}
	double err = 0;
	for (int step = 1; step <= (int)(20 / omega_dt); step++) {
		const double t = step * Dt;
		for (int64_t n = 0; n < num; n++) {
			const mur_t *p = &(*c->list)[n];
			if (coord(p, b) != (side ? nb : -1)) continue;
			h[NA(p->i1, p->j1, p->k1)] = (real_t)sin(omega * (t + d / 2 / v));
		}
		murH(num, *c->list, h);
		for (int64_t n = 0; n < num; n++) {
			const mur_t *p = &(*c->list)[n];
			if (coord(p, b) != (side ? nb : -1)) continue;
			const double exact = sin(omega * (t - d / 2 / v));
			err = fmax(err, fabs(h[NA(p->i, p->j, p->k)] - exact));
		}
	}
	if (m != PEC) *maxerr = fmax(*maxerr, err);

	free(*c->list);
	*c->list = NULL;
}

/* v dt = d (g = 0) の特性線 : 外側 = 1 ステップ前の内側 (厳密) */
static void test_characteristic(void)
{
	mur_t p;
	real_t h[4] = {0, 0, 0, 0};
	const int64_t save_N0 = N0, save_Ni = Ni, save_Nj = Nj, save_Nk = Nk;
	N0 = 0; Ni = 0; Nj = 0; Nk = 1;      /* NA(0,0,k) = k */
	p.i = p.j = 0; p.k = 0;
	p.i1 = p.j1 = 0; p.k1 = 1;
	p.g = 0;
	p.f = 0;
	int exact = 1;
	for (int n = 1; n <= 20; n++) {
		const real_t prev = h[1];
		h[1] = (real_t)(n * n % 7 - 3);  /* 任意の内側の値 */
		murH(1, &p, h);
		exact &= (h[0] == prev);
	}
	CHECK(exact);
	N0 = save_N0; Ni = save_Ni; Nj = save_Nj; Nk = save_Nk;
}

void test_mur(void)
{
	const double dt0 = Dt;
	const int abc0 = iABC;

	mur_setup();
	test_factor();

	comp_t comps[3] = {
		{"Hx", 0, setupMurHx, &fMurHx, &numMurHx, &Hx},
		{"Hy", 1, setupMurHy, &fMurHy, &numMurHy, &Hy},
		{"Hz", 2, setupMurHz, &fMurHz, &numMurHz, &Hz},
	};

	/* 外向き波を 2 つの角周波数で与える。誤差は 1 ステップ分の局所誤差
	   O(h^3) なので ω を半分にすると約 1/8 (先頭のコメント参照) */
	double err[2] = {0, 0};
	const double odt[2] = {0.04, 0.02};
	for (int r = 0; r < 2; r++) {
		omega_dt = odt[r];
		for (int t = 0; t < 3; t++) {
			const int a = comps[t].a;
			for (int b = 0; b < 3; b++) {
				if (b == a) continue;
				for (int side = 0; side < 2; side++) {
					check_face(&comps[t], b, side, 2, &err[r]);
					check_face(&comps[t], b, side, 3, &err[r]);
				}
			}
			/* PEC に接する面 (g = -1) を 1 つ */
			check_face(&comps[t], (a + 1) % 3, 0, PEC, &err[r]);
		}
	}
	printf("    outgoing wave : max |h_out - exact| = %.2e (omega dt = %.2f), %.2e (%.2f), ratio %.2f\n",
		err[0], odt[0], err[1], odt[1], err[1] / err[0]);
	CHECK(err[0] < 0.01);
	CHECK((err[1] / err[0] > 0.09) && (err[1] / err[0] < 0.18));

	test_characteristic();

	mur_teardown();
	Dt = dt0;
	iABC = abc0;
}
