/*
test_feedload.c

給電 (feed) と集中定数負荷 (load = R / C / L) のテスト。

対象:
  setup_feed() / setup_load()   sol/input2.c   格子点と辺の長さ (getindex)、
                                               R/C の材料化、L の係数
  efeed()                       sol/efeed.c    内部抵抗つき電圧源
  eload()                       sol/eload.c    集中インダクタ

期待値の出所 (いずれも回路・電磁気の関係から):
  - 辺を流れる電流は、辺を囲む H の周回積分 I = ∮H·dl (Ampère)。
    テストでは辺の向き â に沿った線電流の作る磁界の向き (右ねじ, â × r̂)
    で辺の周りの 4 つの H を与え、周回の各辺の長さ (双対セル幅) を掛けた
    和から I を求める。実装の差分の書き方は写さない。H は η0 倍で格納
    されている (E と同じ単位) ので I = ∮H·dl / η0。
  - 給電は内部抵抗 R の Thevenin 電源: V = Vs - R I、辺の電界 E = -V / d。
    Vs = volt × vfeed(t)、E は整数ステップの時刻 t = (n+1) dt。
    パルスが終わった後 (|vfeed| <= 1e-6) は、PEC の辺 (線上の給電) だけ
    電源として保持し、そうでない辺は通常の更新に任せる (書き換えない)。
  - 双対セル幅: 内部の節点は隣り合うセル中心の距離、境界の節点は境界セル幅
    (setup_mesh の RXn[0] / RXn[Nx] と同じ。修正前はここで配列外を読んでいた)。
  - R / C 負荷は、長さ d・断面積 A の辺を満たす材料として
        R = d / (σ A)  → σ = (d/A) / R,     C = ε A / d  → εr = C (d/A) / ε0
  - L 負荷: 辺の容量 C_edge = ε0 A / d と並列の LC 共振器になるので、
    H を 0 に固定すると E は ω0 = 1 / √(L C_edge) で振動する (回路の固有値)。
    一定の電流 I0 (一定の H) で駆動すると、E は
        E(t) = (c ∮H·dl / A) / ω0 · sin(ω0 t)   (振幅 = I0 / (ω0 C_edge) / d)
    で、1 ステップ目は E = c dt · (∮H·dl / A) に厳密に一致する。
*/

#include "ofd.h"
#include "ofd_prototype.h"
#include "ofd_test.h"

#include <stdint.h>
#include <stdlib.h>
#include <string.h>

static const double DX[] = {1.0e-3, 1.3e-3, 0.8e-3, 1.6e-3};
static const double DY[] = {0.7e-3, 0.9e-3, 1.2e-3, 1.0e-3, 1.25e-3};
static const double DZ[] = {1.5e-3, 1.0e-3, 0.9e-3, 1.1e-3, 1.3e-3, 0.6e-3};

static int     n_[3];
static double *node_[3], *cell_[3];
static real_t **efield[3], **hfield[3];
static id_t   **eids[3];

static void make_axis(const double *d, int n, double **nodes, double **cells)
{
	/* 配列の前後に番兵を置かない (配列外読みはここで検出したい値になる) */
	double *xn = (double *)malloc((n + 1) * sizeof(double));
	double *xc = (double *)malloc(n * sizeof(double));
	xn[0] = -0.4e-3;
	for (int i = 0; i < n; i++) {
		xn[i + 1] = xn[i] + d[i];
		xc[i] = (xn[i] + xn[i + 1]) / 2;
	}
	*nodes = xn;
	*cells = xc;
}

static uint32_t rng = 777;
static double urand(void)
{
	rng = rng * 1664525u + 1013904223u;
	return ((rng >> 8) & 0xffff) / 65536.0 - 0.5;
}

static void fl_setup(void)
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

	Dt = 1.5e-12;
	Tw = 40 * Dt;

	Nk = 1;
	Nj = (kMax - kMin) + 2 + 1;
	Ni = ((jMax - jMin) + 2 + 1) * Nj;
	N0 = -((iMin - 1) * Ni + (jMin - 1) * Nj + (kMin - 1) * Nk);
	NN = NA(iMax + 1, jMax + 1, kMax + 1) + 1;

	Ex = (real_t *)calloc(NN, sizeof(real_t));
	Ey = (real_t *)calloc(NN, sizeof(real_t));
	Ez = (real_t *)calloc(NN, sizeof(real_t));
	Hx = (real_t *)calloc(NN, sizeof(real_t));
	Hy = (real_t *)calloc(NN, sizeof(real_t));
	Hz = (real_t *)calloc(NN, sizeof(real_t));
	iEx = (id_t *)calloc(NN, sizeof(id_t));
	iEy = (id_t *)calloc(NN, sizeof(id_t));
	iEz = (id_t *)calloc(NN, sizeof(id_t));
	efield[0] = &Ex; efield[1] = &Ey; efield[2] = &Ez;
	hfield[0] = &Hx; hfield[1] = &Hy; hfield[2] = &Hz;
	eids[0] = &iEx; eids[1] = &iEy; eids[2] = &iEz;
}

static void fl_teardown(void)
{
	free(Xn); free(Yn); free(Zn); free(Xc); free(Yc); free(Zc);
	free(Ex); free(Ey); free(Ez); free(Hx); free(Hy); free(Hz);
	free(iEx); free(iEy); free(iEz);
	Xn = Yn = Zn = Xc = Yc = Zc = NULL;
	Ex = Ey = Ez = Hx = Hy = Hz = NULL;
	iEx = iEy = iEz = NULL;
}

/* 節点 i の双対セル幅 (境界は境界セル幅) */
static double dual_width(int a, int i)
{
	const int n = n_[a];
	const double *x = node_[a];
	if (i <= 0) return x[1] - x[0];
	if (i >= n) return x[n] - x[n - 1];
	return (x[i + 1] - x[i - 1]) / 2;
}

/* ---- 辺を囲む H : 線電流の右ねじ ---- */

/* 辺 (向き a, 添字 idx) を囲む 4 つの H を「â に沿って流れる電流の磁界」の向きに
   大きさ h で置き、∮H·dl を返す。周りの H はあらかじめ乱数で埋めておく
   (別の位置を読むと値が変わって検出される)。 */
static double set_loop_h(int a, const int idx[3], double h)
{
	double loop = 0;
	for (int s = 1; s <= 2; s++) {
		const int b = (a + s) % 3;           /* 辺に垂直な軸 */
		const int c = 3 - a - b;             /* そこに置く H の成分 */
		for (int side = -1; side <= 1; side += 2) {
			/* H_c のセル中心は b 方向に +1/2 (添字 idx_b) と -1/2 (idx_b - 1) */
			int q[3] = {idx[0], idx[1], idx[2]};
			if (side < 0) q[b] -= 1;
			/* 向き: (â × side b̂) の c 成分。â × b̂ は (b == a+1) で +ĉ、それ以外 -ĉ */
			const double sgn = side * ((b == (a + 1) % 3) ? +1.0 : -1.0);
			(*hfield[c])[NA(q[0], q[1], q[2])] = (real_t)(sgn * h);
			/* 周回の向き (â まわりに右回り) と H が平行なので H·dl = h × 長さ。
			   H_c の区間は c 方向に並ぶので長さは c 方向の双対セル幅 */
			loop += h * dual_width(c, idx[c]);
		}
	}
	return loop;
}

static void fill_random_h(void)
{
	for (int64_t q = 0; q < NN; q++) {
		Hx[q] = (real_t)(10 * urand());
		Hy[q] = (real_t)(10 * urand());
		Hz[q] = (real_t)(10 * urand());
	}
}

/* ---- setup_feed : 格子点と辺の長さ ---- */

static void test_getindex(void)
{
	/* X 給電は x がセル中心、y,z が節点。Y, Z も同様 */
	const char dirs[5] = {'X', 'Y', 'Z', 'Z', 'Z'};
	const double fx[5] = {Xc[1], Xn[2], Xn[3], Xn[0], Xn[Nx]};
	const double fy[5] = {Yn[3], Yc[2], Yn[1], Yn[2], Yn[Ny]};
	const double fz[5] = {Zn[4], Zn[2], Zc[3], Zc[0], Zc[Nz - 1]};
	const int want[5][3] = {{1, 3, 4}, {2, 2, 2}, {3, 1, 3}, {0, 2, 0}, {Nx, Ny, Nz - 1}};

	NFeed = 5;
	Feed = (feed_t *)calloc(NFeed, sizeof(feed_t));
	for (int n = 0; n < NFeed; n++) Feed[n].dir = dirs[n];
	setup_feed(fx, fy, fz);

	for (int n = 0; n < NFeed; n++) {
		const feed_t *f = &Feed[n];
		const int a = f->dir - 'X';
		const int idx[3] = {f->i, f->j, f->k};
		const double len[3] = {f->dx, f->dy, f->dz};
		CHECK(idx[0] == want[n][0] && idx[1] == want[n][1] && idx[2] == want[n][2]);
		for (int b = 0; b < 3; b++) {
			const double expect = (b == a)
				? (node_[b][idx[b] + 1] - node_[b][idx[b]])   /* 辺そのもの : セル幅 */
				: dual_width(b, idx[b]);                       /* 周回の辺 : 双対幅 */
			CHECK_NEAR(len[b], expect, 1e-15);
		}
	}

	free(Feed);
	Feed = NULL;
	NFeed = 0;
}

/* ---- efeed ---- */

static void test_efeed(void)
{
	const double R = 50;
	const double h = 0.37;

	NFeed = 3;
	Feed = (feed_t *)calloc(NFeed, sizeof(feed_t));
	const char dirs[3] = {'X', 'Y', 'Z'};
	const double fx[3] = {Xc[1], Xn[2], Xn[3]};
	const double fy[3] = {Yn[3], Yc[2], Yn[1]};
	const double fz[3] = {Zn[4], Zn[2], Zc[3]};
	for (int n = 0; n < NFeed; n++) Feed[n].dir = dirs[n];
	setup_feed(fx, fy, fz);
	Feed[0].volt = 1.0;  Feed[0].delay = 0;
	Feed[1].volt = -2.5; Feed[1].delay = 7 * Dt;
	Feed[2].volt = 0.8;  Feed[2].delay = 3 * Dt;
	rFeed = R;

	Solver.maxiter = 400;
	VFeed = (double *)calloc(NFeed * (Solver.maxiter + 1), sizeof(double));
	IFeed = (double *)calloc(NFeed * (Solver.maxiter + 1), sizeof(double));

	/* 全給電が空隙 (pec = 0) / 全給電が線上 (pec = 1) の両方で、
	   パルス中とパルス後を調べる (3 方向とも両方の分岐を踏む) */
	const int times[2] = {30, 390};   /* パルス中 / パルス後 */
	for (int run = 0; run < 4; run++) {
		const int pec = run / 2;
		const int tt = run % 2;
		const int itime = times[tt];
		for (int n = 0; n < NFeed; n++) {
			(*eids[Feed[n].dir - 'X'])[NA(Feed[n].i, Feed[n].j, Feed[n].k)] = (id_t)(pec ? PEC : 0);
		}
		fill_random_h();
		double cur[3];
		real_t before[3];
		for (int n = 0; n < NFeed; n++) {
			const feed_t *f = &Feed[n];
			const int idx[3] = {f->i, f->j, f->k};
			cur[n] = set_loop_h(f->dir - 'X', idx, h * (n + 1)) / ETA0;
			real_t *e = &(*efield[f->dir - 'X'])[NA(f->i, f->j, f->k)];
			*e = (real_t)(5 + n);         /* 通常更新で入っていた値 */
			before[n] = *e;
		}

		efeed(itime);

		for (int n = 0; n < NFeed; n++) {
			const feed_t *f = &Feed[n];
			const int a = f->dir - 'X';
			const double d = (a == 0) ? f->dx : ((a == 1) ? f->dy : f->dz);
			const double t = (itime + 1) * Dt;
			const double v0 = vfeed(t, Tw, f->delay);
			const double vs = f->volt * v0;
			const double v = vs - R * cur[n];
			const int64_t id = n * (Solver.maxiter + 1) + itime;
			const real_t e = (*efield[a])[NA(f->i, f->j, f->k)];

			CHECK_NEAR(IFeed[id], cur[n], 1e-6 * fabs(cur[n]));   /* H は float で格納 */
			CHECK_NEAR(VFeed[id], v, 1e-6 * fmax(fabs(v), 1));
			if ((fabs(v0) > 1e-6) || pec) {
				CHECK_NEAR(e, -v / d, 1e-5 * fabs(v / d));
			}
			else {
				CHECK(e == before[n]);     /* パルス後の空隙は書き換えない */
			}
		}
		if (tt == 0) {
			CHECK(fabs(vfeed((itime + 1) * Dt, Tw, 0)) > 0.1);   /* 本当にパルス中 */
		}
		else {
			CHECK(fabs(vfeed((itime + 1) * Dt, Tw, 7 * Dt)) < 1e-6);
		}
	}

	free(VFeed); free(IFeed);
	VFeed = IFeed = NULL;
	free(Feed);
	Feed = NULL;
	NFeed = 0;
	rFeed = 0;
	memset(iEx, 0, NN * sizeof(id_t));
	memset(iEy, 0, NN * sizeof(id_t));
	memset(iEz, 0, NN * sizeof(id_t));
}

/* ---- setup_load : R / C は材料、L は係数 ---- */

static void test_setup_load(void)
{
	const int inc = 4;
	NMaterial = 2;
	Material = (material_t *)calloc(NMaterial, sizeof(material_t));
	NGeometry = 0;
	Geometry = NULL;
	NInductor = 0;
	Inductor = NULL;

	char   dl[3] = {'X', 'Y', 'Z'};
	char   cl[3] = {'R', 'C', 'L'};
	double xl[3] = {Xc[2], Xn[1], Xn[2]};
	double yl[3] = {Yn[2], Yc[3], Yn[3]};
	double zl[3] = {Zn[3], Zn[4], Zc[1]};
	double pl[3] = {75, 2e-12, 3e-9};
	setup_load(3, dl, xl, yl, zl, cl, pl, inc);

	CHECK(NMaterial == 4);
	CHECK(NGeometry == 2);
	CHECK(NInductor == 1);

	/* R : 辺 x (長さ dx = セル幅, 断面 A = 双対 dy × dz) を満たす導体 */
	{
		const double d = Xn[3] - Xn[2];
		const double A = dual_width(1, 2) * dual_width(2, 3);
		CHECK_NEAR(Material[2].esgm, (d / A) / 75, 1e-9 * (d / A) / 75);
		CHECK_NEAR(Material[2].epsr, 1, 1e-15);
		CHECK(Geometry[0].m == 2);
		CHECK_NEAR(Geometry[0].g[0], Xn[2], 1e-15);
		CHECK_NEAR(Geometry[0].g[1], Xn[3], 1e-15);
		CHECK_NEAR(Geometry[0].g[2], Yn[2], 1e-15);
		CHECK_NEAR(Geometry[0].g[4], Zn[3], 1e-15);
	}
	/* C : εr = C (d/A) / ε0 */
	{
		const double d = Yn[4] - Yn[3];
		const double A = dual_width(2, 4) * dual_width(0, 1);
		const double epsr = 2e-12 * (d / A) / EPS0;
		CHECK_NEAR(Material[3].epsr, epsr, 1e-9 * epsr);
		CHECK_NEAR(Material[3].esgm, 0, 1e-30);
		CHECK(Geometry[1].m == 3);
		CHECK_NEAR(Geometry[1].g[2], Yn[3], 1e-15);
		CHECK_NEAR(Geometry[1].g[3], Yn[4], 1e-15);
	}
	/* L : 位置と辺の長さ */
	{
		const inductor_t *p = &Inductor[0];
		CHECK(p->dir == 'Z' && p->i == 2 && p->j == 3 && p->k == 1);
		CHECK_NEAR(p->dz, Zn[2] - Zn[1], 1e-15);
		CHECK_NEAR(p->dx, dual_width(0, 2), 1e-15);
		CHECK_NEAR(p->dy, dual_width(1, 3), 1e-15);
		CHECK(p->e == 0 && p->esum == 0);
	}

	free(Material); Material = NULL; NMaterial = 0;
	free(Geometry); Geometry = NULL; NGeometry = 0;
	/* Inductor は次のテストで使う */
}

/* ---- eload : LC 共振器 ---- */

static void test_eload(void)
{
	/* 3 方向それぞれに 1 つずつ。L = 100 nH でこのセル (C_edge ≈ 5〜20 fF)
	   なら ω0 dt ≈ 0.03〜0.07 で、離散化による周波数の誤差 (θ²/24) は 2e-4 以下 */
	const char dirs[3] = {'X', 'Y', 'Z'};
	const double fx[3] = {Xc[2], Xn[1], Xn[2]};
	const double fy[3] = {Yn[2], Yc[3], Yn[3]};
	const double fz[3] = {Zn[3], Zn[4], Zc[1]};
	const double L = 100e-9;

	free(Inductor);
	Inductor = NULL;
	NInductor = 0;
	NMaterial = 2;
	Material = (material_t *)calloc(NMaterial, sizeof(material_t));

	for (int n = 0; n < 3; n++) {
		char dl = dirs[n], cl = 'L';
		double xl = fx[n], yl = fy[n], zl = fz[n], pl = L;
		setup_load(1, &dl, &xl, &yl, &zl, &cl, &pl, 4);
	}
	CHECK(NInductor == 3);

	for (int n = 0; n < 3; n++) {
		inductor_t *p = &Inductor[n];
		const int a = p->dir - 'X';
		const int idx[3] = {p->i, p->j, p->k};
		const double d = (a == 0) ? p->dx : ((a == 1) ? p->dy : p->dz);
		const double A = p->dx * p->dy * p->dz / d;
		const double cedge = EPS0 * A / d;
		const double w0 = 1 / sqrt(L * cedge);
		CHECK((w0 * Dt > 0.02) && (w0 * Dt < 0.1));

		/* (1) H = 0 で E0 から自由振動 : 周期 2π/ω0、振幅一定 */
		fill_random_h();
		set_loop_h(a, idx, 0);
		p->e = 1.0;
		p->esum = 0;
		const int nstep = 4000;
		double last = 1.0, tcross0 = -1, tcross = -1, amax = 0;
		int ncross = 0;
		for (int s = 1; s <= nstep; s++) {
			eload();
			const double e = (*efield[a])[NA(p->i, p->j, p->k)];
			if ((last > 0) && (e <= 0)) {
				/* 下向きの零交差 (線形補間) */
				const double tc = (s - 1 + last / (last - e)) * Dt;
				if (ncross == 0) tcross0 = tc;
				tcross = tc;
				ncross++;
			}
			if (s > nstep - 200) amax = fmax(amax, fabs(e));
			last = e;
		}
		const double period = (tcross - tcross0) / (ncross - 1);
		CHECK(ncross > 20);
		CHECK_NEAR(period, 2 * PI / w0, 1e-3 * 2 * PI / w0);
		CHECK_NEAR(amax, 1.0, 0.05);

		/* (2) 一定電流で駆動 : 1 ステップ目は E = c dt ∮H·dl / A (厳密)、
		   振幅は (c ∮H·dl / A) / ω0 */
		fill_random_h();
		const double loop = set_loop_h(a, idx, 0.2);
		const double curl = loop / A;
		p->e = 0;
		p->esum = 0;
		eload();
		const double e1 = (*efield[a])[NA(p->i, p->j, p->k)];
		CHECK_NEAR(e1, C * Dt * curl, 1e-5 * fabs(C * Dt * curl));
		double emax = fabs(e1), esum = e1;
		const int nper = (int)(2 * PI / w0 / Dt + 0.5);
		for (int s = 2; s <= 5 * nper; s++) {
			eload();
			const double e = (*efield[a])[NA(p->i, p->j, p->k)];
			emax = fmax(emax, fabs(e));
			esum += e;
		}
		const double amp = C * curl / w0;
		CHECK_NEAR(emax, amp, 0.04 * amp);
		/* E は 0 を中心に振れる (インダクタ電流が駆動電流の周りで振動し、
		   E の平均は 0)。整数周期の平均がほぼ 0 */
		CHECK(fabs(esum / (5 * nper)) < 0.05 * amp);
	}

	free(Inductor); Inductor = NULL; NInductor = 0;
	free(Material); Material = NULL; NMaterial = 0;
	free(Geometry); Geometry = NULL; NGeometry = 0;
}

void test_feedload(void)
{
	const double dt0 = Dt, tw0 = Tw;

	fl_setup();
	test_getindex();
	test_efeed();
	test_setup_load();
	test_eload();
	fl_teardown();

	Dt = dt0;
	Tw = tw0;
}
