/*
test_observe.c

観測と周波数領域変換のテスト。

対象:
  setupDft()   sol/setupDft.c   波源スペクトルによる正規化 cFdft と DFT 係数
  dftNear3d()  sol/dftNear3d.c  近傍界の DFT の累積
  vpoint()     sol/vpoint.c     観測点の電圧波形
  average()    sol/average.c    収束判定に使う平均の |E|, |H|

期待値の出所:
  - cFdft は波源波形の離散フーリエ和 Σ f(n dt) e^{-jω n dt}。波源は
    ガウス微分 f = A x e^{-x²} (x = (t - t0)/τ) なので連続フーリエ変換が
    閉じた形で書ける:
        F(ω) = A τ e^{-jω t0} · (-j √π ωτ / 2) e^{-(ωτ)²/4}
    パルスが dt に比べて滑らかで、時間窓に収まっていれば和 ≈ F / dt。
    給電 (vfeed: A = √2 e^{1/2}, τ = tw/4, t0 = tw + td) と
    平面波 (finc: A = 1, τ = 1/ai, t0 = 領域中心への到達時刻) の両方。
  - DFT 係数: H は半整数時刻 (n + 1/2) dt、E は整数時刻 (n + 1) dt
    (sol/solve.c の時間進行と同じ) の e^{-jωt} を cFdft で割ったもの。
  - 近傍界の DFT は Σ_n field(n) × 係数(n) の累積 (実部・虚部、周波数ごとの
    オフセット)。範囲は E が自身の格子点、H は補間用のゴースト 1 層を含む。
  - 観測点の電圧は V = -∫E·dl = -E d (辺に沿った電界と辺の長さ)。平面波では
    全電界 = 散乱界 + 入射界 (入射界は E の格子位置・時刻 (n+1) dt)。
  - 平均: 各セルで辺の 4 本 (H は面の 2 枚) の平均をとってから絶対値を
    足し、セル数で割る。一様な場 (e1,e2,e3) なら |e1|+|e2|+|e3|、
    市松模様 (隣の辺で符号反転) なら 0、一次関数ならセル中心の値
    (3 成分とも 3 軸すべてに傾きを持たせる。傾きの無い軸では、その軸方向に
    ずれた辺を読んでも値が変わらず検出できない)。
  - setupDft の時刻を 1 ステップずらす変更 (t = (n+1) dt) は、標本と位相が
    同時にずれるだけで和が変わらない (端の項は 0) 等価な書き換えなので
    検出対象にしない。
*/

#include "ofd.h"
#include "complex.h"
#include "ofd_prototype.h"
#include "finc.h"
#include "ofd_test.h"

#include <stdlib.h>
#include <string.h>

#define N 4

static void ob_setup(void)
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
		Xn[i] = 1e-3 * i;
		Yn[i] = 1.2e-3 * i - 0.5e-3;
		Zn[i] = 0.9e-3 * i + 0.2e-3;
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
	Hx = (real_t *)calloc(NN, sizeof(real_t));
	Hy = (real_t *)calloc(NN, sizeof(real_t));
	Hz = (real_t *)calloc(NN, sizeof(real_t));
	Dt = 1e-12;
	Tw = 60 * Dt;
	NFeed = 0;
	IPlanewave = 0;
	WaveAmp = 0;
	memset(&Planewave, 0, sizeof(Planewave));
}

static void ob_teardown(void)
{
	free(Xn); free(Yn); free(Zn); free(Xc); free(Yc); free(Zc);
	free(Ex); free(Ey); free(Ez); free(Hx); free(Hy); free(Hz);
	Xn = Yn = Zn = Xc = Yc = Zc = NULL;
	Ex = Ey = Ez = Hx = Hy = Hz = NULL;
	NFeed = 0;
	IPlanewave = 0;
	memset(&Planewave, 0, sizeof(Planewave));
}

/* ガウス微分 A x e^{-x²}, x = (t - t0)/τ のフーリエ変換 */
static d_complex_t ft_gauss_deriv(double A, double tau, double t0, double omega)
{
	const double k = omega * tau;
	const double mag = A * tau * sqrt(PI) * k / 2 * exp(-k * k / 4);
	/* -j e^{-jω t0} */
	const double ph = -omega * t0 - PI / 2;
	return d_complex(mag * cos(ph), mag * sin(ph));
}

static void test_setupdft(void)
{
	const double f[3] = {0.5e9, 2e9, 6e9};   /* スペクトルの山 (≈ 1/(π τ)) の前後 */
	NFreq2 = 3;
	Freq2 = (double *)malloc(NFreq2 * sizeof(double));
	memcpy(Freq2, f, sizeof(f));
	Solver.maxiter = 600;
	cFdft = (d_complex_t *)malloc(NFreq2 * sizeof(d_complex_t));
	cEdft = (d_complex_t *)malloc((size_t)NFreq2 * (Solver.maxiter + 1) * sizeof(d_complex_t));
	cHdft = (d_complex_t *)malloc((size_t)NFreq2 * (Solver.maxiter + 1) * sizeof(d_complex_t));

	for (int src = 0; src < 2; src++) {
		double A, tau, t0;
		if (src == 0) {
			/* 給電 : vfeed(t, Tw, 0) */
			NFeed = 1;
			IPlanewave = 0;
			A = sqrt(2.0) * exp(0.5);
			tau = Tw / 4;
			t0 = Tw;
		}
		else {
			/* 平面波 : 領域中心での入射波形 (ei = 1) */
			NFeed = 0;
			IPlanewave = 1;
			Planewave.ri[0] = 0.6; Planewave.ri[1] = 0; Planewave.ri[2] = -0.8;
			Planewave.r0[0] = -2e-3; Planewave.r0[1] = 0; Planewave.r0[2] = 4e-3;
			Planewave.ai = 1 / (12 * Dt);
			A = 1;
			tau = 1 / Planewave.ai;
			const double xc = (Xn[0] + Xn[N]) / 2, yc = (Yn[0] + Yn[N]) / 2, zc = (Zn[0] + Zn[N]) / 2;
			t0 = ((xc - Planewave.r0[0]) * Planewave.ri[0]
			    + (yc - Planewave.r0[1]) * Planewave.ri[1]
			    + (zc - Planewave.r0[2]) * Planewave.ri[2]) / C;
		}
		if (src == 1) {
			/* 平面波は t = 0 で立ち上がり途中だと窓の打ち切り誤差が出るので、
			   基準点を後ろへずらして到達を 150 dt 遅らせる */
			const double shift = 150 * Dt * C;
			for (int a = 0; a < 3; a++) Planewave.r0[a] -= shift * Planewave.ri[a];
			t0 += 150 * Dt;
		}

		setupDft();

		for (int fi = 0; fi < NFreq2; fi++) {
			const double omega = 2 * PI * Freq2[fi];
			const d_complex_t F = ft_gauss_deriv(A, tau, t0, omega);
			const d_complex_t S = d_rmul(Dt, cFdft[fi]);    /* Σ dt ≈ ∫ */
			CHECK(d_abs(d_sub(S, F)) < 1e-4 * d_abs(F));
			if (d_abs(d_sub(S, F)) >= 1e-4 * d_abs(F)) {
				printf("    src %d f=%.2e : |S-F|/|F| = %.2e\n", src, Freq2[fi], d_abs(d_sub(S, F)) / d_abs(F));
			}

			/* DFT 係数 : H は (n+1/2)dt、E は (n+1)dt の e^{-jωt} / cFdft */
			int ok = 1;
			for (int n = 0; n <= Solver.maxiter; n += 37) {
				const int id = n * NFreq2 + fi;
				const d_complex_t eh = d_exp(-omega * (n + 0.5) * Dt);
				const d_complex_t ee = d_exp(-omega * (n + 1.0) * Dt);
				ok &= (d_abs(d_sub(d_mul(cHdft[id], cFdft[fi]), eh)) < 1e-12);
				ok &= (d_abs(d_sub(d_mul(cEdft[id], cFdft[fi]), ee)) < 1e-12);
			}
			CHECK(ok);
		}
	}

	/* dftNear3d : Σ field × 係数。2 ステップ分を累積して全点を確かめる */
	const int64_t nall = (int64_t)NFreq2 * NN;
	float **cr[6] = {&cEx_r, &cEy_r, &cEz_r, &cHx_r, &cHy_r, &cHz_r};
	float **ci[6] = {&cEx_i, &cEy_i, &cEz_i, &cHx_i, &cHy_i, &cHz_i};
	for (int c = 0; c < 6; c++) {
		*cr[c] = (float *)calloc(nall, sizeof(float));
		*ci[c] = (float *)calloc(nall, sizeof(float));
	}
	real_t *fld[6] = {Ex, Ey, Ez, Hx, Hy, Hz};
	const int steps[2] = {5, 123};
	for (int s = 0; s < 2; s++) {
		for (int c = 0; c < 6; c++) {
			for (int64_t q = 0; q < NN; q++) fld[c][q] = (real_t)(0.001 * ((q * 7 + c * 13 + s * 5) % 97) - 0.05);
		}
		dftNear3d(steps[s]);
	}
	for (int c = 0; c < 6; c++) {
		const int isE = (c < 3);
		const int a = c % 3;
		int bad = 0;
		for (int fq = 0; fq < NFreq2; fq++) {
		for (int i = -1; i <= N + 1; i++) {
		for (int j = -1; j <= N + 1; j++) {
		for (int k = -1; k <= N + 1; k++) {
			const int idx[3] = {i, j, k};
			const int64_t q = NA(i, j, k);
			if ((q < 0) || (q >= NN)) continue;
			/* 範囲 : E_a は a 方向セル [0,N)、他は節点 [0,N]。
			   H_a は a 方向節点 [-1,N+1]、他はセル [-1,N] (ゴースト 1 層込み) */
			int in = 1;
			for (int b = 0; b < 3; b++) {
				int lo, hi;
				if (isE) { lo = 0;  hi = (b == a) ? N - 1 : N; }
				else     { lo = -1; hi = (b == a) ? N + 1 : N; }
				in &= (idx[b] >= lo) && (idx[b] <= hi);
			}
			double er = 0, ei = 0;
			if (in) {
				for (int s = 0; s < 2; s++) {
					const double v = (real_t)(0.001 * ((q * 7 + c * 13 + s * 5) % 97) - 0.05);
					const d_complex_t w = isE ? cEdft[steps[s] * NFreq2 + fq] : cHdft[steps[s] * NFreq2 + fq];
					er += v * w.r;
					ei += v * w.i;
				}
			}
			const double scale = 1e-5 * (fabs(er) + fabs(ei)) + 1e-30;
			bad += (fabs((*cr[c])[fq * NN + q] - er) > scale + 1e-6 * d_abs(cEdft[fq]));
			bad += (fabs((*ci[c])[fq * NN + q] - ei) > scale + 1e-6 * d_abs(cEdft[fq]));
		}
		}
		}
		}
		CHECK(bad == 0);
		if (bad) printf("    dftNear3d component %d : %d mismatches\n", c, bad);
	}
	for (int c = 0; c < 6; c++) {
		free(*cr[c]); free(*ci[c]);
		*cr[c] = *ci[c] = NULL;
	}

	free(cFdft); free(cEdft); free(cHdft); free(Freq2);
	cFdft = cEdft = cHdft = NULL;
	Freq2 = NULL;
	NFreq2 = 0;
	NFeed = 0;
	IPlanewave = 0;
	memset(&Planewave, 0, sizeof(Planewave));
}

static void test_vpoint(void)
{
	Solver.maxiter = 50;
	NPoint = 1;
	/* 観測点 1 つ + 差分用の 2 点 (setup_point が足す) = 3 */
	Point = (point_t *)calloc(NPoint + 2, sizeof(point_t));
	const char dirs[3] = {'X', 'Y', 'Z'};
	const int pos[3][3] = {{1, 2, 3}, {2, 1, 3}, {3, 2, 1}};
	for (int n = 0; n < 3; n++) {
		Point[n].dir = dirs[n];
		Point[n].i = pos[n][0]; Point[n].j = pos[n][1]; Point[n].k = pos[n][2];
		Point[n].dx = 1.1e-3 * (n + 1); Point[n].dy = 0.7e-3 * (n + 1); Point[n].dz = 1.3e-3 * (n + 1);
	}
	VPoint = (double *)calloc((NPoint + 2) * (Solver.maxiter + 1), sizeof(double));
	EX(1, 2, 3) = 0.25f; EY(2, 1, 3) = -0.5f; EZ(3, 2, 1) = 0.75f;

	for (int pw = 0; pw < 2; pw++) {
		IPlanewave = pw;
		if (pw) {
			Planewave.ri[0] = -0.48; Planewave.ri[1] = -0.6; Planewave.ri[2] = -0.64;
			Planewave.ei[0] = 0.8;   Planewave.ei[1] = -0.36; Planewave.ei[2] = -0.48;
			Planewave.r0[0] = Planewave.r0[1] = Planewave.r0[2] = 2e-3;
			Planewave.ai = 1 / (5 * Dt);
		}
		const int itime = 3;   /* 入射パルスが観測点付近の山にある時刻 (at ≈ 0.8) */
		vpoint(itime);
		for (int n = 0; n < 3; n++) {
			const point_t *p = &Point[n];
			double e, d, r[3];
			if (n == 0) { e = EX(p->i, p->j, p->k); d = p->dx; r[0] = Xc[p->i]; r[1] = Yn[p->j]; r[2] = Zn[p->k]; }
			else if (n == 1) { e = EY(p->i, p->j, p->k); d = p->dy; r[0] = Xn[p->i]; r[1] = Yc[p->j]; r[2] = Zn[p->k]; }
			else { e = EZ(p->i, p->j, p->k); d = p->dz; r[0] = Xn[p->i]; r[1] = Yn[p->j]; r[2] = Zc[p->k]; }
			if (pw) {
				real_t fi, dfi;
				finc(r[0], r[1], r[2], (itime + 1) * Dt, Planewave.r0, Planewave.ri, Planewave.ei[n], Planewave.ai, Dt, &fi, &dfi);
				CHECK(fabs(fi) > 0.05);
				e += fi;
			}
			const int64_t id = n * (Solver.maxiter + 1) + itime;
			CHECK_NEAR(VPoint[id], -e * d, 1e-6 * fabs(e * d));
		}
	}

	free(Point); free(VPoint);
	Point = NULL; VPoint = NULL;
	NPoint = 0;
	IPlanewave = 0;
	memset(&Planewave, 0, sizeof(Planewave));
}

static void set_all(double (*fe)(int, int, int, int), double (*fh)(int, int, int, int))
{
	real_t *E[3] = {Ex, Ey, Ez}, *H[3] = {Hx, Hy, Hz};
	for (int i = -1; i <= N + 1; i++) {
	for (int j = -1; j <= N + 1; j++) {
	for (int k = -1; k <= N + 1; k++) {
		const int64_t q = NA(i, j, k);
		if ((q < 0) || (q >= NN)) continue;
		for (int a = 0; a < 3; a++) {
			E[a][q] = (real_t)fe(a, i, j, k);
			H[a][q] = (real_t)fh(a, i, j, k);
		}
	}
	}
	}
}

static double uniform_e(int a, int i, int j, int k) { (void)i; (void)j; (void)k; return (a == 0) ? 0.3 : ((a == 1) ? -0.5 : 0.2); }
static double uniform_h(int a, int i, int j, int k) { (void)i; (void)j; (void)k; return (a == 0) ? -0.1 : ((a == 1) ? 0.4 : 0.25); }
/* 辺の並び (E_a の他の 2 軸) で符号が交互 → 4 本の平均は 0 */
static double checker_e(int a, int i, int j, int k)
{
	const int s = (a == 0) ? (j + k) : ((a == 1) ? (k + i) : (i + j));
	return (s % 2) ? 1.0 : -1.0;
}
/* 面の並び (H_a の自分の軸) で符号が交互 → 2 枚の平均は 0 */
static double checker_h(int a, int i, int j, int k)
{
	const int s = (a == 0) ? i : ((a == 1) ? j : k);
	return ((s + 4) % 2) ? 1.0 : -1.0;
}
/* 一次関数 E_a = 1 + c_a·(i, j, k) (係数はすべて正、成分ごとに別)。
   4 本の平均はセル中心の値 (自分の軸はセル番号のまま、他の 2 軸は +1/2) */
static const double LIN[3][3] = {{0.5, 1, 2}, {3, 0.25, 1.5}, {0.75, 2.5, 0.4}};
static double linear_e(int a, int i, int j, int k)
{
	return 1 + LIN[a][0] * i + LIN[a][1] * j + LIN[a][2] * k;
}
static double zero_h(int a, int i, int j, int k) { (void)a; (void)i; (void)j; (void)k; return 0; }

static void test_average(void)
{
	double f[2];

	set_all(uniform_e, uniform_h);
	average(f);
	CHECK_NEAR(f[0], 0.3 + 0.5 + 0.2, 1e-6);
	CHECK_NEAR(f[1], 0.1 + 0.4 + 0.25, 1e-6);

	set_all(checker_e, checker_h);
	average(f);
	CHECK_NEAR(f[0], 0, 1e-12);
	CHECK_NEAR(f[1], 0, 1e-12);

	set_all(linear_e, zero_h);
	average(f);
	double expect = 0;
	for (int i = 0; i < N; i++) {
	for (int j = 0; j < N; j++) {
	for (int k = 0; k < N; k++) {
		for (int a = 0; a < 3; a++) {
			double c[3] = {i + 0.5, j + 0.5, k + 0.5};
			c[a] -= 0.5;                      /* 自分の軸はセル番号 */
			expect += 1 + LIN[a][0] * c[0] + LIN[a][1] * c[1] + LIN[a][2] * c[2];
		}
	}
	}
	}
	expect /= (double)N * N * N;
	CHECK_NEAR(f[0], expect, 1e-5 * expect);
	CHECK_NEAR(f[1], 0, 1e-12);
}

void test_observe(void)
{
	const double dt0 = Dt, tw0 = Tw;

	ob_setup();
	test_setupdft();
	test_vpoint();
	test_average();
	ob_teardown();

	Dt = dt0;
	Tw = tw0;
}
