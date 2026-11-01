/*
test_post.c

後処理 ofd_post の計算 (post/calcFar1d.c, post/calcFar2d.c, post/calcNear1d.c,
post/calcNear2d.c, post/planewave.c) のテスト。

ここで作るのは far1d.log / far2d.log / near1d.log / near2d.log と作図の値で、
壊れても ofd_post は普通に終わり、図の向きや値がおかしいことにしか現れない。

## 遠方界 (calcFar1d / calcFar2d)

x 方向の微小電流 (等価面 1 枚 : n̂ = ŷ, H = ẑ -> J = n̂×H = x̂) の放射電力は
    |Eθ|^2 + |Eφ|^2  ∝  |x̂ の視線に垂直な成分|^2 = 1 - (sinθ cosφ)^2
で、面ごとに形が違う。入力の約束 (README / sol/input_data の plotfar1d) は
    X 面 (yz 面) : θ = 角度, φ = 90   -> 1                (一定)
    Y 面 (xz 面) : θ = 角度, φ = 0    -> cos^2 θ
    Z 面 (xy 面) : θ = 90,   φ = 角度 -> sin^2 φ
    V 面 : θ = 角度, φ = 指定         -> 1 - sin^2θ cos^2φ0
    H 面 : θ = 指定, φ = 角度         -> 1 - sin^2θ0 cos^2φ
角度は 0..360 を div 等分 (両端を含む div+1 点)。far1dnorm では各面の最大が 1。
far2d は θ = 0..180 を divtheta、φ = 0..360 を divphi 等分する。

## 近傍界 (calcNear1d / calcNear2d)

線・面の上の節点 (i,j,k) を周波数ごとに並べ、節点補間 NodeE_c / NodeH_c
(tests/test_farfield.c で検証済み) の値を入れる。ここで確かめるのは
**並びと (i,j,k) の割り当て** : 線の位置は (README の plotnear1d の約束)
    X 方向の線 : (y, z)    Y 方向の線 : (z, x)    Z 方向の線 : (x, y)
面は X 面 = (j,k)、Y 面 = (i,k)、Z 面 = (i,j) の順で、外側が 1 番目。
場は成分ごとに x, y, z で傾きの違う一次関数にして、添字の取り違えが値に出るようにする。

## 平面波の入射界 (planewave)

e^{jωt} の約束で +r̂i 方向へ進む波は exp(-jk (r - r0)·r̂i) (r0 は領域の中心)。
中心で Ei そのもの、進行方向に λ/4 進むと -j 倍、垂直方向には位相が変わらない。
*/

#include "ofd.h"
#include "complex.h"
#include "ofd_prototype.h"
#include "ofd_test.h"

#include <stdlib.h>
#include <string.h>

void   alloc_farfield(void);
void   NodeE_c(int ifreq, int i, int j, int k, d_complex_t *cex, d_complex_t *cey, d_complex_t *cez);
void   NodeH_c(int ifreq, int i, int j, int k, d_complex_t *chx, d_complex_t *chy, d_complex_t *chz);

/* ---------------------------------------------------------------- */
/* 遠方界                                                           */
/* ---------------------------------------------------------------- */

static void setup_xdipole(void)
{
	NFreq2 = 1;
	Freq2 = (double *)malloc(sizeof(double));
	Freq2[0] = C;          /* λ = 1 m */
	NSurface = 1;
	Surface = (surface_t *)calloc(1, sizeof(surface_t));
	SurfaceEx = (d_complex_t **)malloc(sizeof(d_complex_t *));
	SurfaceEy = (d_complex_t **)malloc(sizeof(d_complex_t *));
	SurfaceEz = (d_complex_t **)malloc(sizeof(d_complex_t *));
	SurfaceHx = (d_complex_t **)malloc(sizeof(d_complex_t *));
	SurfaceHy = (d_complex_t **)malloc(sizeof(d_complex_t *));
	SurfaceHz = (d_complex_t **)malloc(sizeof(d_complex_t *));
	SurfaceEx[0] = (d_complex_t *)calloc(1, sizeof(d_complex_t));
	SurfaceEy[0] = (d_complex_t *)calloc(1, sizeof(d_complex_t));
	SurfaceEz[0] = (d_complex_t *)calloc(1, sizeof(d_complex_t));
	SurfaceHx[0] = (d_complex_t *)calloc(1, sizeof(d_complex_t));
	SurfaceHy[0] = (d_complex_t *)calloc(1, sizeof(d_complex_t));
	SurfaceHz[0] = (d_complex_t *)calloc(1, sizeof(d_complex_t));
	Surface[0].ny = 1;
	Surface[0].ds = 1;
	SurfaceHz[0][0] = d_complex(1, 0);
}

static void free_xdipole(void)
{
	free(SurfaceEx[0]); free(SurfaceEy[0]); free(SurfaceEz[0]);
	free(SurfaceHx[0]); free(SurfaceHy[0]); free(SurfaceHz[0]);
	free(SurfaceEx); free(SurfaceEy); free(SurfaceEz);
	free(SurfaceHx); free(SurfaceHy); free(SurfaceHz);
	free(Surface); free(Freq2);
	Surface = NULL; Freq2 = NULL;
	SurfaceEx = SurfaceEy = SurfaceEz = SurfaceHx = SurfaceHy = SurfaceHz = NULL;
	NSurface = 0; NFreq2 = 0;
}

/* x 方向の微小電流の相対放射電力 */
static double xpower(double theta, double phi)
{
	const double s = sin(theta * DTOR) * cos(phi * DTOR);
	return 1 - (s * s);
}

static double power(d_complex_t et, d_complex_t ep) { return d_norm(et) + d_norm(ep); }

static void test_far1d(void)
{
	setup_xdipole();
	const double ffctr[1] = {3.0};   /* 1 だと基準電力が 1 になり、正規化の有無が値に出ない */

	/* 基準 : θ = 90, φ = 90 (x 方向に垂直) の電力。xpower = 1 に対応する */
	d_complex_t et0, ep0;
	farfield(0, 90, 90, ffctr[0], &et0, &ep0);
	const double p0 = power(et0, ep0);
	CHECK(p0 > 0);

	const char dirs[5] = {'X', 'Y', 'Z', 'V', 'H'};
	const double fixed[5] = {0, 0, 0, 60, 30};
	const int div = 12;              /* 30 度刻み */
	NFar1d = 5;
	Far1d = (far1d_t *)calloc(NFar1d, sizeof(far1d_t));
	for (int n = 0; n < NFar1d; n++) {
		Far1d[n].dir = dirs[n];
		Far1d[n].div = div;
		Far1d[n].angle = fixed[n];
	}
	d_complex_t ***et = (d_complex_t ***)malloc(NFar1d * sizeof(d_complex_t **));
	d_complex_t ***ep = (d_complex_t ***)malloc(NFar1d * sizeof(d_complex_t **));
	for (int n = 0; n < NFar1d; n++) {
		et[n] = (d_complex_t **)malloc(sizeof(d_complex_t *));
		ep[n] = (d_complex_t **)malloc(sizeof(d_complex_t *));
		et[n][0] = (d_complex_t *)calloc(div + 1, sizeof(d_complex_t));
		ep[n][0] = (d_complex_t *)calloc(div + 1, sizeof(d_complex_t));
	}

	for (int norm = 0; norm < 2; norm++) {
		Far1dNorm = norm;
		calcFar1d(ffctr, et, ep);
		for (int n = 0; n < NFar1d; n++) {
			double pmax = 0;
			for (int i = 0; i <= div; i++) {
				const double a = 360.0 * i / div;
				double th = 0, ph = 0;
				switch (dirs[n]) {
				case 'X': th = a;        ph = 90;       break;
				case 'Y': th = a;        ph = 0;        break;
				case 'Z': th = 90;       ph = a;        break;
				case 'V': th = a;        ph = fixed[n]; break;
				case 'H': th = fixed[n]; ph = a;        break;
				}
				pmax = MAX(pmax, xpower(th, ph));
			}
			int bad = 0;
			for (int i = 0; i <= div; i++) {
				const double a = 360.0 * i / div;
				double th = 0, ph = 0;
				switch (dirs[n]) {
				case 'X': th = a;        ph = 90;       break;
				case 'Y': th = a;        ph = 0;        break;
				case 'Z': th = 90;       ph = a;        break;
				case 'V': th = a;        ph = fixed[n]; break;
				case 'H': th = fixed[n]; ph = a;        break;
				}
				/* 正規化なし : 基準との比。あり : 面の最大で割った値 */
				const double expect = norm ? (xpower(th, ph) / pmax) : (xpower(th, ph) * p0);
				const double got = power(et[n][0][i], ep[n][0][i]);
				if (fabs(got - expect) > 1e-9 * (expect + p0)) bad++;
			}
			CHECK(bad == 0);
			if (bad) printf("    far1d dir=%c norm=%d : %d points differ\n", dirs[n], norm, bad);
		}
	}

	for (int n = 0; n < NFar1d; n++) {
		free(et[n][0]); free(ep[n][0]); free(et[n]); free(ep[n]);
	}
	free(et); free(ep);
	free(Far1d); Far1d = NULL; NFar1d = 0; Far1dNorm = 0;
	free_xdipole();
}

static void test_far2d(void)
{
	setup_xdipole();
	const double ffctr[1] = {3.0};   /* 1 だと基準電力が 1 になり、正規化の有無が値に出ない */
	d_complex_t et0, ep0;
	farfield(0, 90, 90, ffctr[0], &et0, &ep0);
	const double p0 = power(et0, ep0);

	Far2d.divtheta = 6;     /* 30 度 */
	Far2d.divphi = 8;       /* 45 度 */
	d_complex_t ***et = (d_complex_t ***)malloc(sizeof(d_complex_t **));
	d_complex_t ***ep = (d_complex_t ***)malloc(sizeof(d_complex_t **));
	et[0] = (d_complex_t **)malloc((Far2d.divtheta + 1) * sizeof(d_complex_t *));
	ep[0] = (d_complex_t **)malloc((Far2d.divtheta + 1) * sizeof(d_complex_t *));
	for (int it = 0; it <= Far2d.divtheta; it++) {
		et[0][it] = (d_complex_t *)calloc(Far2d.divphi + 1, sizeof(d_complex_t));
		ep[0][it] = (d_complex_t *)calloc(Far2d.divphi + 1, sizeof(d_complex_t));
	}
	calcFar2d(ffctr, et, ep);
	int bad = 0;
	for (int it = 0; it <= Far2d.divtheta; it++) {
	for (int ip = 0; ip <= Far2d.divphi; ip++) {
		const double expect = xpower(180.0 * it / Far2d.divtheta, 360.0 * ip / Far2d.divphi) * p0;
		if (fabs(power(et[0][it][ip], ep[0][it][ip]) - expect) > 1e-9 * p0) bad++;
	}
	}
	CHECK(bad == 0);
	for (int it = 0; it <= Far2d.divtheta; it++) { free(et[0][it]); free(ep[0][it]); }
	free(et[0]); free(ep[0]); free(et); free(ep);
	free_xdipole();
}

/* ---------------------------------------------------------------- */
/* 近傍界                                                           */
/* ---------------------------------------------------------------- */

#define NXP (4)
#define NYP (3)
#define NZP (5)
static float *nbuf[12];

static void setup_near_grid(void)
{
	Nx = NXP; Ny = NYP; Nz = NZP;
	iMin = jMin = kMin = 0;
	iMax = Nx; jMax = Ny; kMax = Nz;
	PBCx = PBCy = PBCz = 0;
	Nk = 1;
	Nj = (Nz + 3) * Nk;
	Ni = (Ny + 3) * Nj;
	N0 = Ni + Nj + Nk;
	NN = NA(Nx + 1, Ny + 1, Nz + 1) + 1;
	Xn = (double *)malloc((Nx + 1) * sizeof(double)); Xc = (double *)malloc(Nx * sizeof(double));
	Yn = (double *)malloc((Ny + 1) * sizeof(double)); Yc = (double *)malloc(Ny * sizeof(double));
	Zn = (double *)malloc((Nz + 1) * sizeof(double)); Zc = (double *)malloc(Nz * sizeof(double));
	for (int i = 0; i <= Nx; i++) Xn[i] = 0.1 * i;
	for (int j = 0; j <= Ny; j++) Yn[j] = 0.1 * j;
	for (int k = 0; k <= Nz; k++) Zn[k] = 0.1 * k;
	for (int i = 0; i < Nx; i++) Xc[i] = 0.1 * (i + 0.5);
	for (int j = 0; j < Ny; j++) Yc[j] = 0.1 * (j + 0.5);
	for (int k = 0; k < Nz; k++) Zc[k] = 0.1 * (k + 0.5);

	NFreq2 = 2;
	Freq2 = (double *)malloc(2 * sizeof(double));
	Freq2[0] = 1e9; Freq2[1] = 2e9;
	for (int b = 0; b < 12; b++) nbuf[b] = (float *)calloc((size_t)NFreq2 * NN, sizeof(float));
	cEx_r = nbuf[0]; cEx_i = nbuf[1]; cEy_r = nbuf[2]; cEy_i = nbuf[3]; cEz_r = nbuf[4]; cEz_i = nbuf[5];
	cHx_r = nbuf[6]; cHx_i = nbuf[7]; cHy_r = nbuf[8]; cHy_i = nbuf[9]; cHz_r = nbuf[10]; cHz_i = nbuf[11];
	/* 成分・周波数ごとに傾きの違う値 (配列の全域、余白も含む) */
	for (int f = 0; f < NFreq2; f++) {
	for (int i = -1; i <= Nx + 1; i++) {
	for (int j = -1; j <= Ny + 1; j++) {
	for (int k = -1; k <= Nz + 1; k++) {
		const int64_t n = (int64_t)f * NN + NA(i, j, k);
		for (int b = 0; b < 12; b++) {
			nbuf[b][n] = (float)((b + 1) * 0.5 + (f + 1) * 7.0 + i * 1.0 + j * 10.0 + k * 100.0 + b * i * 0.25);
		}
	}
	}
	}
	}
}

static void free_near_grid(void)
{
	for (int b = 0; b < 12; b++) free(nbuf[b]);
	cEx_r = cEx_i = cEy_r = cEy_i = cEz_r = cEz_i = NULL;
	cHx_r = cHx_i = cHy_r = cHy_i = cHz_r = cHz_i = NULL;
	free(Xn); free(Yn); free(Zn); free(Xc); free(Yc); free(Zc);
	Xn = Yn = Zn = Xc = Yc = Zc = NULL;
	free(Freq2); Freq2 = NULL;
	NFreq2 = 0;
	Nx = Ny = Nz = 0; NN = 0;
}

static int same_c(d_complex_t a, d_complex_t b) { return (fabs(a.r - b.r) < 1e-6) && (fabs(a.i - b.i) < 1e-6); }

/* 節点 (i,j,k) の 6 成分が配列の adr 番目と一致するか */
static int node_ok(int f, int i, int j, int k, int64_t adr,
	const d_complex_t *ex, const d_complex_t *ey, const d_complex_t *ez,
	const d_complex_t *hx, const d_complex_t *hy, const d_complex_t *hz)
{
	d_complex_t e[3], h[3];
	NodeE_c(f, i, j, k, &e[0], &e[1], &e[2]);
	NodeH_c(f, i, j, k, &h[0], &h[1], &h[2]);
	return same_c(ex[adr], e[0]) && same_c(ey[adr], e[1]) && same_c(ez[adr], e[2])
	    && same_c(hx[adr], h[0]) && same_c(hy[adr], h[1]) && same_c(hz[adr], h[2]);
}

static void test_near1d(void)
{
	setup_near_grid();
	NNear1d = 3;
	Near1d = (near1d_t *)calloc(NNear1d, sizeof(near1d_t));
	/* X 線 : (y,z) = (j=1, k=4)   Y 線 : (z,x) = (k=3, i=2)   Z 線 : (x,y) = (i=4, j=2) */
	Near1d[0].dir = 'X'; Near1d[0].id1 = 1; Near1d[0].id2 = 4;
	Near1d[1].dir = 'Y'; Near1d[1].id1 = 3; Near1d[1].id2 = 2;
	Near1d[2].dir = 'Z'; Near1d[2].id1 = 4; Near1d[2].id2 = 2;
	calcNear1d();

	int64_t adr = 0;
	int ok[3] = {1, 1, 1};
	for (int n = 0; n < NNear1d; n++) {
		const int len = (n == 0) ? Nx : (n == 1) ? Ny : Nz;
		for (int f = 0; f < NFreq2; f++) {
			for (int l = 0; l <= len; l++) {
				int i, j, k;
				if (n == 0)      { i = l; j = 1; k = 4; }
				else if (n == 1) { i = 2; j = l; k = 3; }
				else             { i = 4; j = 2; k = l; }
				ok[n] &= node_ok(f, i, j, k, adr, Near1dEx, Near1dEy, Near1dEz, Near1dHx, Near1dHy, Near1dHz);
				adr++;
			}
		}
	}
	CHECK(ok[0]);
	CHECK(ok[1]);
	CHECK(ok[2]);

	free(Near1dEx); free(Near1dEy); free(Near1dEz); free(Near1dHx); free(Near1dHy); free(Near1dHz);
	Near1dEx = Near1dEy = Near1dEz = Near1dHx = Near1dHy = Near1dHz = NULL;
	free(Near1d); Near1d = NULL; NNear1d = 0;
	free_near_grid();
}

static void test_near2d(void)
{
	setup_near_grid();
	NNear2d = 3;
	Near2d = (near2d_t *)calloc(NNear2d, sizeof(near2d_t));
	Near2d[0].dir = 'X'; Near2d[0].id0 = 3;
	Near2d[1].dir = 'Y'; Near2d[1].id0 = 1;
	Near2d[2].dir = 'Z'; Near2d[2].id0 = 2;
	calcNear2d();

	int64_t adr = 0;
	int ok[3] = {1, 1, 1};
	for (int n = 0; n < NNear2d; n++) {
		const int d1 = (n == 0) ? Ny : Nx;
		const int d2 = (n == 2) ? Ny : Nz;
		for (int f = 0; f < NFreq2; f++) {
			for (int l1 = 0; l1 <= d1; l1++) {
			for (int l2 = 0; l2 <= d2; l2++) {
				int i, j, k;
				if (n == 0)      { i = 3;  j = l1; k = l2; }
				else if (n == 1) { i = l1; j = 1;  k = l2; }
				else             { i = l1; j = l2; k = 2;  }
				ok[n] &= node_ok(f, i, j, k, adr, Near2dEx, Near2dEy, Near2dEz, Near2dHx, Near2dHy, Near2dHz);
				adr++;
			}
			}
		}
	}
	CHECK(ok[0]);
	CHECK(ok[1]);
	CHECK(ok[2]);

	free(Near2dEx); free(Near2dEy); free(Near2dEz); free(Near2dHx); free(Near2dHy); free(Near2dHz);
	Near2dEx = Near2dEy = Near2dEz = Near2dHx = Near2dHy = Near2dHz = NULL;
	free(Near2d); Near2d = NULL; NNear2d = 0;
	free_near_grid();
}

/* ---------------------------------------------------------------- */
/* 平面波の入射界                                                    */
/* ---------------------------------------------------------------- */

static void test_planewave(void)
{
	Nx = 2; Ny = 2; Nz = 2;
	Xn = (double *)malloc(3 * sizeof(double));
	Yn = (double *)malloc(3 * sizeof(double));
	Zn = (double *)malloc(3 * sizeof(double));
	/* 中心 r0 = (0.3, -0.1, 0.5) */
	Xn[0] = 0.1; Xn[1] = 0.3; Xn[2] = 0.5;
	Yn[0] = -0.4; Yn[1] = -0.1; Yn[2] = 0.2;
	Zn[0] = 0.0; Zn[1] = 0.5; Zn[2] = 1.0;
	const double r0[3] = {0.3, -0.1, 0.5};

	/* 進行方向 r̂i = (0.6, 0, 0.8)、Ei ⊥ r̂i、Hi = r̂i × Ei / η0 相当の値 (ここでは任意) */
	planewave_t pw0 = Planewave;
	const int ipw0 = IPlanewave;
	IPlanewave = 1;
	Planewave.ri[0] = 0.6; Planewave.ri[1] = 0.0; Planewave.ri[2] = 0.8;
	Planewave.ei[0] = 0.8; Planewave.ei[1] = 0.0; Planewave.ei[2] = -0.6;
	Planewave.hi[0] = 0.0; Planewave.hi[1] = 2.5; Planewave.hi[2] = 0.0;

	const double freq = C / 2.0;    /* λ = 2 m */
	const double lam = 2.0;
	d_complex_t e[3], h[3];

	/* 中心 : Ei, Hi そのもの */
	planewave(freq, r0[0], r0[1], r0[2], e, h);
	CHECK(same_c(e[0], d_complex(0.8, 0)) && same_c(e[2], d_complex(-0.6, 0)) && same_c(h[1], d_complex(2.5, 0)));

	/* 進行方向に λ/4 : -j 倍 */
	const double q = lam / 4;
	planewave(freq, r0[0] + 0.6 * q, r0[1], r0[2] + 0.8 * q, e, h);
	CHECK(same_c(e[0], d_complex(0, -0.8)) && same_c(e[2], d_complex(0, 0.6)) && same_c(h[1], d_complex(0, -2.5)));

	/* 進行方向に λ : 元に戻る。垂直方向 (ŷ と r̂i×ŷ) には位相が変わらない */
	planewave(freq, r0[0] + 0.6 * lam, r0[1], r0[2] + 0.8 * lam, e, h);
	CHECK(same_c(e[0], d_complex(0.8, 0)));
	planewave(freq, r0[0] + 0.8 * 0.37, r0[1] + 0.9, r0[2] - 0.6 * 0.37, e, h);
	CHECK(same_c(e[0], d_complex(0.8, 0)) && same_c(h[1], d_complex(2.5, 0)));

	/* 平面波でなければ何も書かない */
	IPlanewave = 0;
	e[0] = d_complex(123, 456);
	planewave(freq, 0, 0, 0, e, h);
	CHECK(same_c(e[0], d_complex(123, 456)));

	IPlanewave = ipw0;
	Planewave = pw0;
	free(Xn); free(Yn); free(Zn);
	Xn = Yn = Zn = NULL;
	Nx = Ny = Nz = 0;
}

void test_post(void)
{
	test_far1d();
	test_far2d();
	test_near1d();
	test_near2d();
	test_planewave();
}
