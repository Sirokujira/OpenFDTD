/*
test_outputcross.c

sol/outputCross.c (散乱断面積 BCS / FCS) のユニットテスト。

calcCross() は平面波の入射方向パラメータ (θ, φ) について
    後方 BCS = |E|^2 at (θ, φ)            … 波源へ戻る方向
    前方 FCS = |E|^2 at (180-θ, 180+φ)    … 対蹠点 = 波の進む方向
を求める。(θ, φ) は**波が来る方向**を表す (sol/input2.c で伝搬ベクトルを
ri = -r1 としている)。

検証の根拠:
 1. 対蹠点 … 前方の方向は単位ベクトル -r̂ の方向として、
    θ' = acos(-cosθ), φ' = atan2(-sinθ sinφ, -sinθ cosφ) から**独立に**
    求めて遠方界を計算し、calcCross の結果と比べる。角度の式
    (180-θ, 180+φ) を写していない。
 2. 物理 … 面に Huygens 源 (Ey = 1, Z0 Hz = 1, 法線 x̂) を置くと
    +x 方向にだけ放射し、-x 方向には打ち消しで 0 になる
    (tests/test_farfield.c で確認済みの性質: 前方 |E| = 2, 後方 0)。
    波が -x 側から来て +x へ進む配置 (θ, φ) = (90, 180) にすると
        BCS (-x 方向) = 0,   FCS (+x 方向) = ffctr^2 * 2^2
    になる。平面波の正規化 ffctr = k / sqrt(4π) で、λ = 1 m なら
    k = 2π より ffctr^2 = π、FCS = 4π [m^2]。
 3. 正規化 … BCS = ffctr^2 * |E(ffctr=1)|^2 であること。
*/

#include "ofd.h"
#include "complex.h"
#include "ofd_test.h"

#include <stdlib.h>

void calcCross(int ifreq, double *bcs, double *fcs);
void farfield(int ifreq, double theta, double phi, double ffctr, d_complex_t *etheta, d_complex_t *ephi);

#define NS (4)

static void alloc_surface(void)
{
	NFreq2 = 1;
	Freq2 = (double *)malloc(sizeof(double));
	Freq2[0] = C;                         /* λ = 1 m, k = 2π */
	NSurface = NS;
	Surface = (surface_t *)calloc(NS, sizeof(surface_t));
	SurfaceEx = (d_complex_t **)malloc(sizeof(d_complex_t *));
	SurfaceEy = (d_complex_t **)malloc(sizeof(d_complex_t *));
	SurfaceEz = (d_complex_t **)malloc(sizeof(d_complex_t *));
	SurfaceHx = (d_complex_t **)malloc(sizeof(d_complex_t *));
	SurfaceHy = (d_complex_t **)malloc(sizeof(d_complex_t *));
	SurfaceHz = (d_complex_t **)malloc(sizeof(d_complex_t *));
	SurfaceEx[0] = (d_complex_t *)calloc(NS, sizeof(d_complex_t));
	SurfaceEy[0] = (d_complex_t *)calloc(NS, sizeof(d_complex_t));
	SurfaceEz[0] = (d_complex_t *)calloc(NS, sizeof(d_complex_t));
	SurfaceHx[0] = (d_complex_t *)calloc(NS, sizeof(d_complex_t));
	SurfaceHy[0] = (d_complex_t *)calloc(NS, sizeof(d_complex_t));
	SurfaceHz[0] = (d_complex_t *)calloc(NS, sizeof(d_complex_t));

	/* 平面波入射の正規化 (farfactor) を使う */
	NFeed = 0;
	IPlanewave = 1;
}

static void free_surface(void)
{
	free(SurfaceEx[0]); free(SurfaceEy[0]); free(SurfaceEz[0]);
	free(SurfaceHx[0]); free(SurfaceHy[0]); free(SurfaceHz[0]);
	free(SurfaceEx); free(SurfaceEy); free(SurfaceEz);
	free(SurfaceHx); free(SurfaceHy); free(SurfaceHz);
	free(Surface); free(Freq2);
	IPlanewave = 0; NSurface = 0; NFreq2 = 0;
}

static void patch(int n, double nx, double ny, double nz,
                  double x, double y, double z,
                  double ex, double ey, double ez,
                  double hx, double hy, double hz)
{
	Surface[n].nx = nx; Surface[n].ny = ny; Surface[n].nz = nz;
	Surface[n].x = x;   Surface[n].y = y;   Surface[n].z = z;
	Surface[n].ds = 1.0;
	SurfaceEx[0][n] = d_complex(ex, 0); SurfaceEy[0][n] = d_complex(ey, 0);
	SurfaceEz[0][n] = d_complex(ez, 0); SurfaceHx[0][n] = d_complex(hx, 0);
	SurfaceHy[0][n] = d_complex(hy, 0); SurfaceHz[0][n] = d_complex(hz, 0);
}

static double e2(double th, double ph, double ffctr)
{
	d_complex_t et, ep;
	farfield(0, th, ph, ffctr, &et, &ep);
	return d_norm(et) + d_norm(ep);
}

/* 1. 物理 : 前方にだけ放射する Huygens 源で BCS = 0, FCS = 4π */
static void test_huygens(void)
{
	alloc_surface();
	NSurface = 1;
	patch(0, 1, 0, 0, 0, 0, 0,  0, 1, 0,  0, 0, 1);

	/* 波は -x 側から来て +x へ進む = 来る方向 (θ, φ) = (90, 180) */
	Planewave.theta = 90;
	Planewave.phi = 180;

	double bcs, fcs;
	calcCross(0, &bcs, &fcs);
	CHECK_NEAR(bcs, 0.0, 1e-12);          /* 後方 (-x) は打ち消し */
	CHECK_NEAR(fcs, 4 * PI, 1e-9);        /* 前方 (+x) は ffctr^2 * 4 = 4π */

	free_surface();
}

/* 2. 対蹠点と正規化 : 一般の (非対称な) 源で、前方の方向を -r̂ から求める */
static void test_antipode(void)
{
	alloc_surface();
	/* 位置・法線・成分をばらした 4 枚。どの方向にも対称性が無い */
	patch(0, 1, 0, 0,  0.10, 0.20, 0.30,  0, 0.3, 1,    0, 1, -0.2);
	patch(1, 0, 1, 0, -0.15, 0.05, 0.25,  1, 0, -0.4,   0.5, 0, 1);
	patch(2, 0, 0, 1,  0.22, -0.3, 0.05,  0.2, 0.7, 0,  1, 0.3, 0);
	patch(3, -1, 0, 0, -0.4, 0.1, -0.1,   0, -0.6, 0.5, 0, 0.2, 0.9);

	const double angles[3][2] = {{50, 35}, {120, 300}, {20, 200}};
	for (int a = 0; a < 3; a++) {
		const double th = angles[a][0], ph = angles[a][1];
		Planewave.theta = th;
		Planewave.phi = ph;

		double bcs, fcs;
		calcCross(0, &bcs, &fcs);

		/* 正規化 : ffctr = k/sqrt(4π), k = 2π -> ffctr^2 = π */
		CHECK_NEAR(bcs, PI * e2(th, ph, 1.0), 1e-9 * (1 + bcs));

		/* 前方 = 単位ベクトル -r̂ の方向を独立に求める */
		const double st = sin(th * DTOR), ct = cos(th * DTOR);
		const double sp = sin(ph * DTOR), cp = cos(ph * DTOR);
		const double rx = -st * cp, ry = -st * sp, rz = -ct;
		const double th2 = acos(rz) / DTOR;
		const double ph2 = atan2(ry, rx) / DTOR;
		CHECK_NEAR(fcs, PI * e2(th2, ph2, 1.0), 1e-9 * (1 + fcs));

		/* 前方と後方は別物であること (源が非対称なので) */
		CHECK(fabs(fcs - bcs) > 1e-6);
	}

	free_surface();
}

void test_outputcross(void)
{
	test_huygens();
	test_antipode();
}
