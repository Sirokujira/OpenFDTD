/*
test_outputcoupling.c

sol/outputCoupling.c (給電と観測点の結合度) のユニットテスト。

    C_pf = DFT(V_point) / DFT(V_feed)

GUI が読む /metadata/coupling と ofd.log の結合度表のもと。
振幅[dB] と位相[deg] への変換は dbDeg() (sol/utils.c) で、
S パラメータと共通 (tests/test_outputspara.c で検証している)。

## 時系列の作り方

2 つの周波数 f1, f2 の正弦波を重ねる。どちらも整数周期になるように
取ると両者は DFT 上で直交するので、
    DFT_f(V) = (Ntime/2) * (その周波数成分の複素振幅)
となり、結合度は複素振幅の比として解析的に決まる。
周波数ごとに違う比を与えておけば、周波数の添字を取り違えたときに
値が変わる。給電・観測点もそれぞれ 2 本にして、行の取り違えも捕まえる。
*/

#include "ofd.h"
#include "complex.h"
#include "ofd_test.h"

#include <stdlib.h>

d_complex_t coupling(int ifeed, int ipoint, int ifreq);

#define NT (2000)
#define F1 (1.0e9)      /* 1 周期 40 点 -> 50 周期 */
#define F2 (2.0e9)      /* 1 周期 20 点 -> 100 周期 */

/* 2 周波数の成分を持つ時系列を書く */
static void put(double *v, d_complex_t c1, d_complex_t c2)
{
	for (int n = 0; n <= Solver.maxiter; n++) {
		const double t = n * Dt;
		v[n] = d_abs(c1) * cos((2 * PI * F1 * t) + atan2(c1.i, c1.r))
		     + d_abs(c2) * cos((2 * PI * F2 * t) + atan2(c2.i, c2.r));
	}
}

void test_outputcoupling(void)
{
	NFeed = 2;
	NPoint = 2;
	NFreq1 = 2;
	Ntime = NT;
	Solver.maxiter = NT;
	Dt = 1.0 / (F1 * 40);

	Freq1 = (double *)malloc(2 * sizeof(double));
	Freq1[0] = F1;
	Freq1[1] = F2;

	const size_t len = (size_t)(Solver.maxiter + 1);
	VFeed  = (double *)malloc(2 * len * sizeof(double));
	VPoint = (double *)malloc((2 + 2) * len * sizeof(double));

	/* 給電 2 本・観測点 2 本、それぞれ周波数ごとに別の複素振幅 */
	const d_complex_t vf[2][2] = {               /* [feed][freq] */
		{{1.0, 0.5}, {0.8, -0.3}},
		{{-0.6, 0.9}, {0.4, 0.7}},
	};
	const d_complex_t vp[2][2] = {               /* [point][freq] */
		{{0.2, -0.1}, {-0.35, 0.15}},
		{{0.05, 0.45}, {0.6, -0.25}},
	};
	for (int i = 0; i < 2; i++) {
		put(&VFeed[i * len],  vf[i][0], vf[i][1]);
		put(&VPoint[i * len], vp[i][0], vp[i][1]);
	}

	/* 全組み合わせで C = Vp/Vf を独立に計算して比較 */
	for (int f = 0; f < 2; f++) {
		for (int p = 0; p < 2; p++) {
			for (int k = 0; k < 2; k++) {
				const d_complex_t want = d_div(vp[p][k], vf[f][k]);
				const d_complex_t got  = coupling(f, p, k);
				CHECK_NEAR(got.r, want.r, 1e-9);
				CHECK_NEAR(got.i, want.i, 1e-9);
			}
		}
	}

	/* 観測点の電圧が給電と同じなら結合度は 1 (0 dB, 0 deg) */
	put(&VPoint[0], vf[1][0], vf[1][1]);
	{
		const d_complex_t c = coupling(1, 0, 0);
		CHECK_NEAR(c.r, 1.0, 1e-9);
		CHECK_NEAR(c.i, 0.0, 1e-9);
	}

	free(Freq1); free(VFeed); free(VPoint);
	Freq1 = NULL; VFeed = NULL; VPoint = NULL;
	NFeed = NPoint = NFreq1 = 0;
}
