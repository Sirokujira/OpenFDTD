/*
test_vfeed.c

sol/vfeed.c の給電電圧波形のユニットテスト。

    V(t) = sqrt(2) * exp(0.5) * a * exp(-a^2),  a = (t - tw - td) / (tw/4)

前置係数 sqrt(2)*exp(0.5) はピーク振幅がちょうど 1 になるための正規化。
給電波形が変わると入力インピーダンス Zin が丸ごとずれるので、
係数と時間軸の定義を式レベルで固定しておく。
*/

#include "ofd.h"
#include "ofd_test.h"

double vfeed(double t, double tw, double td);

void test_vfeed(void)
{
	const double tw = 1.0e-10;
	const double td = 2.0e-11;

	/* 中心 t = tw + td でゼロ交差する */
	CHECK_NEAR(vfeed(tw + td, tw, td), 0.0, 1e-15);

	/* ピークは a = 1/sqrt(2) すなわち t = tw + td + (tw/4)/sqrt(2) で 1 */
	{
		const double tp = tw + td + (tw / 4) / sqrt(2.0);
		CHECK_NEAR(vfeed(tp, tw, td), 1.0, 1e-12);
	}

	/* ピークが全体の最大値であること */
	{
		const double tp = tw + td + (tw / 4) / sqrt(2.0);
		const double peak = vfeed(tp, tw, td);
		int i;
		for (i = 0; i <= 2000; i++) {
			const double t = i * (4 * tw) / 2000;
			CHECK(vfeed(t, tw, td) <= peak + 1e-12);
		}
	}

	/* 中心に関して奇対称 : V(c+d) = -V(c-d) */
	{
		const double c = tw + td;
		const double d = 1.5e-11;
		CHECK_NEAR(vfeed(c + d, tw, td), -vfeed(c - d, tw, td), 1e-14);
	}

	/* 定義式と一致すること */
	{
		const double t = 7.0e-11;
		const double a = (t - tw - td) / (tw / 4);
		CHECK_NEAR(vfeed(t, tw, td), sqrt(2.0) * exp(0.5) * a * exp(-a * a), 1e-14);
	}

	/* 遅延 td は波形を時間軸上で平行移動するだけであること */
	{
		const double t = 5.0e-11;
		CHECK_NEAR(vfeed(t + td, tw, td), vfeed(t, tw, 0.0), 1e-14);
	}

	/* 十分離れた時刻では減衰して 0 に近づく */
	CHECK(fabs(vfeed(tw + td + 10 * tw, tw, td)) < 1e-10);
	CHECK(fabs(vfeed(tw + td - 10 * tw, tw, td)) < 1e-10);
}
