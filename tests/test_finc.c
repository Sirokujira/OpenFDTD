/*
test_finc.c

include/finc.h (入射波形) のユニットテスト。

finc() は平面波入射の全実装 (CPU / MPI / CUDA) が共有する波源であり、
ここが壊れると散乱界が丸ごと狂うが "normal end" では終わってしまう
(実際に長期間見逃されていた。CI の sphere_rcs_check.sh はその再発防止)。
波形そのものを式レベルで固定しておく。

2 つのモードを持つ:
  WaveAmp <= 0 : ガウス微分パルス (既定)
  WaveAmp >  0 : CW 正弦波 + raised-cosine ランプ (TPA 等の非線形解析用)
*/

#include "ofd.h"
#include "finc.h"
#include "ofd_test.h"

/* real_t は既定で float なので、比較は単精度の分解能に合わせる */
#define FTOL (1e-6)

static void test_finc_pulse(void)
{
	const double r0[3] = {0.0, 0.0, 0.0};
	const double ri[3] = {1.0, 0.0, 0.0};   /* +X 伝搬 */
	const double ai = 1.0e10;
	const double dt = 1.0e-12;
	const double fc = 2.0;
	real_t fi, dfi;

	WaveAmp = 0;      /* パルスモード */
	WaveOmega = 0;

	/* t=0, 原点 : at=0 -> fi=0, dfi = dt*ai*fc */
	finc(0.0, 0.0, 0.0, 0.0, r0, ri, fc, ai, dt, &fi, &dfi);
	CHECK_NEAR(fi, 0.0, FTOL);
	CHECK_NEAR(dfi, dt * ai * fc, FTOL);

	/* 一般の時刻で定義式と一致すること */
	{
		const double t = 3.0e-11;
		const double at = ai * t;
		const double ex = exp(-at * at);
		finc(0.0, 0.0, 0.0, t, r0, ri, fc, ai, dt, &fi, &dfi);
		CHECK_NEAR(fi, at * ex * fc, FTOL);
		CHECK_NEAR(dfi, dt * ai * (1 - 2 * at * at) * ex * fc, FTOL);
	}

	/* (a*t)^2 >= 16 では 0 に打ち切られる */
	{
		const double t_big = 100.0 / ai;
		finc(0.0, 0.0, 0.0, t_big, r0, ri, fc, ai, dt, &fi, &dfi);
		CHECK(fi == (real_t)0);
		CHECK(dfi == (real_t)0);
	}

	/* 伝搬遅延 : +X に x 進んだ点では t を x/c だけ遅らせると同じ値になる */
	{
		const double c = 2.99792458e8;
		const double t = 2.0e-11;
		const double x = 0.05;
		real_t fi0, dfi0, fi1, dfi1;
		finc(0.0, 0.0, 0.0, t, r0, ri, fc, ai, dt, &fi0, &dfi0);
		finc(x, 0.0, 0.0, t + x / c, r0, ri, fc, ai, dt, &fi1, &dfi1);
		CHECK_NEAR(fi1, fi0, FTOL);
		CHECK_NEAR(dfi1, dfi0, FTOL);
	}

	/* 伝搬方向に直交する向きには遅延が生じない */
	{
		const double t = 2.0e-11;
		real_t fi0, dfi0, fi1, dfi1;
		finc(0.0, 0.0, 0.0, t, r0, ri, fc, ai, dt, &fi0, &dfi0);
		finc(0.0, 0.3, -0.2, t, r0, ri, fc, ai, dt, &fi1, &dfi1);
		CHECK_NEAR(fi1, fi0, FTOL);
		CHECK_NEAR(dfi1, dfi0, FTOL);
	}

	/* dfi = dt * d(fi)/dt であること (中心差分で確認) */
	{
		const double t = 1.5e-11;
		const double h = 1.0e-14;
		real_t fp, fm, dummy, d0;
		finc(0.0, 0.0, 0.0, t + h, r0, ri, fc, ai, dt, &fp, &dummy);
		finc(0.0, 0.0, 0.0, t - h, r0, ri, fc, ai, dt, &fm, &dummy);
		finc(0.0, 0.0, 0.0, t, r0, ri, fc, ai, dt, &dummy, &d0);
		CHECK_REL(d0, dt * (fp - fm) / (2 * h), 1e-3);
	}

	/* 振幅係数 fc について線形であること */
	{
		const double t = 2.0e-11;
		real_t f1, d1, f2, d2;
		finc(0.0, 0.0, 0.0, t, r0, ri, 1.0, ai, dt, &f1, &d1);
		finc(0.0, 0.0, 0.0, t, r0, ri, 3.0, ai, dt, &f2, &d2);
		CHECK_REL(f2, 3.0 * f1, 1e-5);
		CHECK_REL(d2, 3.0 * d1, 1e-5);
	}
}

static void test_finc_cw(void)
{
	const double r0[3] = {0.0, 0.0, 0.0};
	const double ri[3] = {0.0, 0.0, 1.0};   /* +Z 伝搬 */
	const double ai = 1.0e10;               /* CW では使われない */
	const double dt = 1.0e-14;
	const double fc = 1.0;
	const double e0 = 5.0e6;
	const double omega = 2 * PI * 1.0e14;
	real_t fi, dfi;

	WaveAmp = e0;
	WaveOmega = omega;

	{
		const double tr = FINC_CW_RAMP_CYCLES * (2 * PI / omega);

		/* t <= 0 は 0 (立ち上がり前) */
		finc(0.0, 0.0, 0.0, 0.0, r0, ri, fc, ai, dt, &fi, &dfi);
		CHECK(fi == (real_t)0);
		CHECK(dfi == (real_t)0);

		/* ランプ区間の途中で定義式と一致すること */
		{
			const double t = 0.3 * tr;
			const double s = 0.5 * (1 - cos(PI * t / tr));
			const double ds = 0.5 * (PI / tr) * sin(PI * t / tr);
			finc(0.0, 0.0, 0.0, t, r0, ri, fc, ai, dt, &fi, &dfi);
			CHECK_REL(fi, e0 * fc * s * sin(omega * t), 1e-5);
			CHECK_REL(dfi, dt * e0 * fc * ((ds * sin(omega * t))
				+ (s * omega * cos(omega * t))), 1e-5);
		}

		/* ランプ終了後は包絡線が 1 で飽和する : |fi| <= E0 かつ
		   sin が 1 になる位相でちょうど E0 になる */
		{
			const double t = 3 * tr + (0.25 * (2 * PI / omega));
			finc(0.0, 0.0, 0.0, t, r0, ri, fc, ai, dt, &fi, &dfi);
			CHECK_REL(fi, e0 * sin(omega * t), 1e-4);
			CHECK(fabs((double)fi) <= e0 * (1 + 1e-5));
		}

		/* ランプの接続点 t=tr で連続であること (不連続だと広帯域雑音になる) */
		{
			const double h = 1.0e-6 * tr;
			real_t fa, fb, dummy;
			finc(0.0, 0.0, 0.0, tr - h, r0, ri, fc, ai, dt, &fa, &dummy);
			finc(0.0, 0.0, 0.0, tr + h, r0, ri, fc, ai, dt, &fb, &dummy);
			CHECK_NEAR(fa, fb, 1e-3 * e0);
		}

		/* 立ち上がり直後は振幅が小さい (ランプが効いている) */
		{
			const double t = 0.02 * tr;
			finc(0.0, 0.0, 0.0, t, r0, ri, fc, ai, dt, &fi, &dfi);
			CHECK(fabs((double)fi) < 0.1 * e0);
		}

		/* dfi = dt * d(fi)/dt (CW 分岐でも成立すること) */
		{
			const double t = 2 * tr;
			const double h = 1.0e-4 / omega;
			real_t fp, fm, dummy, d0;
			finc(0.0, 0.0, 0.0, t + h, r0, ri, fc, ai, dt, &fp, &dummy);
			finc(0.0, 0.0, 0.0, t - h, r0, ri, fc, ai, dt, &fm, &dummy);
			finc(0.0, 0.0, 0.0, t, r0, ri, fc, ai, dt, &dummy, &d0);
			CHECK_REL(d0, dt * (fp - fm) / (2 * h), 1e-2);
		}

		/* +Z 伝搬の遅延が z 方向にのみ効くこと */
		{
			const double c = 2.99792458e8;
			const double t = 2 * tr;
			const double z = 1.0e-6;
			real_t f0, f1, dummy;
			finc(0.0, 0.0, 0.0, t, r0, ri, fc, ai, dt, &f0, &dummy);
			finc(0.0, 0.0, z, t + z / c, r0, ri, fc, ai, dt, &f1, &dummy);
			CHECK_REL(f1, f0, 1e-4);
		}
	}

	/* 後続のテストに影響しないよう既定値へ戻す */
	WaveAmp = 0;
	WaveOmega = 0;
}

void test_finc(void)
{
	test_finc_pulse();
	test_finc_cw();
}
