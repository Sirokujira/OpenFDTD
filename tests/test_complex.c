/*
test_complex.c

include/complex.h の複素数ヘルパーのユニットテスト。

このヘッダは遠方界・近傍界・Zin など後処理の数値計算すべての土台なので、
四則演算・絶対値・偏角・初等関数が定義どおりであることを個別に確認する。
期待値は複素数の定義から独立に導出している (実装式の写しではない)。
*/

#include "ofd.h"
#include "complex.h"
#include "ofd_test.h"

#define TOL (1e-12)

void test_complex(void)
{
	const d_complex_t a = d_complex(3.0, 4.0);
	const d_complex_t b = d_complex(1.0, -2.0);

	/* 生成 */
	CHECK_NEAR(a.r, 3.0, TOL);
	CHECK_NEAR(a.i, 4.0, TOL);

	/* 加減算 */
	CHECK_NEAR(d_add(a, b).r, 4.0, TOL);
	CHECK_NEAR(d_add(a, b).i, 2.0, TOL);
	CHECK_NEAR(d_sub(a, b).r, 2.0, TOL);
	CHECK_NEAR(d_sub(a, b).i, 6.0, TOL);

	/* 乗算 : (3+4i)(1-2i) = 3-6i+4i+8 = 11-2i */
	CHECK_NEAR(d_mul(a, b).r, 11.0, TOL);
	CHECK_NEAR(d_mul(a, b).i, -2.0, TOL);

	/* 除算 : a/b の検算は (a/b)*b == a で行う */
	{
		const d_complex_t q = d_div(a, b);
		const d_complex_t p = d_mul(q, b);
		CHECK_NEAR(p.r, a.r, 1e-10);
		CHECK_NEAR(p.i, a.i, 1e-10);
	}

	/* 0 除算は (0,0) を返す (実装が定めた安全側の挙動) */
	{
		const d_complex_t q = d_div(a, d_complex(0.0, 0.0));
		CHECK_NEAR(q.r, 0.0, TOL);
		CHECK_NEAR(q.i, 0.0, TOL);
	}

	/* 実数倍 */
	CHECK_NEAR(d_rmul(2.0, a).r, 6.0, TOL);
	CHECK_NEAR(d_rmul(2.0, a).i, 8.0, TOL);

	/* 多項加算 */
	{
		const d_complex_t s3 = d_add3(a, b, d_complex(1.0, 1.0));
		const d_complex_t s4 = d_add4(a, b, d_complex(1.0, 1.0), d_complex(-5.0, -3.0));
		CHECK_NEAR(s3.r, 5.0, TOL);
		CHECK_NEAR(s3.i, 3.0, TOL);
		CHECK_NEAR(s4.r, 0.0, TOL);
		CHECK_NEAR(s4.i, 0.0, TOL);
	}

	/* 絶対値・ノルム : |3+4i| = 5 */
	CHECK_NEAR(d_abs(a), 5.0, TOL);
	CHECK_NEAR(d_norm(a), 25.0, TOL);

	/* 偏角 : ラジアンと度が一致すること */
	CHECK_NEAR(d_rad(d_complex(0.0, 1.0)), 2 * atan(1.0), TOL);
	CHECK_NEAR(d_deg(d_complex(0.0, 1.0)), 90.0, 1e-10);
	CHECK_NEAR(d_deg(d_complex(-1.0, 0.0)), 180.0, 1e-10);

	/* 単位複素指数 : |e^{ix}| = 1, e^{i0} = 1 */
	CHECK_NEAR(d_abs(d_exp(0.7)), 1.0, TOL);
	CHECK_NEAR(d_exp(0.0).r, 1.0, TOL);
	CHECK_NEAR(d_exp(0.0).i, 0.0, TOL);

	/* 三角関数 : sin(x+iy) = sin x cosh y + i cos x sinh y */
	{
		const double x = 0.6, y = 0.3;
		const d_complex_t s = d_sin(d_complex(x, y));
		const d_complex_t c = d_cos(d_complex(x, y));
		CHECK_NEAR(s.r, sin(x) * cosh(y), 1e-12);
		CHECK_NEAR(s.i, cos(x) * sinh(y), 1e-12);
		CHECK_NEAR(c.r, cos(x) * cosh(y), 1e-12);
		CHECK_NEAR(c.i, -sin(x) * sinh(y), 1e-12);

		/* 恒等式 sin^2 + cos^2 = 1 */
		{
			const d_complex_t id = d_add(d_mul(s, s), d_mul(c, c));
			CHECK_NEAR(id.r, 1.0, 1e-12);
			CHECK_NEAR(id.i, 0.0, 1e-12);
		}
	}

	/* 平方根 : sqrt(z)^2 == z */
	{
		const d_complex_t r = d_sqrt(a);
		const d_complex_t r2 = d_mul(r, r);
		CHECK_NEAR(r2.r, a.r, 1e-10);
		CHECK_NEAR(r2.i, a.i, 1e-10);
		/* sqrt(-1) = i (主値) */
		{
			const d_complex_t j = d_sqrt(d_complex(-1.0, 0.0));
			CHECK_NEAR(j.r, 0.0, 1e-12);
			CHECK_NEAR(j.i, 1.0, 1e-12);
		}
	}

	/* 逆数 : z * (1/z) == 1 */
	{
		const d_complex_t p = d_mul(a, d_inv(a));
		CHECK_NEAR(p.r, 1.0, 1e-12);
		CHECK_NEAR(p.i, 0.0, 1e-12);
	}

	/* 符号反転 : z + (-z) == 0 */
	{
		const d_complex_t s = d_add(a, d_neg(a));
		CHECK_NEAR(s.r, 0.0, TOL);
		CHECK_NEAR(s.i, 0.0, TOL);
	}

	/* float 版との相互変換で値が保たれること (単精度の丸めの範囲で) */
	{
		const f_complex_t f = f_cast(a);
		const d_complex_t d = d_cast(f);
		CHECK_NEAR(d.r, a.r, 1e-6);
		CHECK_NEAR(d.i, a.i, 1e-6);
	}

	/* float 版の四則も double 版と一致すること */
	{
		const f_complex_t fa = f_complex(3.0, 4.0);
		const f_complex_t fb = f_complex(1.0, -2.0);
		CHECK_NEAR(f_add(fa, fb).r, 4.0, 1e-6);
		CHECK_NEAR(f_sub(fa, fb).i, 6.0, 1e-6);
		CHECK_NEAR(f_mul(fa, fb).r, 11.0, 1e-5);
		CHECK_NEAR(f_mul(fa, fb).i, -2.0, 1e-5);
		CHECK_NEAR(f_rmul(2.0f, fa).r, 6.0, 1e-6);
		CHECK_NEAR(f_add3(fa, fb, f_complex(1.0, 1.0)).r, 5.0, 1e-6);
		CHECK_NEAR(f_add4(fa, fb, f_complex(1.0, 1.0), f_complex(-5.0, -3.0)).r, 0.0, 1e-6);
		/* f_div も (a/b)*b == a で検算する */
		{
			const f_complex_t q = f_div(fa, fb);
			const f_complex_t p = f_mul(q, fb);
			CHECK_NEAR(p.r, 3.0, 1e-5);
			CHECK_NEAR(p.i, 4.0, 1e-5);
		}
	}
}
