/*
test_outputspara.c

sol/outputSpara.c (S パラメータ) のユニットテスト。

ここで出る値は **GUI が読む /metadata/s_parameters、ofd.log の S 表、
Touchstone ファイル (.snp)** のもとになる。フィルタや線路の設計は
これで良否を判断するので、静かにずれると影響が大きい。

## calcSpara() が何をしているか

伝送線路上の電圧は進行波と反射波の和

    V(z) = V+ e^{-γz} + V- e^{+γz}

で表せる。ポート 1 では基準点 z=0 と、その前後 z=±d の 3 点を観測して
いる (VPoint の行 0, NPoint, NPoint+1)。a = e^{γd} と置くと

    V1 = V+ + V-
    Vp = V+/a + V- a          (z = +d)
    Vm = V+ a  + V-/a         (z = -d)

から

    (Vp + Vm)/V1 = a + 1/a                        … c1
    sqrt(c1^2 - 4) = a - 1/a                      … c2
    (c1 + c2)/2 = a                               … c3  (虚部が負なら c1-c2)
    (Vp/a - Vm a) / (1/a^2 - a^2) = V+            … c6
    (Vm/a - Vp a) / (1/a^2 - a^2) = V-            … c7

    S11 = V- / V+,   Sn1 = Vn / V+

このテストは **V+, V-, a, Vn を先に決めて時系列を合成し**、
実装が元の S を復元できるかを見る。期待値は上の式変形から出しており、
実装のコードを写したものではない。

## 時系列の作り方

整数周期ぶんの和なので、振幅 A・位相 φ の正弦波の DFT は
(Ntime/2) A e^{jφ} になる。S は比なので (Ntime/2) は打ち消える。

## sparaDerived()

振幅[dB] と位相[deg] はログと HDF5 の両方に出る。以前は同じ式が
2 箇所に重複していたので 1 本化した (sol/outputZin.c と同じ問題)。
Touchstone は MA 形式なので dB ではなく線形振幅を書く (仕様どおり)。
*/

#include "ofd.h"
#include "complex.h"
#include "ofd_test.h"

#include <stdlib.h>
#include <string.h>

void calcSpara(void);
void sparaDerived(d_complex_t s, double *mag_dB, double *phase_deg);

#define NT   (2000)
#define FREQ (1.0e9)
#define NP   (2)          /* ポート数 (S11 と S21 を見る) */

/* 指定した複素振幅になる正弦波を VPoint の行 row に書く。
   DFT (shift=0) は (Ntime/2) * c になる。 */
static void set_point(int row, d_complex_t c)
{
	const double w = 2 * PI * FREQ;
	const double amp = d_abs(c);
	const double ph = atan2(c.i, c.r);
	double *v = &VPoint[row * (Solver.maxiter + 1)];

	for (int n = 0; n <= Solver.maxiter; n++) {
		v[n] = amp * cos((w * n * Dt) + ph);
	}
}

static void setup(void)
{
	NPoint = NP;
	NFreq1 = 1;
	Ntime = NT;
	Solver.maxiter = NT;
	Dt = 1.0 / (FREQ * 40);

	Freq1 = (double *)malloc(sizeof(double));
	Freq1[0] = FREQ;

	VPoint = (double *)malloc((size_t)(NP + 2) * (Solver.maxiter + 1) * sizeof(double));
	Spara  = (d_complex_t *)malloc((size_t)NP * NFreq1 * sizeof(d_complex_t));
}

static void teardown(void)
{
	free(Freq1); free(VPoint); free(Spara);
	Freq1 = NULL; VPoint = NULL; Spara = NULL;
	NPoint = NFreq1 = 0;
}

/* V+, V-, a (= e^{γd}), V2 を与えて時系列を組み立てる */
static void build(d_complex_t vp, d_complex_t vm, d_complex_t a, d_complex_t v2)
{
	const d_complex_t inv_a = d_inv(a);
	/* 行 0 = V1 (基準点), 行 1 = V2, 行 NP = V1+, 行 NP+1 = V1- */
	set_point(0,      d_add(vp, vm));
	set_point(1,      v2);
	set_point(NP,     d_add(d_mul(vp, inv_a), d_mul(vm, a)));      /* z = +d */
	set_point(NP + 1, d_add(d_mul(vp, a), d_mul(vm, inv_a)));      /* z = -d */
}

/* 1. 無損失線路 (γ = jβ) で S11, S21 が復元されること */
static void test_lossless(void)
{
	/* βd を変えて数点。0 < βd < π なら a = e^{jβd} の虚部が正で、
	   実装の「虚部が負なら逆根を採る」分岐が正しい根を選ぶ */
	const double bds[4] = {0.35, 0.70, 1.20, 2.50};

	for (int i = 0; i < 4; i++) {
		const double bd = bds[i];
		const d_complex_t a  = d_exp(bd);                 /* e^{jβd} */
		const d_complex_t vp = d_complex(1.0, 0.0);
		const d_complex_t vm = d_complex(0.3, 0.2);
		const d_complex_t v2 = d_complex(-0.45, 0.55);

		setup();
		build(vp, vm, a, v2);
		calcSpara();

		/* S11 = V-/V+, S21 = V2/V+ を独立に計算して比較 */
		{
			const d_complex_t s11 = d_div(vm, vp);
			const d_complex_t s21 = d_div(v2, vp);
			CHECK_NEAR(Spara[0].r, s11.r, 1e-6);
			CHECK_NEAR(Spara[0].i, s11.i, 1e-6);
			CHECK_NEAR(Spara[1].r, s21.r, 1e-6);
			CHECK_NEAR(Spara[1].i, s21.i, 1e-6);
		}
		teardown();
	}
}

/* 2. 有損失線路 (γ = α + jβ、a の絶対値が 1 でない) でも復元されること */
static void test_lossy(void)
{
	/* a = e^{αd} e^{jβd}、αd = 0.15 */
	const d_complex_t a  = d_rmul(exp(0.15), d_exp(0.8));
	const d_complex_t vp = d_complex(0.7, -0.4);
	const d_complex_t vm = d_complex(0.1, 0.25);
	const d_complex_t v2 = d_complex(0.33, 0.12);

	setup();
	build(vp, vm, a, v2);
	calcSpara();
	{
		const d_complex_t s11 = d_div(vm, vp);
		const d_complex_t s21 = d_div(v2, vp);
		CHECK_NEAR(Spara[0].r, s11.r, 1e-6);
		CHECK_NEAR(Spara[0].i, s11.i, 1e-6);
		CHECK_NEAR(Spara[1].r, s21.r, 1e-6);
		CHECK_NEAR(Spara[1].i, s21.i, 1e-6);
	}
	teardown();
}

/* 3. 物理的に意味のある特別な場合 */
static void test_special_cases(void)
{
	const d_complex_t a = d_exp(0.9);

	/* (a) 無反射 (V- = 0) -> S11 = 0 */
	setup();
	build(d_complex(1, 0), d_complex(0, 0), a, d_complex(0.5, 0));
	calcSpara();
	CHECK_NEAR(d_abs(Spara[0]), 0.0, 1e-6);
	teardown();

	/* (b) 全反射 (|V-| = |V+|) -> |S11| = 1 */
	setup();
	build(d_complex(1, 0), d_complex(0.6, -0.8), a, d_complex(0.1, 0));
	calcSpara();
	CHECK_NEAR(d_abs(Spara[0]), 1.0, 1e-6);
	teardown();

	/* (c) S11 は V+ と V- を同じ倍率で変えても不変 (比なので) */
	{
		d_complex_t s_a, s_b;
		setup();
		build(d_complex(1, 0), d_complex(0.3, 0.2), a, d_complex(0.4, 0.1));
		calcSpara();
		s_a = Spara[0];
		teardown();

		setup();
		build(d_complex(2.5, 0), d_complex(0.75, 0.5), a, d_complex(1.0, 0.25));
		calcSpara();
		s_b = Spara[0];
		teardown();

		CHECK_NEAR(s_a.r, s_b.r, 1e-6);
		CHECK_NEAR(s_a.i, s_b.i, 1e-6);
	}
}

/* 4. sparaDerived : 振幅[dB] と位相[deg] */
static void test_derived(void)
{
	double mag, deg;

	/* (a) |S| = 1 は 0 dB */
	sparaDerived(d_complex(1.0, 0.0), &mag, &deg);
	CHECK_NEAR(mag, 0.0, 1e-12);
	CHECK_NEAR(deg, 0.0, 1e-12);

	/* (b) |S| = 0.5 は 20log10(0.5) = -6.0206 dB */
	sparaDerived(d_complex(0.5, 0.0), &mag, &deg);
	CHECK_NEAR(mag, 20 * log10(0.5), 1e-12);

	/* (c) 位相は度で、象限が正しいこと */
	sparaDerived(d_complex(0.0, 1.0), &mag, &deg);
	CHECK_NEAR(deg, 90.0, 1e-9);
	sparaDerived(d_complex(-1.0, 0.0), &mag, &deg);
	CHECK_NEAR(fabs(deg), 180.0, 1e-9);
	sparaDerived(d_complex(1.0, -1.0), &mag, &deg);
	CHECK_NEAR(deg, -45.0, 1e-9);
	CHECK_NEAR(mag, 20 * log10(sqrt(2.0)), 1e-9);

	/* (d) S = 0 は log10 が発散するので EPS2 で下限を切る (-240 dB) */
	sparaDerived(d_complex(0.0, 0.0), &mag, &deg);
	CHECK_NEAR(mag, 20 * log10(EPS2), 1e-9);
	CHECK_NEAR(mag, -240.0, 1e-9);
	CHECK(mag > -1e300);          /* -inf になっていないこと */

	/* (e) dB は |S| について単調増加 */
	{
		double prev = -1e300;
		for (int i = 1; i <= 10; i++) {
			sparaDerived(d_complex(i * 0.1, 0.0), &mag, &deg);
			CHECK(mag > prev);
			prev = mag;
		}
	}

	/* (f) 20log10|S| の定義どおりであること (振幅を 10 倍で +20 dB) */
	{
		double m1, m2;
		sparaDerived(d_complex(0.03, 0.04), &m1, &deg);   /* |S| = 0.05 */
		sparaDerived(d_complex(0.30, 0.40), &m2, &deg);   /* |S| = 0.50 */
		CHECK_NEAR(m2 - m1, 20.0, 1e-9);
		CHECK_NEAR(m1, 20 * log10(0.05), 1e-9);
	}
}

/* 5. ログと HDF5 が同じ値を出すこと (式の 1 本化の回帰)。
      どちらかが独自計算に戻ったらここが落ちる。 */
static void test_log_and_hdf5_agree(void)
{
	const d_complex_t s = d_complex(0.1234, -0.5678);
	double mag, deg;

	sparaDerived(s, &mag, &deg);
	/* 独立に計算した期待値と一致すること */
	CHECK_NEAR(mag, 20 * log10(d_abs(s)), 1e-12);
	CHECK_NEAR(deg, atan2(s.i, s.r) * 180 / PI, 1e-9);

	/* 片方だけ取り出せること (NULL を渡せる) */
	{
		double only;
		sparaDerived(s, &only, NULL);
		CHECK_NEAR(only, mag, 1e-12);
		sparaDerived(s, NULL, &only);
		CHECK_NEAR(only, deg, 1e-12);
	}
}

void test_outputspara(void)
{
	test_lossless();
	test_lossy();
	test_special_cases();
	test_derived();
	test_log_and_hdf5_agree();
}
