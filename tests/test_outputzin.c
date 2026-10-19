/*
test_outputzin.c

sol/outputZin.c (入力インピーダンスと反射・入力電力) のユニットテスト。

ここで出る値は **GUI が読む /metadata/input_impedance と ofd.log の
インピーダンス表そのもの**で、アンテナ・線路の良否をこれで判断する。
壊れても実行は normal end で終わるので、値がおかしいことにしか
現れない。

検証する量と、それぞれの独立な根拠:

  Zin  = DFT(V) / DFT(I)
         給電点の電圧と電流の比。既知の V(t), I(t) を与えれば
         答えは解析的に決まる。電流の DFT には -0.5 のシフトが
         入る (leapfrog で E と H が半ステップずれるため) ので、
         **半ステップずらした電流を与えたときに位相が戻ること**で
         このシフトを固定する。
  Ref  = 10 log10 |Γ|^2 = 20 log10 |Γ|,  Γ = (Zin - Z0) / (Zin + Z0)
         整合時 (Zin = Z0) は Γ = 0、開放・短絡では |Γ| = 1 (0 dB)。
  Pin  = Re(V conj(I))
         Pin[0] が入力電力、Pin[1] は不整合損を戻した値で
         Pin[1] = Pin[0] / (1 - |Γ|^2)。整合時は両者が一致する。
  Y    = 1/Zin  (mS 単位なので 1e3 倍)
  VSWR = (1 + |Γ|) / (1 - |Γ|)

VSWR と Y は zinDerived() にまとめてある。以前は同じ式が HDF5 書き出しと
ログ出力の 2 箇所に重複していて、片方だけ直すと GUI とログが食い違う
状態だった。ここではその 1 本化した関数を検証する。

期待値は**実装式の写しではなく**、上の定義から独立に計算している
(.claude/rules/features.md の方針)。
*/

#include "ofd.h"
#include "complex.h"
#include "ofd_test.h"

#include <stdlib.h>
#include <string.h>

void calcZin(void);
void calcPin(void);
void zinDerived(d_complex_t zin, double ref_db, double *gin_mS, double *bin_mS, double *vswr);

#define NT   (2000)      /* 時間ステップ数 */
#define FREQ (1.0e9)     /* 試験周波数 */

/* 給電の時系列を作る。
   V(t) = cos(ω t)、I(t) = |Y| cos(ω t + φ) とすると Zin = 1/(|Y| e^{jφ})。

   電流は leapfrog により**半ステップ前**の時刻に置く。calcdft() は
   Σ f[n] exp(-jω(n + shift)Δt) を計算し、実装は電流に shift = -0.5 を
   渡すので、サンプル I[n] は時刻 (n - 0.5)Δt のものとして扱われる。
   ここで同じ時刻に値を置いておけば、Zin は与えた負荷そのものになる。 */
#define VPHASE (37.0 * PI / 180)   /* 電圧の初期位相 */

static void set_feed(double zr, double zi)
{
	const double w = 2 * PI * FREQ;
	const d_complex_t z = d_complex(zr, zi);
	const d_complex_t y = d_inv(z);
	const double ya = d_abs(y);
	const double yp = atan2(y.i, y.r);

	/* 電圧に 0 でない位相を与えるのは必須。位相 0 だと DFT(V) が実数に
	   なり Im(V) = 0 なので、Pin = Re(V conj(I)) = Re(V)Re(I) + Im(V)Im(I)
	   の**第 2 項が常に 0** になる。実際その状態では、第 2 項を落とす変異も
	   符号を反転する変異も検出できなかった。 */
	for (int n = 0; n <= Solver.maxiter; n++) {
		const double t = n * Dt;
		VFeed[n] = cos((w * t) + VPHASE);
		/* 電流は t - Dt/2 の時刻の値 (calcdft の shift = -0.5 に対応) */
		IFeed[n] = ya * cos((w * (t - (0.5 * Dt))) + VPHASE + yp);
	}
}

/* 整数周期ぶんの和なので、DFT は解析的に決まる:
     DFT(V) = (Ntime/2) e^{jθ},  DFT(I) = (Ntime/2) e^{jθ} Y
   したがって
     Pin = Re(V conj(I)) = (Ntime/2)^2 * Re(Y) = (Ntime/2)^2 * G
   電圧の位相 θ は打ち消えるので、期待値は Z だけで決まる。
   calcPin は V と I の両方を正規化係数 C で割るので、Pin は |C|^2 で
   割られる (C の位相も打ち消える)。 */
static double expected_pin(double zr, double zi)
{
	const d_complex_t y = d_inv(d_complex(zr, zi));
	const double half = NT / 2.0;
	return half * half * y.r / d_norm(cFdft[0]);
}

static void setup(void)
{
	NFeed = 1;
	NFreq1 = 1;
	NFreq2 = 1;
	Ntime = NT;
	Solver.maxiter = NT;
	Dt = 1.0 / (FREQ * 40);          /* 1 周期あたり 40 点 */

	Freq1 = (double *)malloc(sizeof(double));
	Freq2 = (double *)malloc(sizeof(double));
	Freq1[0] = Freq2[0] = FREQ;

	Feed = (feed_t *)calloc(1, sizeof(feed_t));
	Feed[0].z0 = 50.0;

	VFeed = (double *)malloc((size_t)(Solver.maxiter + 1) * sizeof(double));
	IFeed = (double *)malloc((size_t)(Solver.maxiter + 1) * sizeof(double));

	/* 入射スペクトルによる正規化 (calcPin が V と I の両方を割る)。
	   **1 にしてはいけない** — 1 だと割り算が恒等変換になり、正規化を
	   丸ごと外す変異を検出できない。|C| != 1 の複素数にしておくと
	   Pin が |C|^2 で割られるので、外れたことが値に出る。 */
	cFdft = (d_complex_t *)malloc(sizeof(d_complex_t));
	cFdft[0] = d_complex(2, 1);      /* |C|^2 = 5 */

	Zin = NULL; Ref = NULL; Pin[0] = Pin[1] = NULL;
}

static void teardown(void)
{
	free(Freq1); free(Freq2); free(Feed);
	free(VFeed); free(IFeed); free(cFdft);
	free(Zin); free(Ref); free(Pin[0]); free(Pin[1]);
	Zin = NULL; Ref = NULL; Pin[0] = Pin[1] = NULL;
	NFeed = NFreq1 = NFreq2 = 0;
}

/* 1. Zin が与えた負荷と一致すること */
static void test_zin(void)
{
	/* 抵抗・誘導性・容量性・整合の 4 通り */
	const double zs[4][2] = {{75.0, 0.0}, {30.0, 40.0}, {60.0, -25.0}, {50.0, 0.0}};

	for (int i = 0; i < 4; i++) {
		setup();
		set_feed(zs[i][0], zs[i][1]);
		calcZin();

		/* DFT は有限長なので厳密ではない。相対 1e-3 で見る */
		CHECK_REL(Zin[0].r, zs[i][0], 1e-3);
		CHECK_REL(Zin[0].i, zs[i][1], 1e-3);

		/* Ref[dB] = 20 log10 |Γ|, Γ = (Z - Z0)/(Z + Z0) を独立に計算 */
		{
			const d_complex_t z = d_complex(zs[i][0], zs[i][1]);
			const d_complex_t z0 = d_complex(50.0, 0);
			const d_complex_t g = d_div(d_sub(z, z0), d_add(z, z0));
			const double ga = d_abs(g);
			if (ga > 1e-12) {
				CHECK_NEAR(Ref[0], 20 * log10(ga), 1e-2);
			}
			else {
				CHECK(Ref[0] < -100);   /* 整合時は限りなく小さい */
			}
		}
		teardown();
	}
}

/* 2. 電流の半ステップシフトが効いていること。
      実装は電流の DFT に -0.5 のシフトを掛けて E と H の時刻を揃える。
      これが無い/符号が逆だと、純抵抗を与えても Zin に虚部が出る。 */
static void test_half_step_shift(void)
{
	setup();
	set_feed(75.0, 0.0);      /* 純抵抗 */
	calcZin();
	/* 半ステップの位相は ω Dt/2 = 2π/80 = 4.5°。シフトが無ければ
	   虚部が |Z| tan(4.5°) ≈ 5.9 出るので、0.1 で見れば十分に区別できる */
	CHECK_NEAR(Zin[0].i, 0.0, 0.1);
	CHECK(fabs(Zin[0].i) < 1.0);
	teardown();
}

/* 3. 入力電力 Pin */
static void test_pin(void)
{
	/* (a) 整合負荷 : 不整合損の補正が入らない (Pin[1] == Pin[0]) */
	setup();
	set_feed(50.0, 0.0);
	calcPin();
	CHECK(Pin[0][0] > 0);
	CHECK_REL(Pin[1][0], Pin[0][0], 1e-6);
	/* 絶対値も見る。Z = 50 なら G = 20 mS、|C|^2 = 5 なので
	   Pin = 1000^2 * 0.02 / 5 = 4000 */
	CHECK_REL(Pin[0][0], expected_pin(50.0, 0.0), 1e-3);
	CHECK_REL(Pin[0][0], 4000.0, 1e-3);
	teardown();

	/* (a2) 複素負荷でも Pin = (Ntime/2)^2 * Re(1/Z) になること。
	   Re(V conj(I)) の 2 項が両方効く配置 (電圧位相が 0 でない) */
	setup();
	set_feed(30.0, 40.0);
	calcPin();
	CHECK_REL(Pin[0][0], expected_pin(30.0, 40.0), 1e-3);
	teardown();

	/* (b) 不整合負荷 : Pin[1] = Pin[0] / (1 - |Γ|^2) */
	setup();
	set_feed(150.0, 0.0);
	calcPin();
	{
		const d_complex_t z = d_complex(150.0, 0);
		const d_complex_t z0 = d_complex(50.0, 0);
		const d_complex_t g = d_div(d_sub(z, z0), d_add(z, z0));
		const double gn = d_norm(g);          /* |Γ|^2 = 0.25 */
		CHECK_NEAR(gn, 0.25, 1e-9);
		CHECK_REL(Pin[1][0], Pin[0][0] / (1 - gn), 1e-3);
		CHECK(Pin[1][0] > Pin[0][0]);         /* 補正で必ず増える */
		CHECK_REL(Pin[0][0], expected_pin(150.0, 0.0), 1e-3);
	}
	teardown();
}

/* 4. zinDerived : アドミタンスと VSWR */
static void test_derived(void)
{
	double gin, bin, vswr;

	/* (a) 純抵抗 50Ω、整合 (Ref = -inf 相当の十分小さい値) */
	zinDerived(d_complex(50.0, 0.0), -300.0, &gin, &bin, &vswr);
	CHECK_NEAR(gin, 1000.0 / 50.0, 1e-9);   /* Y = 1/50 S = 20 mS */
	CHECK_NEAR(bin, 0.0, 1e-9);
	CHECK_NEAR(vswr, 1.0, 1e-6);            /* 整合なら VSWR = 1 */

	/* (b) 複素インピーダンス : Y = 1/Z を独立に計算して比較 */
	{
		const d_complex_t z = d_complex(30.0, 40.0);   /* |Z| = 50 */
		zinDerived(z, -6.0, &gin, &bin, &vswr);
		/* 1/(30+40j) = (30-40j)/2500 = 0.012 - 0.016j  S */
		CHECK_NEAR(gin, 12.0, 1e-9);
		CHECK_NEAR(bin, -16.0, 1e-9);
	}

	/* (c) VSWR が |Γ| = 10^(Ref/20) から作られること。
	   Ref = -20 dB -> |Γ| = 0.1 -> VSWR = 1.1/0.9 = 1.2222... */
	zinDerived(d_complex(50.0, 0.0), -20.0, &gin, &bin, &vswr);
	CHECK_NEAR(vswr, 1.1 / 0.9, 1e-9);
	/* Ref = -6.0206 dB -> |Γ| = 0.5 -> VSWR = 3 */
	zinDerived(d_complex(50.0, 0.0), 20 * log10(0.5), &gin, &bin, &vswr);
	CHECK_NEAR(vswr, 3.0, 1e-9);

	/* (d) 全反射 (|Γ| = 1, Ref = 0 dB) では発散するので頭打ちにすること */
	zinDerived(d_complex(1e9, 0.0), 0.0, &gin, &bin, &vswr);
	CHECK_NEAR(vswr, 1000.0, 1e-9);

	/* (e) VSWR は |Γ| について単調増加で、常に 1 以上 */
	{
		double prev = 0;
		for (int i = 0; i <= 9; i++) {
			const double ga = i * 0.1;               /* |Γ| = 0.0 .. 0.9 */
			const double ref = (ga > 0) ? 20 * log10(ga) : -300.0;
			zinDerived(d_complex(50.0, 0.0), ref, &gin, &bin, &vswr);
			CHECK(vswr >= 1.0 - 1e-9);
			CHECK(vswr >= prev - 1e-9);
			prev = vswr;
		}
	}
}

/* 5. ログと HDF5 が同じ値を出すこと (式の 1 本化の回帰)。
      zinDerived() を経由していれば同じなのは自明だが、どちらかが
      独自計算に戻ったらここが落ちる。GUI とログの食い違いは
      利用者からは「どちらが正しいのか分からない」形で現れる。 */
static void test_log_and_hdf5_agree(void)
{
	const d_complex_t z = d_complex(46.736, 1.365);   /* 実測値 (45 GHz のマイクロストリップ) */
	const d_complex_t z0 = d_complex(50.0, 0);
	const d_complex_t g = d_div(d_sub(z, z0), d_add(z, z0));
	const double ref = 20 * log10(d_abs(g));

	double gin, bin, vswr;
	zinDerived(z, ref, &gin, &bin, &vswr);

	/* 独立に計算した期待値と一致すること
	   (この 4 値は ofd.log の表と HDF5 の input_impedance の両方に出る) */
	{
		const d_complex_t y = d_inv(z);
		CHECK_NEAR(gin, y.r * 1e3, 1e-9);
		CHECK_NEAR(bin, y.i * 1e3, 1e-9);
		const double ga = d_abs(g);
		CHECK_NEAR(vswr, (1 + ga) / (1 - ga), 1e-9);
	}
	/* 実測ログの値と桁が合うこと (Gin 21.379 mS, Bin -0.624 mS, VSWR 1.076) */
	CHECK_NEAR(gin, 21.379, 1e-3);
	CHECK_NEAR(bin, -0.624, 1e-3);
	CHECK_NEAR(vswr, 1.076, 1e-3);
	CHECK_NEAR(ref, -28.737, 1e-3);
}

void test_outputzin(void)
{
	test_zin();
	test_half_step_shift();
	test_pin();
	test_derived();
	test_log_and_hdf5_agree();
}
