/*
test_utils.c

sol/utils.c の汎用ロジックのユニットテスト。

  nearest  … 最近傍の格子点番号
  getspan  … 区間 [p1,p2] を覆う格子点番号の範囲
  calcdft  … 1 周波数の離散フーリエ変換 (Zin / S 行列 / 結合度の土台)
  tokenize … 入力ファイルの行分割

utils.c は ofd.h を include しない自己完結した翻訳単位なので、
グローバル状態の準備なしにそのまま呼べる。
*/

#include "ofd.h"
#include "ofd_test.h"

/* utils.c の関数 (utils.c は独自に同じ形の d_complex_t を typedef している) */
int  tokenize(char *str, const char *tokensep, char *token[], int maxtoken);
int  nearest(double x, int n1, int n2, const double *p);
void getspan(const double p[], int n, int i1, int i2, double p1, double p2, int *n1, int *n2, double eps);
d_complex_t calcdft(int ntime, const double f[], double freq, double dt, double shift);

#define GEPS (1e-6)

static void test_nearest(void)
{
	const double p[8] = {0, 1, 2, 3, 4, 5, 6, 7};

	CHECK(nearest(0.0, 0, 7, p) == 0);
	CHECK(nearest(3.2, 0, 7, p) == 3);
	CHECK(nearest(3.8, 0, 7, p) == 4);
	CHECK(nearest(7.0, 0, 7, p) == 7);

	/* 範囲外は端に張り付く */
	CHECK(nearest(-10.0, 0, 7, p) == 0);
	CHECK(nearest(100.0, 0, 7, p) == 7);

	/* 探索範囲を狭めるとその中の最近傍になる */
	CHECK(nearest(0.0, 3, 5, p) == 3);
	CHECK(nearest(100.0, 3, 5, p) == 5);

	/* 等距離のときは先に見つかった (小さい) 番号を返す
	   (実装は d < dmin の狭義比較なので更新されない) */
	CHECK(nearest(3.5, 0, 7, p) == 3);

	/* 単一点の範囲 */
	CHECK(nearest(99.0, 2, 2, p) == 2);
}

static void test_getspan(void)
{
	const double p[8] = {0, 1, 2, 3, 4, 5, 6, 7};
	int n1, n2;

	/* 格子点に一致する区間 */
	getspan(p, 8, 0, 6, 2.0, 5.0, &n1, &n2, GEPS);
	CHECK(n1 == 2);
	CHECK(n2 == 5);

	/* 格子点の間にある区間は外側に広がる */
	getspan(p, 8, 0, 6, 2.3, 5.7, &n1, &n2, GEPS);
	CHECK(n1 == 2);
	CHECK(n2 == 6);
	/* 不変量 : [p[n1], p[n2]] が [p1, p2] を包含する */
	CHECK(p[n1] - GEPS <= 2.3);
	CHECK(5.7 <= p[n2] + GEPS);

	/* 格子の下端より小さい区間は空 (n1 > n2) を返す */
	getspan(p, 8, 0, 6, -5.0, -3.0, &n1, &n2, GEPS);
	CHECK(n1 > n2);

	/* 格子の上端より大きい区間も空を返す */
	getspan(p, 8, 0, 6, 10.0, 12.0, &n1, &n2, GEPS);
	CHECK(n1 > n2);

	/* i1 == i2 のときはそのまま返す */
	getspan(p, 8, 3, 3, 0.0, 100.0, &n1, &n2, GEPS);
	CHECK(n1 == 3);
	CHECK(n2 == 3);

	/* p1 > p2 で渡しても内部で正規化される */
	{
		int m1, m2;
		getspan(p, 8, 0, 6, 5.0, 2.0, &n1, &n2, GEPS);
		getspan(p, 8, 0, 6, 2.0, 5.0, &m1, &m2, GEPS);
		CHECK(n1 == m1);
		CHECK(n2 == m2);
	}

	/* i1, i2 が範囲外でも [0, n-1] にクランプされる */
	getspan(p, 8, -10, 100, 2.0, 5.0, &n1, &n2, GEPS);
	CHECK(n1 == 2);
	CHECK(n2 == 5);

	/* i1 > i2 で渡しても入れ替えて扱う */
	getspan(p, 8, 6, 0, 2.0, 5.0, &n1, &n2, GEPS);
	CHECK(n1 == 2);
	CHECK(n2 == 5);
}

static void test_calcdft(void)
{
	const double dt = 1.0e-12;
	double f[64];
	int n;

	for (n = 0; n < 64; n++) {
		f[n] = sin(0.3 * n) + 0.5 * cos(0.11 * n);
	}

	/* 周波数 0 では単純な総和 (虚部は 0) */
	{
		double sum = 0;
		d_complex_t c;
		for (n = 0; n < 64; n++) sum += f[n];
		c = calcdft(64, f, 0.0, dt, 0.0);
		CHECK_REL(c.r, sum, 1e-12);
		CHECK_NEAR(c.i, 0.0, 1e-12);
	}

	/* 定義式 sum f[n] * exp(-i w (n+shift) dt) と一致すること */
	{
		const double freq = 3.0e9;
		const double shift = -0.5;
		const double omega = 2 * PI * freq;
		double sr = 0, si = 0;
		d_complex_t c;
		for (n = 0; n < 64; n++) {
			const double ot = omega * (n + shift) * dt;
			sr += cos(ot) * f[n];
			si -= sin(ot) * f[n];
		}
		c = calcdft(64, f, freq, dt, shift);
		CHECK_REL(c.r, sr, 1e-12);
		CHECK_REL(c.i, si, 1e-12);
	}

	/* 線形性 : DFT(2a + 3b) = 2 DFT(a) + 3 DFT(b) */
	{
		const double freq = 5.0e9;
		double a[32], b[32], ab[32];
		d_complex_t ca, cb, cab;
		for (n = 0; n < 32; n++) {
			a[n] = sin(0.2 * n);
			b[n] = cos(0.5 * n);
			ab[n] = 2 * a[n] + 3 * b[n];
		}
		ca = calcdft(32, a, freq, dt, 0.0);
		cb = calcdft(32, b, freq, dt, 0.0);
		cab = calcdft(32, ab, freq, dt, 0.0);
		CHECK_REL(cab.r, 2 * ca.r + 3 * cb.r, 1e-12);
		CHECK_REL(cab.i, 2 * ca.i + 3 * cb.i, 1e-12);
	}

	/* ntime = 0 は 0 を返す */
	{
		const d_complex_t c = calcdft(0, f, 1.0e9, dt, 0.0);
		CHECK_NEAR(c.r, 0.0, 1e-15);
		CHECK_NEAR(c.i, 0.0, 1e-15);
	}

	/* 単一周波数の正弦波はその周波数で最大の応答を持つ */
	{
		const double f0 = 4.0e9;
		double g[256];
		double amp_on, amp_off;
		d_complex_t c;
		for (n = 0; n < 256; n++) {
			g[n] = cos(2 * PI * f0 * n * dt);
		}
		c = calcdft(256, g, f0, dt, 0.0);
		amp_on = sqrt(c.r * c.r + c.i * c.i);
		c = calcdft(256, g, 3 * f0, dt, 0.0);
		amp_off = sqrt(c.r * c.r + c.i * c.i);
		CHECK(amp_on > amp_off);
	}
}

static void test_tokenize(void)
{
	char *token[8];
	char buf[64];
	int n;

	strcpy(buf, "abc def  ghi");
	n = tokenize(buf, " ", token, 8);
	CHECK(n == 3);
	CHECK(!strcmp(token[0], "abc"));
	CHECK(!strcmp(token[1], "def"));
	CHECK(!strcmp(token[2], "ghi"));
	CHECK(token[3] == NULL);   /* 終端が NULL であること */

	/* 区切り文字だけの行はトークン 0 個 */
	strcpy(buf, "   ");
	n = tokenize(buf, " ", token, 8);
	CHECK(n == 0);
	CHECK(token[0] == NULL);

	/* maxtoken で打ち切られること */
	strcpy(buf, "a b c d e");
	n = tokenize(buf, " ", token, 3);
	CHECK(n == 3);

	/* NULL / maxtoken=0 は 0 を返す (防御的な入口条件) */
	CHECK(tokenize(NULL, " ", token, 8) == 0);
	strcpy(buf, "a b");
	CHECK(tokenize(buf, " ", token, 0) == 0);

	/* 複数の区切り文字 (空白とタブ) */
	strcpy(buf, "x\ty z");
	n = tokenize(buf, " \t", token, 8);
	CHECK(n == 3);
	CHECK(!strcmp(token[1], "y"));
}

void test_utils(void)
{
	test_nearest();
	test_getspan();
	test_calcdft();
	test_tokenize();
}
