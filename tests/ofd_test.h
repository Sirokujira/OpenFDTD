/*
ofd_test.h

C 実装のユニットテスト用の最小アサートフレームワーク。

外部依存を持たない (MSVC / gcc / clang の 3 OS で同じように動く)。
テストは失敗しても即座に abort せず、全件走らせてから件数を集計する
(1 回の CI 実行でできるだけ多くの情報を出すため)。

使い方:
    #include "ofd_test.h"

    void test_something(void)
    {
        CHECK(1 + 1 == 2);
        CHECK_NEAR(sqrt(2.0), 1.41421356, 1e-8);
    }

テスト関数は tests/test_main.c のテーブルに登録する。
*/

#ifndef OFD_TEST_H
#define OFD_TEST_H

#include <stdio.h>
#include <math.h>

/* 集計用のカウンタ (実体は test_main.c) */
extern int ofd_test_checks;
extern int ofd_test_failures;

void ofd_test_fail(const char *expr, const char *file, int line);

/* 条件が偽なら失敗として記録する */
#define CHECK(cond) \
	do { \
		ofd_test_checks++; \
		if (!(cond)) ofd_test_fail(#cond, __FILE__, __LINE__); \
	} while (0)

/* |a - b| <= tol を検査する (NaN は必ず失敗する) */
#define CHECK_NEAR(a, b, tol) \
	do { \
		const double ofd_a_ = (double)(a); \
		const double ofd_b_ = (double)(b); \
		const double ofd_t_ = (double)(tol); \
		ofd_test_checks++; \
		if (!(fabs(ofd_a_ - ofd_b_) <= ofd_t_)) { \
			fprintf(stderr, "    got=%.17g want=%.17g tol=%.17g\n", \
				ofd_a_, ofd_b_, ofd_t_); \
			ofd_test_fail(#a " ~= " #b, __FILE__, __LINE__); \
		} \
	} while (0)

/* 相対誤差で比較する (値の桁が大きい/未知のとき) */
#define CHECK_REL(a, b, rel) \
	do { \
		const double ofd_a_ = (double)(a); \
		const double ofd_b_ = (double)(b); \
		const double ofd_r_ = (double)(rel); \
		const double ofd_s_ = fabs(ofd_b_) > 1 ? fabs(ofd_b_) : 1.0; \
		ofd_test_checks++; \
		if (!(fabs(ofd_a_ - ofd_b_) <= ofd_r_ * ofd_s_)) { \
			fprintf(stderr, "    got=%.17g want=%.17g rel=%.17g\n", \
				ofd_a_, ofd_b_, ofd_r_); \
			ofd_test_fail(#a " ~= " #b, __FILE__, __LINE__); \
		} \
	} while (0)

#endif  /* OFD_TEST_H */
