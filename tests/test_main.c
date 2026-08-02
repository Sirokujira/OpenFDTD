/*
test_main.c

C 実装ユニットテストのランナー。

このファイルだけが MAIN を定義して ofd.h のグローバル変数の実体を作る
(src/sol_Main.c と同じ方式)。finc.h など一部のロジックが WaveAmp /
WaveOmega といったグローバルを読むため、それらの置き場所が必要になる。

実行方法:
    ofd_tests            … 全テストを実行する
    ofd_tests <name>     … 指定した 1 件だけ実行する (CTest はこの形で呼ぶ)
    ofd_tests --list     … テスト名を一覧表示する
*/

#define MAIN
#include "ofd.h"
#undef MAIN

#include "ofd_test.h"

#include <string.h>

/* 集計カウンタ (ofd_test.h の extern 宣言の実体) */
int ofd_test_checks = 0;
int ofd_test_failures = 0;

void ofd_test_fail(const char *expr, const char *file, int line)
{
	ofd_test_failures++;
	fprintf(stderr, "  FAIL %s:%d : %s\n", file, line, expr);
}

/* 各テストファイルが提供する関数 */
void test_complex(void);
void test_finc(void);
void test_utils(void);
void test_vfeed(void);
void test_ingeometry(void);
void test_geomlines(void);

typedef struct {
	const char *name;
	void (*func)(void);
} testcase_t;

static const testcase_t tests[] = {
	{"complex",    test_complex},
	{"finc",       test_finc},
	{"utils",      test_utils},
	{"vfeed",      test_vfeed},
	{"ingeometry", test_ingeometry},
	{"geomlines",  test_geomlines},
};

static const int ntest = (int)(sizeof(tests) / sizeof(tests[0]));

static int run(const testcase_t *t)
{
	const int checks0 = ofd_test_checks;
	const int fails0 = ofd_test_failures;

	printf("[ RUN      ] %s\n", t->name);
	t->func();

	const int nc = ofd_test_checks - checks0;
	const int nf = ofd_test_failures - fails0;
	if (nf == 0) {
		printf("[       OK ] %s (%d checks)\n", t->name, nc);
	}
	else {
		printf("[  FAILED  ] %s (%d/%d checks failed)\n", t->name, nf, nc);
	}

	return nf;
}

int main(int argc, char *argv[])
{
	int i;

	if ((argc > 1) && !strcmp(argv[1], "--list")) {
		for (i = 0; i < ntest; i++) {
			printf("%s\n", tests[i].name);
		}
		return 0;
	}

	if (argc > 1) {
		for (i = 0; i < ntest; i++) {
			if (!strcmp(argv[1], tests[i].name)) {
				return run(&tests[i]) ? 1 : 0;
			}
		}
		fprintf(stderr, "*** unknown test : %s\n", argv[1]);
		return 2;
	}

	for (i = 0; i < ntest; i++) {
		run(&tests[i]);
	}

	printf("\n%d checks, %d failures\n", ofd_test_checks, ofd_test_failures);

	return ofd_test_failures ? 1 : 0;
}
