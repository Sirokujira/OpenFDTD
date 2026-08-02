/*
test_geomlines.c

sol/geomlines.c の 3D 形状ワイヤフレーム生成のユニットテスト。

geomlines() は mode=0 で必要な線分数を返し、mode=1 でその配列を埋める
2 段階の API になっている。呼び出し側は mode=0 の戻り値で配列を確保するので、
「mode=1 が mode=0 の見積もりを超えて書き込まない」ことが
バッファオーバーランを防ぐ最重要の不変量になる。
*/

#include "ofd.h"
#include "ofd_test.h"

int geomlines(int mode, int ngeometry, int *shape, int *id,
	double (*geometry)[8], double (*gline)[2][3], int *mline, double eps);

#define GEPS (1e-6)

/* mode=0 の見積もりで配列を確保し mode=1 を実行して共通の不変量を検査する */
static void check_shape(int shape_code, double g[8], int matid)
{
	int shape[1];
	int id[1];
	double geometry[1][8];
	double (*gline)[2][3];
	int *mline;
	int n0, n1, i, j;

	shape[0] = shape_code;
	id[0] = matid;
	for (i = 0; i < 8; i++) geometry[0][i] = g[i];

	/* mode=0 : 線分数の見積もり */
	n0 = geomlines(0, 1, shape, id, geometry, NULL, NULL, GEPS);
	CHECK(n0 > 0);
	if (n0 <= 0) return;

	/* MSVC は VLA を持たないので malloc で確保する */
	gline = (double (*)[2][3])malloc((size_t)n0 * sizeof(*gline));
	mline = (int *)malloc((size_t)n0 * sizeof(*mline));
	CHECK(gline != NULL);
	CHECK(mline != NULL);
	if ((gline == NULL) || (mline == NULL)) {
		free(gline);
		free(mline);
		return;
	}
	for (i = 0; i < n0; i++) mline[i] = -1;

	/* mode=1 : 実際の線分データ */
	n1 = geomlines(1, 1, shape, id, geometry, gline, mline, GEPS);

	/* 最重要 : 見積もりを超えて書き込まないこと */
	CHECK(n1 <= n0);
	CHECK(n1 > 0);

	/* 書き込まれた線分はすべて指定した物性値番号を持つこと */
	for (i = 0; i < n1; i++) {
		CHECK(mline[i] == matid);
	}

	/* 線分の端点がすべて有限であること (NaN/Inf を作らない) */
	for (i = 0; i < n1; i++) {
		for (j = 0; j < 3; j++) {
			CHECK(gline[i][0][j] == gline[i][0][j]);   /* NaN でない */
			CHECK(gline[i][1][j] == gline[i][1][j]);
		}
	}

	free(gline);
	free(mline);
}

void test_geomlines(void)
{
	/* 直方体 : 12 稜線 */
	{
		double g[8] = {0, 1, 0, 2, 0, 3, 0, 0};
		int shape[1];
		int id[1];
		double geometry[1][8];
		int i;
		shape[0] = 1;
		id[0] = 5;
		for (i = 0; i < 8; i++) geometry[0][i] = g[i];
		CHECK(geomlines(0, 1, shape, id, geometry, NULL, NULL, GEPS) == 12);
		check_shape(1, g, 5);
	}

	/* 楕円体 : 3 断面 x rdiv(=72) */
	{
		double g[8] = {-1, 1, -1, 1, -1, 1, 0, 0};
		int shape[1];
		int id[1];
		double geometry[1][8];
		int i;
		shape[0] = 2;
		id[0] = 2;
		for (i = 0; i < 8; i++) geometry[0][i] = g[i];
		CHECK(geomlines(0, 1, shape, id, geometry, NULL, NULL, GEPS) == 3 * 72);
		check_shape(2, g, 2);
	}

	/* 円柱 (X/Y/Z) : 2 x rdiv + 4 */
	{
		double gx[8] = {0, 2, -1, 1, -1, 1, 0, 0};
		double gy[8] = {-1, 1, 0, 2, -1, 1, 0, 0};
		double gz[8] = {-1, 1, -1, 1, 0, 2, 0, 0};
		check_shape(11, gx, 3);
		check_shape(12, gy, 3);
		check_shape(13, gz, 3);
	}

	/* 三角柱 : 9 線分 */
	{
		double g[8] = {0, 1, 0, 1, 0, 0, 0, 1};
		check_shape(31, g, 4);
		check_shape(32, g, 4);
		check_shape(33, g, 4);
	}

	/* 角錐台 : 12 線分 */
	{
		double g[8] = {0, 1, 0, 0, 2, 2, 1, 1};
		check_shape(41, g, 6);
		check_shape(42, g, 6);
		check_shape(43, g, 6);
	}

	/* 円錐台 : 2 x rdiv + 4 */
	{
		double g[8] = {0, 1, 0, 0, 2, 2, 1, 1};
		check_shape(51, g, 7);
		check_shape(52, g, 7);
		check_shape(53, g, 7);
	}

	/* 未知の形状コードは線分を生成しない */
	{
		int shape[1];
		int id[1];
		double geometry[1][8];
		int i;
		shape[0] = 999;
		id[0] = 1;
		for (i = 0; i < 8; i++) geometry[0][i] = 0;
		CHECK(geomlines(0, 1, shape, id, geometry, NULL, NULL, GEPS) == 0);
	}

	/* 複数形状の合計 : 直方体 2 個で 24 */
	{
		int shape[2];
		int id[2];
		double geometry[2][8];
		int i;
		shape[0] = shape[1] = 1;
		id[0] = 1; id[1] = 2;
		for (i = 0; i < 8; i++) {
			geometry[0][i] = 0;
			geometry[1][i] = 0;
		}
		CHECK(geomlines(0, 2, shape, id, geometry, NULL, NULL, GEPS) == 24);
	}
}
