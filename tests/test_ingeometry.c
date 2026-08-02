/*
test_ingeometry.c

sol/ingeometry.c の形状内外判定のユニットテスト。

ingeometry() は物性値番号を Yee 格子に焼き付ける setupId の中核で、
ここが狂うと解析対象の形状そのものが変わる。形状コードごとに
「確実に内側の点」「確実に外側の点」を置いて判定を固定する。

形状コード : 1=直方体, 2=楕円体, 11/12/13=X/Y/Z 円柱,
             31/32/33=X/Y/Z 三角柱, 41/42/43=X/Y/Z 角錐台,
             51/52/53=X/Y/Z 円錐台
*/

#include "ofd.h"
#include "ofd_test.h"

int ingeometry(double x, double y, double z, int shape, double *g, double eps);

#define GEPS (1e-6)

static void test_box(void)
{
	/* [0,1] x [0,2] x [0,3] の直方体 */
	double g[8] = {0, 1, 0, 2, 0, 3, 0, 0};

	CHECK(ingeometry(0.5, 1.0, 1.5, 1, g, GEPS) == 1);   /* 中心 */
	CHECK(ingeometry(0.0, 0.0, 0.0, 1, g, GEPS) == 1);   /* 角 (境界は内側) */
	CHECK(ingeometry(1.0, 2.0, 3.0, 1, g, GEPS) == 1);   /* 反対の角 */

	CHECK(ingeometry(-0.5, 1.0, 1.5, 1, g, GEPS) == 0);  /* X 方向に外 */
	CHECK(ingeometry(0.5, -0.5, 1.5, 1, g, GEPS) == 0);  /* Y 方向に外 */
	CHECK(ingeometry(0.5, 1.0, 3.5, 1, g, GEPS) == 0);   /* Z 方向に外 */

	/* 厚さ 0 の板 (平面) も判定できること */
	{
		double p[8] = {0.5, 0.5, 0, 2, 0, 3, 0, 0};
		CHECK(ingeometry(0.5, 1.0, 1.5, 1, p, GEPS) == 1);
		CHECK(ingeometry(0.6, 1.0, 1.5, 1, p, GEPS) == 0);
	}
}

static void test_ellipsoid(void)
{
	/* 原点中心・半径 1 の球 (bounding box は [-1,1]^3) */
	double g[8] = {-1, 1, -1, 1, -1, 1, 0, 0};

	CHECK(ingeometry(0.0, 0.0, 0.0, 2, g, GEPS) == 1);    /* 中心 */
	CHECK(ingeometry(0.9, 0.0, 0.0, 2, g, GEPS) == 1);    /* 内側 */
	CHECK(ingeometry(1.0, 0.0, 0.0, 2, g, GEPS) == 1);    /* 表面 */
	CHECK(ingeometry(1.1, 0.0, 0.0, 2, g, GEPS) == 0);    /* 外側 */

	/* 角は bounding box の中だが球の外 */
	CHECK(ingeometry(1.0, 1.0, 1.0, 2, g, GEPS) == 0);

	/* 対角方向 : |r| = 1 の点は内側、それを少し超えると外側 */
	{
		const double a = 1.0 / sqrt(3.0);
		CHECK(ingeometry(a, a, a, 2, g, GEPS) == 1);
		CHECK(ingeometry(1.01 * a, 1.01 * a, 1.01 * a, 2, g, GEPS) == 0);
	}

	/* 扁平な楕円体 : 軸ごとに半径が異なること */
	{
		double e[8] = {-2, 2, -1, 1, -1, 1, 0, 0};
		CHECK(ingeometry(1.9, 0.0, 0.0, 2, e, GEPS) == 1);
		CHECK(ingeometry(0.0, 1.9, 0.0, 2, e, GEPS) == 0);
	}
}

static void test_cylinder(void)
{
	/* Z 円柱 : 底面は原点中心・半径 1、Z は [0,2] */
	double g[8] = {-1, 1, -1, 1, 0, 2, 0, 0};

	CHECK(ingeometry(0.0, 0.0, 1.0, 13, g, GEPS) == 1);   /* 軸上 */
	CHECK(ingeometry(0.9, 0.0, 1.0, 13, g, GEPS) == 1);   /* 側面の内側 */
	CHECK(ingeometry(1.1, 0.0, 1.0, 13, g, GEPS) == 0);   /* 半径方向に外 */
	CHECK(ingeometry(0.0, 0.0, 2.5, 13, g, GEPS) == 0);   /* 軸方向に外 */
	CHECK(ingeometry(0.0, 0.0, 0.0, 13, g, GEPS) == 1);   /* 底面 */
	CHECK(ingeometry(0.0, 0.0, 2.0, 13, g, GEPS) == 1);   /* 上面 */
	/* 角 (半径方向に外) */
	CHECK(ingeometry(1.0, 1.0, 1.0, 13, g, GEPS) == 0);

	/* X 円柱 : X が [0,2]、YZ 断面が半径 1 */
	{
		double x[8] = {0, 2, -1, 1, -1, 1, 0, 0};
		CHECK(ingeometry(1.0, 0.0, 0.0, 11, x, GEPS) == 1);
		CHECK(ingeometry(1.0, 1.1, 0.0, 11, x, GEPS) == 0);
		CHECK(ingeometry(2.5, 0.0, 0.0, 11, x, GEPS) == 0);
	}

	/* Y 円柱 : Y が [0,2]、ZX 断面が半径 1 */
	{
		double y[8] = {-1, 1, 0, 2, -1, 1, 0, 0};
		CHECK(ingeometry(0.0, 1.0, 0.0, 12, y, GEPS) == 1);
		CHECK(ingeometry(0.0, 1.0, 1.1, 12, y, GEPS) == 0);
		CHECK(ingeometry(0.0, 2.5, 0.0, 12, y, GEPS) == 0);
	}
}

static void test_pillar(void)
{
	/* Z 三角柱 : Z が [0,1]、XY 断面は (0,0), (1,0), (0,1) の三角形
	   g = {z1, z2, x1, x2, x3, y1, y2, y3} */
	double g[8] = {0, 1, 0, 1, 0, 0, 0, 1};

	CHECK(ingeometry(0.2, 0.2, 0.5, 33, g, GEPS) == 1);   /* 三角形の内側 */
	CHECK(ingeometry(0.8, 0.8, 0.5, 33, g, GEPS) == 0);   /* 斜辺の外側 */
	CHECK(ingeometry(-0.2, 0.2, 0.5, 33, g, GEPS) == 0);  /* 三角形の外側 */
	CHECK(ingeometry(0.2, 0.2, 1.5, 33, g, GEPS) == 0);   /* 軸方向に外 */
	CHECK(ingeometry(0.0, 0.0, 0.5, 33, g, GEPS) == 1);   /* 頂点 */

	/* 退化した三角形 (面積 0) は常に外側を返す */
	{
		double d[8] = {0, 1, 0, 1, 2, 0, 0, 0};
		CHECK(ingeometry(0.5, 0.0, 0.5, 33, d, GEPS) == 0);
	}
}

static void test_pyramid_cone(void)
{
	/* Z 角錐台 : Z が [0,1]、中心 (0,0)、底面 2x2、上面 1x1
	   g = {z1, z2, x0, y0, h1x, h1y, h2x, h2y} */
	double g[8] = {0, 1, 0, 0, 2, 2, 1, 1};

	CHECK(ingeometry(0.0, 0.0, 0.5, 43, g, GEPS) == 1);   /* 軸上 */
	CHECK(ingeometry(0.9, 0.0, 0.05, 43, g, GEPS) == 1);  /* 底面近くは広い */
	CHECK(ingeometry(0.9, 0.0, 0.95, 43, g, GEPS) == 0);  /* 上面近くは狭い */
	CHECK(ingeometry(0.0, 0.0, 1.5, 43, g, GEPS) == 0);   /* 軸方向に外 */

	/* Z 円錐台 : 同じ寸法だが断面が楕円
	   角の点 (x=y) は角錐台では内側、円錐台では外側になる */
	{
		double c[8] = {0, 1, 0, 0, 2, 2, 2, 2};   /* 実質は円柱 (半径 1) */
		CHECK(ingeometry(0.0, 0.0, 0.5, 53, c, GEPS) == 1);
		CHECK(ingeometry(0.9, 0.0, 0.5, 53, c, GEPS) == 1);
		CHECK(ingeometry(1.1, 0.0, 0.5, 53, c, GEPS) == 0);
		CHECK(ingeometry(0.9, 0.9, 0.5, 53, c, GEPS) == 0);   /* 円の外 */
	}
}

static void test_unknown_shape(void)
{
	/* 未知の形状コードは常に 0 (外側) を返すこと */
	double g[8] = {0, 1, 0, 1, 0, 1, 0, 0};
	CHECK(ingeometry(0.5, 0.5, 0.5, 999, g, GEPS) == 0);
	CHECK(ingeometry(0.5, 0.5, 0.5, 0, g, GEPS) == 0);
}

void test_ingeometry(void)
{
	test_box();
	test_ellipsoid();
	test_cylinder();
	test_pillar();
	test_pyramid_cone();
	test_unknown_shape();
}
