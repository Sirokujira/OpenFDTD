/*
test_plot3d.c

形状の 3D 表示 plot3dGeom() (sol/plot3dGeom.c, sol/ev3d.c) のテスト。

plot3dgeom = 1 で出る geom.ev3 (ASCII) と geom3d.htm を作らせ、ev3 を読んで
描かれた線分を数える。ev3 の形式 (sol/ev3d.c 冒頭):
    -1                                   ページ
    2 I X1 Y1 Z1 X2 Y2 Z2 R G B          線分
    -4 / 文字列                          タイトル

期待値 (表示の約束):
  - 形状の線分 (Gline) はすべて描かれ、材料 id で色が分かれる
    (PEC = 黒 0,0,0 / それ以外 = マゼンタ 255,0,255)。MGline の参照を検査する
  - 格子線 (灰 200,200,200) は各軸の節点ごとに 2 本 : 2 (Nx+1) + 2 (Ny+1) + 2 (Nz+1)
  - 給電 (赤)・インダクタ (橙 255,165,0)・観測点 (緑、差分用の 2 点を含む) は
    その辺 (向きの軸に 1 セル) の両端を結ぶ線分
  - タイトルに Title、形状数、Nx/Ny/Nz
  - HTML 版も書かれること
*/

#include "ofd.h"
#include "ofd_prototype.h"
#include "ofd_test.h"

#include <stdlib.h>
#include <string.h>

typedef struct { double a[3], b[3]; int rgb[3]; } seg_t;

static int same(double x, double y) { return fabs(x - y) <= 1e-6 * (fabs(x) + fabs(y) + 1e-12); }

/* 線分 s が (p, q) を結ぶか (向きは問わない) */
static int connects(const seg_t *s, const double p[3], const double q[3])
{
	int f = 1, r = 1;
	for (int c = 0; c < 3; c++) {
		f &= same(s->a[c], p[c]) && same(s->b[c], q[c]);
		r &= same(s->a[c], q[c]) && same(s->b[c], p[c]);
	}
	return f || r;
}

static int color(const seg_t *s, int r, int g, int b) { return (s->rgb[0] == r) && (s->rgb[1] == g) && (s->rgb[2] == b); }

/* 向き dir の辺 (i,j,k) の両端 */
static void edge(char dir, int i, int j, int k, double p[3], double q[3])
{
	p[0] = Xn[i]; p[1] = Yn[j]; p[2] = Zn[k];
	q[0] = Xn[i + (dir == 'X')]; q[1] = Yn[j + (dir == 'Y')]; q[2] = Zn[k + (dir == 'Z')];
}

void test_plot3d(void)
{
	char title0[256];
	memcpy(title0, Title, sizeof(title0));
	strcpy(Title, "plot3d test");

	Nx = 4; Ny = 3; Nz = 5;
	Xn = (double *)malloc((Nx + 1) * sizeof(double));
	Yn = (double *)malloc((Ny + 1) * sizeof(double));
	Zn = (double *)malloc((Nz + 1) * sizeof(double));
	for (int i = 0; i <= Nx; i++) Xn[i] = 0.01 * i;
	for (int j = 0; j <= Ny; j++) Yn[j] = 0.02 * j - 0.03;
	for (int k = 0; k <= Nz; k++) Zn[k] = 0.015 * k + 0.1;

	/* 形状の線分 : PEC 2 本、誘電体 3 本 */
	NGeometry = 2;
	NGline = 5;
	Gline = (double (*)[2][3])malloc(NGline * sizeof(*Gline));
	MGline = (id_t *)malloc(NGline * sizeof(id_t));
	const id_t mg[5] = {PEC, 2, 3, PEC, 2};
	for (int n = 0; n < NGline; n++) {
		for (int c = 0; c < 3; c++) {
			Gline[n][0][c] = 0.001 * (n + 1) + 0.1 * c;
			Gline[n][1][c] = 0.002 * (n + 1) + 0.1 * c + 0.05;
		}
		MGline[n] = mg[n];
	}

	/* 給電 2 (X, Z)、インダクタ 1 (Y)、観測点 1 + 差分用 2 */
	NFeed = 2;
	Feed = (feed_t *)calloc(NFeed, sizeof(feed_t));
	Feed[0].dir = 'X'; Feed[0].i = 1; Feed[0].j = 2; Feed[0].k = 3;
	Feed[1].dir = 'Z'; Feed[1].i = 3; Feed[1].j = 1; Feed[1].k = 0;
	NInductor = 1;
	Inductor = (inductor_t *)calloc(NInductor, sizeof(inductor_t));
	Inductor[0].dir = 'Y'; Inductor[0].i = 2; Inductor[0].j = 0; Inductor[0].k = 4;
	NPoint = 1;
	Point = (point_t *)calloc(NPoint + 2, sizeof(point_t));
	Point[0].dir = 'Z'; Point[0].i = 2; Point[0].j = 2; Point[0].k = 2;
	Point[1] = Point[0]; Point[1].i = 3;
	Point[2] = Point[0]; Point[2].i = 1;

	remove("geom.ev3");
	remove("geom3d.htm");
	plot3dGeom();

	/* ---- geom.ev3 を読む ---- */
	FILE *fp = fopen("geom.ev3", "r");
	CHECK(fp != NULL);
	seg_t *seg = (seg_t *)malloc(1000 * sizeof(seg_t));
	int nseg = 0, has_title = 0, has_ngeom = 0, has_size = 0, npage = 0;
	char line[1024];
	while ((fp != NULL) && fgets(line, sizeof(line), fp)) {
		int idx, num;
		seg_t s;
		if (!strncmp(line, "-1", 2)) { npage++; continue; }
		if (strstr(line, "plot3d test")) has_title = 1;
		if (strstr(line, "No. of geometries = 2")) has_ngeom = 1;
		if (strstr(line, "Nx=4 Ny=3 Nz=5")) has_size = 1;
		if ((sscanf(line, "%d %d %lf %lf %lf %lf %lf %lf %d %d %d", &idx, &num,
			&s.a[0], &s.a[1], &s.a[2], &s.b[0], &s.b[1], &s.b[2], &s.rgb[0], &s.rgb[1], &s.rgb[2]) == 11)
			&& (idx == 2) && (nseg < 1000)) {
			seg[nseg++] = s;
		}
	}
	if (fp != NULL) fclose(fp);

	CHECK(npage == 1);
	CHECK(has_title && has_ngeom && has_size);

	/* 形状の線分 : 全部、色は材料で */
	int geom_ok = 1;
	for (int n = 0; n < NGline; n++) {
		int found = 0;
		for (int m = 0; m < nseg; m++) {
			if (connects(&seg[m], Gline[n][0], Gline[n][1])) {
				found = (MGline[n] == PEC) ? color(&seg[m], 0, 0, 0) : color(&seg[m], 255, 0, 255);
			}
		}
		geom_ok &= found;
	}
	CHECK(geom_ok);

	/* 格子線 */
	int ngray = 0;
	for (int m = 0; m < nseg; m++) ngray += color(&seg[m], 200, 200, 200);
	CHECK(ngray == 2 * (Nx + 1) + 2 * (Ny + 1) + 2 * (Nz + 1));

	/* 給電・インダクタ・観測点はその辺 */
	double p[3], q[3];
	int ok_f = 1, ok_l = 1, ok_p = 1, nred = 0, norange = 0, ngreen = 0;
	for (int m = 0; m < nseg; m++) {
		nred += color(&seg[m], 255, 0, 0);
		norange += color(&seg[m], 255, 165, 0);
		ngreen += color(&seg[m], 0, 255, 0);
	}
	for (int n = 0; n < NFeed; n++) {
		edge(Feed[n].dir, Feed[n].i, Feed[n].j, Feed[n].k, p, q);
		int found = 0;
		for (int m = 0; m < nseg; m++) found |= connects(&seg[m], p, q) && color(&seg[m], 255, 0, 0);
		ok_f &= found;
	}
	edge(Inductor[0].dir, Inductor[0].i, Inductor[0].j, Inductor[0].k, p, q);
	{
		int found = 0;
		for (int m = 0; m < nseg; m++) found |= connects(&seg[m], p, q) && color(&seg[m], 255, 165, 0);
		ok_l &= found;
	}
	for (int n = 0; n < NPoint + 2; n++) {
		edge(Point[n].dir, Point[n].i, Point[n].j, Point[n].k, p, q);
		int found = 0;
		for (int m = 0; m < nseg; m++) found |= connects(&seg[m], p, q) && color(&seg[m], 0, 255, 0);
		ok_p &= found;
	}
	CHECK(ok_f && nred == NFeed);
	CHECK(ok_l && norange == NInductor);
	CHECK(ok_p && ngreen == NPoint + 2);
	CHECK(nseg == NGline + ngray + nred + norange + ngreen);

	/* HTML 版 */
	fp = fopen("geom3d.htm", "r");
	CHECK(fp != NULL);
	if (fp != NULL) {
		size_t n = fread(line, 1, sizeof(line) - 1, fp);
		line[n] = '\0';
		CHECK(strstr(line, "<") != NULL);
		/* ページ表示はデータのあるページ数 (1 ページ) */
		int has_div = 0;
		rewind(fp);
		while (fgets(line, sizeof(line), fp)) if (strstr(line, "<div id=\"page\">1/1</div>")) has_div = 1;
		CHECK(has_div);
		fclose(fp);
	}
	remove("geom.ev3");
	remove("geom3d.htm");

	free(seg);
	free(Xn); free(Yn); free(Zn);
	free(Gline); free(MGline);
	free(Feed); free(Inductor); free(Point);
	Xn = Yn = Zn = NULL;
	Gline = NULL; MGline = NULL;
	Feed = NULL; Inductor = NULL; Point = NULL;
	NGline = 0; NGeometry = 0; NFeed = 0; NInductor = 0; NPoint = 0;
	Nx = Ny = Nz = 0;
	memcpy(Title, title0, sizeof(title0));
}
