/*
test_ev3dlib.c

遠方界の 3D パターン ev3dlib_func() (post/ev3dlib.c) のテスト。

plotfar2d の 3D 図 (post/plot3dFar2d0.c) は r(θ,φ) を球座標の四角形の面で描く。
約束:
  - θ = th1..th2 を nth 等分、φ = ph1..ph2 を nph 等分し、各 (ith, iph) の
    四角形の 4 頂点は (ith,iph), (ith,iph+1), (ith+1,iph+1), (ith+1,iph) の順
  - 頂点は原点から距離 r + off、向き (sinθ cosφ, sinθ sinφ, cosθ)
  - r に負の値があれば、最小値が 0 になるよう全体を持ち上げる (off = -rmin)。
    dB 表示で負になる値でも形が保たれる
  - 面の色は 4 頂点の平均 / 最大 で決める (最大の面は赤 255,0,0)
  - nth <= 0 / nph <= 0 では何も描かない (nth のガードを外しても面のループは回らないので
    等価。0 割りの値を使う前に止める意味で残っている)
ev3 形式 (sol/ev3d.c 冒頭) : 4 I X1 Y1 Z1 ... X4 Y4 Z4 R G B
*/

#include "ev.h"
#include "ofd_test.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

void ev3dlib_func(int nth, double th1, double th2, int nph, double ph1, double ph2, double **rdata);
extern unsigned char ev3d_rgb[3];

#define TMPEV3 "ofd_test_ev3dlib.ev3"
#define DTOR   (0.017453292519943295)

typedef struct { double v[4][3]; int rgb[3]; } evquad_t;

static int read_quads(evquad_t *q, int maxq)
{
	FILE *fp = fopen(TMPEV3, "r");
	if (fp == NULL) return -1;
	char line[1024];
	int n = 0;
	while (fgets(line, sizeof(line), fp)) {
		int idx, num;
		double a[12];
		int rgb[3];
		if ((sscanf(line, "%d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %d %d %d", &idx, &num,
			&a[0], &a[1], &a[2], &a[3], &a[4], &a[5], &a[6], &a[7], &a[8], &a[9], &a[10], &a[11],
			&rgb[0], &rgb[1], &rgb[2]) == 17) && (idx == 4) && (n < maxq)) {
			for (int k = 0; k < 4; k++) for (int c = 0; c < 3; c++) q[n].v[k][c] = a[3 * k + c];
			memcpy(q[n].rgb, rgb, sizeof(rgb));
			n++;
		}
	}
	fclose(fp);
	return n;
}

static double rfun(double th, double ph) { return 1 + th / 180 + ph / 720; }

/* 頂点 v が (th, ph) の向きで距離 r にあるか */
static int at(const double v[3], double r, double th, double ph)
{
	const double e[3] = {r * sin(th * DTOR) * cos(ph * DTOR), r * sin(th * DTOR) * sin(ph * DTOR), r * cos(th * DTOR)};
	return (fabs(v[0] - e[0]) < 1e-5) && (fabs(v[1] - e[1]) < 1e-5) && (fabs(v[2] - e[2]) < 1e-5);
}

static int run_case(int nth, int nph, double shift, evquad_t *q, int maxq)
{
	double **r = (double **)malloc((nth + 1) * sizeof(double *));
	for (int i = 0; i <= nth; i++) {
		r[i] = (double *)malloc((nph + 1) * sizeof(double));
		for (int j = 0; j <= nph; j++) r[i][j] = rfun((nth > 0) ? 180.0 * i / nth : 0, 360.0 * j / nph) + shift;
	}
	ev3d_init();
	ev3d_newPage();
	ev3d_drawLine(0, 0, 0, 0, 0, 0);       /* ページを空にしない (nth = 0 の場合) */
	ev3dlib_func(nth, 0, 180, nph, 0, 360, r);
	ev3d_file(1, TMPEV3, 0);
	ev3d_output();
	for (int i = 0; i <= nth; i++) free(r[i]);
	free(r);
	const int n = read_quads(q, maxq);
	remove(TMPEV3);
	return n;
}

void test_ev3dlib(void)
{
	const int nth = 6, nph = 8;
	evquad_t *q = (evquad_t *)malloc(1000 * sizeof(evquad_t));

	/* (a) 正の r : 面の数、頂点の位置と順、色 */
	CHECK(run_case(nth, nph, 0, q, 1000) == nth * nph);
	{
		int ok = 1, m = 0, nred = 0;
		double fmax = 0;
		int imax = -1;
		for (int i = 0; i < nth; i++) {
		for (int j = 0; j < nph; j++, m++) {
			const double t1 = 180.0 * i / nth, t2 = 180.0 * (i + 1) / nth;
			const double p1 = 360.0 * j / nph, p2 = 360.0 * (j + 1) / nph;
			ok &= at(q[m].v[0], rfun(t1, p1), t1, p1);
			ok &= at(q[m].v[1], rfun(t1, p2), t1, p2);
			ok &= at(q[m].v[2], rfun(t2, p2), t2, p2);
			ok &= at(q[m].v[3], rfun(t2, p1), t2, p1);
			const double rav = (rfun(t1, p1) + rfun(t1, p2) + rfun(t2, p2) + rfun(t2, p1)) / 4;
			if (rav > fmax) { fmax = rav; imax = m; }
			nred += (q[m].rgb[0] == 255) && (q[m].rgb[1] == 0) && (q[m].rgb[2] == 0);
		}
		}
		CHECK(ok);
		/* 平均が最大の面 (θ, φ とも最後) は最大値 rfun(180,360) = 2.5 に近く赤。
		   最小の面 (θ, φ とも最初、平均 ≈ 1.16 / 2.5) は赤でない */
		CHECK((imax == nth * nph - 1) && (q[imax].rgb[0] == 255) && (q[imax].rgb[2] == 0));
		CHECK(!((q[0].rgb[0] == 255) && (q[0].rgb[1] == 0)));
		CHECK(nred < nth * nph);
	}

	/* (b) 負の値を含む r : 最小値が 0 になるよう持ち上げる (rfun - 1.5 の最小は -0.5) */
	CHECK(run_case(nth, nph, -1.5, q, 1000) == nth * nph);
	{
		int ok = 1, m = 0;
		for (int i = 0; i < nth; i++) {
		for (int j = 0; j < nph; j++, m++) {
			const double t1 = 180.0 * i / nth, t2 = 180.0 * (i + 1) / nth;
			const double p1 = 360.0 * j / nph, p2 = 360.0 * (j + 1) / nph;
			ok &= at(q[m].v[0], rfun(t1, p1) - 1.0, t1, p1);
			ok &= at(q[m].v[2], rfun(t2, p2) - 1.0, t2, p2);
		}
		}
		CHECK(ok);
		/* 色も持ち上げた後の値で決まる : 最初の面の平均 (持ち上げ後) / 最大 (持ち上げ後 1.5) */
		const double t = 180.0 / nth, p = 360.0 / nph;
		const double rav = (rfun(0, 0) + rfun(0, p) + rfun(t, p) + rfun(t, 0)) / 4 - 1.0;
		ev3d_setColorV(rav / 1.5, 1);
		CHECK((q[0].rgb[0] == ev3d_rgb[0]) && (q[0].rgb[1] == ev3d_rgb[1]) && (q[0].rgb[2] == ev3d_rgb[2]));
	}

	/* (c) 分割数 0 は何も描かない */
	CHECK(run_case(0, nph, 0, q, 1000) == 0);

	/* (d) HTML 版のページ表示 (sol/ev3d.c) : 最後の空ページは数えない。
	   以前は ev3d_page + 1 を出していて、空ページを足した "1/2" になっていた */
	{
		ev3d_init();
		ev3d_newPage();
		ev3d_drawLine(0, 0, 0, 1, 1, 1);
		ev3d_newPage();                         /* 空のまま出力する */
		ev3d_file(0, "ofd_test_ev3dlib.htm", 0);
		ev3d_output();
		FILE *fp = fopen("ofd_test_ev3dlib.htm", "r");
		int has = 0;
		char line[4096];
		while ((fp != NULL) && fgets(line, sizeof(line), fp)) if (strstr(line, "<div id=\"page\">1/1</div>")) has = 1;
		if (fp != NULL) fclose(fp);
		remove("ofd_test_ev3dlib.htm");
		CHECK(has);
	}

	free(q);
}
