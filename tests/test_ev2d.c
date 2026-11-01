/*
test_ev2d.c

後処理の 2D 作図 (post/ev2d.c の ev2 / HTML 出力、post/ev2dlib.c の作図部品) のテスト。

ofd_post の図 (周波数特性・遠方界・近傍界) はすべてここを通る。壊れても
ofd_post は普通に終わり、図の軸や値がずれることにしか現れない。
期待値は ev2 形式の定義 (post/ev2d.c 冒頭) と作図部品の約束から作る:

  -1 W H                ページ (幅・高さ)
  -2 R G B              色の変更 (直前と違うときだけ出る。ページの頭で黒に戻る)
  -3 X Y H / 文字列     文字列
  2 X1 Y1 X2 Y2         線分
  3 ... / 4 ...         三角形 / 四角形の塗りつぶし

  ev2dlib_grid  : 黒の外枠 (線分 4 本) と、灰 (192) の内側の線 (xdiv-1) + (ydiv-1) 本
  ev2dlib_func2 : 点 (x[i], y[i]) を枠 [x1,x2]×[y1,y2] へ一次変換で写し、
                  x[0]..x[div] が枠の左端..右端、ymin..ymax が下端..上端。
                  枠からはみ出す y は枠の端に切り詰める。div = 0 は水平線 1 本
                  (枠の外なら描かない)
  ev2dlib_func1 : x を 0..div の等間隔にした func2
  ev2d_setColorV : 白黒 (color = 0) は f = 0 で白、1 で黒、範囲外は切り詰め。
                  カラーは 20 段の表で f = 0 が青 (0,0,255)、1 が赤 (255,0,0)
  データが無いページは出力しない。1 ページ 1 万件を超えても全件出る
*/

#include "ev.h"
#include "ofd_test.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

extern unsigned char ev2d_rgb[3];
extern int ev2d_page;

#define TMPEV2 "ofd_test_ev2d.ev2"
#define TMPHTM "ofd_test_ev2d.htm"
#define MAXREC 40000

typedef struct { int idx; double v[8]; int rgb[3]; int page; } rec_t;

static rec_t *rec;
static int nrec, npage;
static double pw[8], ph[8];
static char str0[256];

/* ev2 を読み、各図形に「そのときの色」と「ページ番号」を付けて並べる */
static void read_ev2(void)
{
	nrec = 0; npage = 0;
	str0[0] = '\0';
	FILE *fp = fopen(TMPEV2, "r");
	if (fp == NULL) return;
	char line[1024];
	int rgb[3] = {0, 0, 0};
	while (fgets(line, sizeof(line), fp)) {
		int idx;
		if (sscanf(line, "%d", &idx) != 1) continue;
		if (idx == -1) {
			if (npage < 8) sscanf(line, "%d %lf %lf", &idx, &pw[npage], &ph[npage]);
			npage++;
			continue;
		}
		if (idx == -2) { sscanf(line, "%d %d %d %d", &idx, &rgb[0], &rgb[1], &rgb[2]); continue; }
		if (nrec >= MAXREC) continue;
		rec_t *r = &rec[nrec++];
		memset(r, 0, sizeof(*r));
		r->idx = idx; r->page = npage - 1;
		memcpy(r->rgb, rgb, sizeof(rgb));
		const int nv = (idx == 2) ? 4 : (idx == 3) ? 6 : (idx == 4) ? 8 : (idx == -3) ? 3 : 4;
		const char *p = line;
		while (*p == ' ' || *p == '-') p++;
		while (*p && *p != ' ') p++;
		for (int k = 0; k < nv; k++) r->v[k] = strtod(p, (char **)&p);
		if (idx == -3) {
			if (fgets(line, sizeof(line), fp) && (str0[0] == '\0')) {
				line[strcspn(line, "\r\n")] = '\0';
				strcpy(str0, line);
			}
		}
	}
	fclose(fp);
}

static int near(double a, double b) { return fabs(a - b) < 1e-6 * (fabs(a) + fabs(b) + 1); }

static int is_line(const rec_t *r, double x1, double y1, double x2, double y2)
{
	return (r->idx == 2) && near(r->v[0], x1) && near(r->v[1], y1) && near(r->v[2], x2) && near(r->v[3], y2);
}

static int rgb_is(const rec_t *r, int a, int b, int c) { return (r->rgb[0] == a) && (r->rgb[1] == b) && (r->rgb[2] == c); }

void test_ev2d(void)
{
	rec = (rec_t *)malloc(MAXREC * sizeof(rec_t));

	ev2d_init(1.5, 1.0);

	/* ---- ページ 1 : 図形と色 ---- */
	ev2d_newPage();
	ev2d_setColor(10, 20, 30);
	ev2d_drawLine(0.1, 0.2, 0.3, 0.4);
	ev2d_drawLine(0.5, 0.6, 0.7, 0.8);          /* 同じ色 : 色の行は 1 回だけ */
	ev2d_setColor(10, 20, 99);                  /* 青だけ違う色も「変更」 */
	ev2d_drawLine(0.9, 0.9, 1.0, 1.0);
	ev2d_setColor(200, 0, 0);
	ev2d_fillTriangle(0, 0, 1, 0, 0, 1);
	ev2d_fillRectangle(0.2, 0.3, 0.6, 0.9);     /* 四角形 (反時計回りの 4 頂点) */
	ev2d_drawRectangle(1, 2, 3, 4);             /* 線分 4 本 */
	ev2d_setColor(0, 0, 255);
	ev2d_drawString(0.05, 0.95, 0.03, "freq [GHz]");

	/* データの無いページは作らない (newPage を続けても空ページは出ない) */
	ev2d_newPage();
	ev2d_newPage();
	CHECK(ev2d_page == 1);

	/* ---- ページ 2 : 格子と関数 ---- */
	ev2dlib_grid(10, 20, 110, 220, 4, 2);
	ev2d_setColor(0, 128, 0);
	{
		/* x = 2, 3, 5 (不等間隔、先頭が 0 でない)、y = 0, 5, 20 を ymin = 0, ymax = 10 で
		   枠 (10,20)-(110,220) に */
		const double x[3] = {2, 3, 5}, y[3] = {0, 5, 20};
		ev2dlib_func2(2, x, y, 0, 10, 10, 20, 110, 220);
		/* div = 0 : y = 2.5 の水平線、y = 50 (枠の外) は描かない */
		const double y0[1] = {2.5}, y1[1] = {50};
		ev2dlib_func2(0, x, y0, 0, 10, 10, 20, 110, 220);
		ev2dlib_func2(0, x, y1, 0, 10, 10, 20, 110, 220);
		/* ymin == ymax は何もしない */
		ev2dlib_func2(2, x, y, 3, 3, 10, 20, 110, 220);
		/* func1 : x は 0..div の等間隔 (func2 が x[0]..x[div] を枠の幅に正規化するので、
		   等間隔の刻みの大きさ自体は結果に効かない) */
		const double f[3] = {10, 0, 10};
		ev2dlib_func1(2, f, 0, 10, 0, 0, 100, 100);
	}

	/* ---- ページ 3 : 1 万件を超える (確保の追加) ---- */
	ev2d_newPage();
	const int nbig = 25003;
	for (int n = 0; n < nbig; n++) ev2d_drawLine(n, 0, n, 1);

	ev2d_file(1, TMPEV2);
	ev2d_output();
	read_ev2();

	CHECK(npage == 3);
	CHECK(ev2d_page == 2);
	CHECK(near(pw[0], 1.5) && near(ph[0], 1.0));

	/* ページ 1 */
	{
		int n1 = 0, ok = 1;
		for (int m = 0; m < nrec; m++) if (rec[m].page == 0) n1++;
		CHECK(n1 == 3 + 1 + 1 + 4 + 1);
		ok &= is_line(&rec[0], 0.1, 0.2, 0.3, 0.4) && rgb_is(&rec[0], 10, 20, 30);
		ok &= is_line(&rec[1], 0.5, 0.6, 0.7, 0.8) && rgb_is(&rec[1], 10, 20, 30);
		ok &= is_line(&rec[2], 0.9, 0.9, 1.0, 1.0) && rgb_is(&rec[2], 10, 20, 99);
		ok &= (rec[3].idx == 3) && rgb_is(&rec[3], 200, 0, 0) && near(rec[3].v[2], 1) && near(rec[3].v[5], 1);
		ok &= (rec[4].idx == 4) && near(rec[4].v[0], 0.2) && near(rec[4].v[1], 0.3)
		   && near(rec[4].v[2], 0.6) && near(rec[4].v[3], 0.3)
		   && near(rec[4].v[4], 0.6) && near(rec[4].v[5], 0.9)
		   && near(rec[4].v[6], 0.2) && near(rec[4].v[7], 0.9);
		ok &= is_line(&rec[5], 1, 2, 3, 2) && is_line(&rec[6], 3, 2, 3, 4)
		   && is_line(&rec[7], 3, 4, 1, 4) && is_line(&rec[8], 1, 4, 1, 2);
		ok &= (rec[9].idx == -3) && rgb_is(&rec[9], 0, 0, 255) && near(rec[9].v[0], 0.05) && near(rec[9].v[2], 0.03);
		CHECK(ok);
		CHECK(strcmp(str0, "freq [GHz]") == 0);
	}
	/* 色の行は変わったときだけ : ファイル中の "-2" は
	     ページ 1 : 頭の黒 + (10,20,30) (10,20,99) (200,0,0) (0,0,255) = 5
	     ページ 2 : 頭の黒 + 灰 + 緑 (格子の後の黒には図形が無いので出ない) = 3
	     ページ 3 : 頭の黒 (newPage で黒に戻る)                       = 1 */
	{
		FILE *fp = fopen(TMPEV2, "r");
		int ncol = 0;
		char line[1024];
		while ((fp != NULL) && fgets(line, sizeof(line), fp)) if (!strncmp(line, "-2 ", 3)) ncol++;
		if (fp != NULL) fclose(fp);
		CHECK(ncol == 9);
	}

	/* ページ 2 : 格子 */
	{
		int nblack = 0, ngray = 0, ok = 1;
		int m = 0;
		while ((m < nrec) && (rec[m].page != 1)) m++;
		const int start = m;
		for (int k = 0; k < 4; k++, m++) nblack += (rec[m].idx == 2) && rgb_is(&rec[m], 0, 0, 0);
		for (int k = 0; k < 3 + 1; k++, m++) ngray += (rec[m].idx == 2) && rgb_is(&rec[m], 192, 192, 192);
		CHECK((nblack == 4) && (ngray == 4));
		/* 内側の縦線は x = 35, 60, 85、横線は y = 120 */
		ok &= is_line(&rec[start + 4], 35, 20, 35, 220) && is_line(&rec[start + 6], 85, 20, 85, 220);
		ok &= is_line(&rec[start + 7], 10, 120, 110, 120);
		CHECK(ok);

		/* func2 : 線分 2 本 (x は 10 + 100/3 * x、y は 20 + 20 * y、上端 220 で切る)、
		   水平線 1 本 (y = 70)。枠の外の水平線・ymin == ymax は描かない */
		const rec_t *f = &rec[start + 8];
		CHECK(is_line(&f[0], 10, 20, 10 + 100.0 / 3, 120) && rgb_is(&f[0], 0, 128, 0));
		CHECK(is_line(&f[1], 10 + 100.0 / 3, 120, 110, 220));
		CHECK(is_line(&f[2], 10, 70, 110, 70));
		/* func1 : (0,100)-(50,0)-(100,100) */
		CHECK(is_line(&f[3], 0, 100, 50, 0) && is_line(&f[4], 50, 0, 100, 100));
		int n2 = 0;
		for (int k = 0; k < nrec; k++) n2 += (rec[k].page == 1);
		CHECK(n2 == 8 + 3 + 2);
	}

	/* ページ 3 : 全件 */
	{
		int n3 = 0, last_ok = 0;
		for (int k = 0; k < nrec; k++) {
			if (rec[k].page == 2) {
				n3++;
				if (is_line(&rec[k], nbig - 1, 0, nbig - 1, 1)) last_ok = 1;
			}
		}
		CHECK(n3 == nbig);
		CHECK(last_ok);
	}

	/* HTML 版 : ページ数と JavaScript のデータ配列 */
	ev2d_file(0, TMPHTM);
	ev2d_output();
	{
		FILE *fp = fopen(TMPHTM, "r");
		CHECK(fp != NULL);
		int has_npage = 0, has_div = 0, nline_items = 0;
		char line[4096];
		while ((fp != NULL) && fgets(line, sizeof(line), fp)) {
			if (strstr(line, "var npage = 3;")) has_npage = 1;
			/* 初期表示のページ数も 3 (以前は内部のページ番号で "1/2" だった) */
			if (strstr(line, "<div id=\"page\">1/3</div>")) has_div = 1;
			if (!strncmp(line, "[2,", 3)) nline_items++;
		}
		if (fp != NULL) fclose(fp);
		CHECK(has_npage);
		CHECK(has_div);
		CHECK(nline_items == (3 + 4) + (8 + 3 + 2) + nbig);
	}
	remove(TMPEV2);
	remove(TMPHTM);

	/* 色の割り当て (ev2d_rgb を直接見る) */
	ev2d_setColorV(0.0, 0); CHECK((ev2d_rgb[0] == 255) && (ev2d_rgb[1] == 255) && (ev2d_rgb[2] == 255));
	ev2d_setColorV(1.0, 0); CHECK((ev2d_rgb[0] == 0) && (ev2d_rgb[2] == 0));
	ev2d_setColorV(-3.0, 0); CHECK(ev2d_rgb[0] == 255);
	ev2d_setColorV(0.5, 0); CHECK(ev2d_rgb[0] == 127);
	ev2d_setColorV(0.0, 1); CHECK((ev2d_rgb[0] == 0) && (ev2d_rgb[1] == 0) && (ev2d_rgb[2] == 255));
	ev2d_setColorV(1.0, 1); CHECK((ev2d_rgb[0] == 255) && (ev2d_rgb[1] == 0) && (ev2d_rgb[2] == 0));
	ev2d_setColorV(7.0, 1); CHECK((ev2d_rgb[0] == 255) && (ev2d_rgb[2] == 0));
	ev2d_setColorV(0.52, 1); CHECK((ev2d_rgb[0] == 0) && (ev2d_rgb[1] == 255) && (ev2d_rgb[2] == 0));   /* 10 段目 */

	free(rec);
	rec = NULL;
}
