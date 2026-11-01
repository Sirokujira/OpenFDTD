/*
test_postdata.c

後処理の入力の解釈 post_data() (post/post_data.c) のテスト。

.ofd はソルバーと ofd_post が同じファイルを読む。ofd_post の側の約束:
  - ソルバーのキー (xmesh, geometry, feed, solver, ...) は黙って無視する
    (.claude/rules/features.md : 新しいキーで ofd_post が落ちてはいけない)
  - '#' の行と空行は読まない、"end" 以降は読まない
  - 改行は LF でも CRLF でもよい (data/sample の多くは CRLF)
  - plotfar1d は X/Y/Z (分割数) と V/H (分割数 + 固定角)、plotnear1d は
    成分 方向 位置 位置、plotnear2d は 成分 方向 位置。plotzin などは
    "1" で自動スケール、"2 min max div" で指定スケール
  - 誤り (方向・分割数・成分名の長さ・先頭行) は 1 を返す
期待値は README / 入力形式の約束から作る。
*/

#include "ofd.h"
#include "ofd_prototype.h"
#include "ofd_test.h"

#include <stdlib.h>
#include <string.h>

int post_data(FILE *fp);

#define TMPIN "ofd_test_postdata.ofd"

static int run(const char *text)
{
	FILE *fp = fopen(TMPIN, "wb");
	if (fp == NULL) return -1;
	fputs(text, fp);
	fclose(fp);
	fp = fopen(TMPIN, "rb");
	if (fp == NULL) return -1;
	const int ret = post_data(fp);
	fclose(fp);
	remove(TMPIN);
	return ret;
}

static void reset(void)
{
	free(Far1d);  Far1d = NULL;
	free(Near1d); Near1d = NULL;
	free(Near2d); Near2d = NULL;
	NFar1d = NNear1d = NNear2d = 0;
}

void test_postdata(void)
{
	const int html0 = HTML;
	HTML = 0;

	/* ---- 正常な入力 : ソルバーのキーを混ぜ、CRLF・コメント・空行を含む ---- */
	const char *ok =
		"OpenFDTD 4 2\r\n"
		"title = post test\r\n"
		"xmesh = -0.1 20 0.1\r\n"
		"geometry = 1 1 0 0 0 0 -0.025 0.025\r\n"
		"feed = Z 0 0 0 1 0 50\r\n"
		"future_solver_key = 1 2 3\r\n"          /* 知らないキー : 無視 */
		"# plotfar1d = Z 99\r\n"                  /* コメント : 読まない */
		"\r\n"
		"\n"
		"matchingloss = 1\r\n"
		"plotiter = 1\n"                          /* LF の行も混ぜる */
		"plotfeed = 1\r\n"
		"plotsmith = 1\r\n"
		"plotzin = 1\r\n"
		"plotref = 2 -30 0 6\r\n"
		"freqdiv = 4\r\n"
		"plotfar0d = 45 30 2 -20 10 3\r\n"
		"plotfar1d = x 72\r\n"                    /* 小文字も可 */
		"plotfar1d = V 36 60\r\n"
		"plotfar1d = H 24 30\r\n"
		"far1dcomponent = 1 1 0\r\n"
		"far1dnorm = 1\r\n"
		"far1dscale = -30 10 4\r\n"
		"plotfar2d = 18 36\r\n"
		"far2dcomponent = 1 0 1 0 0 0 1\r\n"
		"far2dobj = 0.25\r\n"
		"plotnear1d = Ex Y 0.01 -0.02\r\n"
		"plotnear1d = H z 0.03 0.04\r\n"
		"plotnear2d = Ez X 0.005\r\n"
		"near2ddim = 0 1\r\n"
		"near2dnoinc = 1\r\n"
		"near2dzoom = 0.1 -0.1 0.2 -0.2\r\n"
		"end\r\n"
		"plotfar1d = Z 99\r\n";                   /* end の後 : 読まない */
	CHECK(run(ok) == 0);
	CHECK(MatchingLoss == 1);
	CHECK((Piter == 1) && (Pfeed == 1) && (Ppoint == 0));
	CHECK((IFreq[0] == 1) && (IFreq[1] == 1) && (IFreq[2] == 0) && (IFreq[3] == 1));
	CHECK((FreqScale[1].user == 0) && (FreqScale[3].user == 1));
	CHECK((FreqScale[3].min == -30) && (FreqScale[3].max == 0) && (FreqScale[3].div == 6));
	CHECK(Freqdiv == 4);
	CHECK((IFar0d == 1) && (Far0d[0] == 45) && (Far0d[1] == 30) && (Far0dScale.user == 1) && (Far0dScale.div == 3));
	CHECK(NFar1d == 3);
	if (NFar1d == 3) {
		CHECK((Far1d[0].dir == 'X') && (Far1d[0].div == 72));
		CHECK((Far1d[1].dir == 'V') && (Far1d[1].div == 36) && (Far1d[1].angle == 60));
		CHECK((Far1d[2].dir == 'H') && (Far1d[2].div == 24) && (Far1d[2].angle == 30));
	}
	CHECK((Far1dComp[0] == 1) && (Far1dComp[1] == 1) && (Far1dComp[2] == 0) && (Far1dNorm == 1));
	CHECK((Far1dScale.user == 1) && (Far1dScale.min == -30) && (Far1dScale.max == 10) && (Far1dScale.div == 4));
	CHECK((NFar2d == 1) && (Far2d.divtheta == 18) && (Far2d.divphi == 36));
	CHECK((Far2dComp[0] == 1) && (Far2dComp[2] == 1) && (Far2dComp[6] == 1) && (Far2dComp[1] == 0));
	CHECK(Far2dObj == 0.25);
	CHECK(NNear1d == 2);
	if (NNear1d == 2) {
		CHECK(!strcmp(Near1d[0].cmp, "Ex") && (Near1d[0].dir == 'Y') && (Near1d[0].pos1 == 0.01) && (Near1d[0].pos2 == -0.02));
		CHECK(!strcmp(Near1d[1].cmp, "H") && (Near1d[1].dir == 'Z'));
	}
	CHECK(NNear2d == 1);
	if (NNear2d == 1) {
		CHECK(!strcmp(Near2d[0].cmp, "Ez") && (Near2d[0].dir == 'X') && (Near2d[0].pos0 == 0.005));
	}
	CHECK((Near2dDim[0] == 0) && (Near2dDim[1] == 1) && (Near2dNoinc == 1));
	/* zoom は min/max に並べ替える */
	CHECK((Near2dIzoom == 1) && (Near2dHzoom[0] == -0.1) && (Near2dHzoom[1] == 0.1));
	CHECK((Near2dVzoom[0] == -0.2) && (Near2dVzoom[1] == 0.2));
	reset();

	/* ---- 既定値 : 前の呼び出しの値が残らない ---- */
	CHECK(run("OpenFDTD 4 2\nend\n") == 0);
	CHECK((MatchingLoss == 0) && (Piter == 0) && (Pfeed == 0) && (Freqdiv == 10));
	CHECK((IFreq[0] == 0) && (IFreq[1] == 0) && (IFreq[3] == 0) && (IFar0d == 0));
	CHECK((NFar1d == 0) && (NFar2d == 0) && (NNear1d == 0) && (NNear2d == 0));
	CHECK((Far1dNorm == 0) && (Far1dScale.user == 0) && (Far1dScale.db == 1));
	CHECK((Far2dComp[0] == 1) && (Far2dComp[2] == 0) && (Far2dObj == 0.5));
	CHECK((Near2dDim[0] == 1) && (Near2dDim[1] == 1) && (Near2dNoinc == 0) && (Near2dIzoom == 0));
	reset();

	/* ---- 誤り ---- */
	const char *bad[] = {
		"NotOpenFDTD 4 2\nend\n",                         /* 先頭行 */
		"OpenFDTD 4 2\nplotfar1d = W 36\nend\n",          /* 方向 */
		"OpenFDTD 4 2\nplotfar1d = V 36\nend\n",          /* V に固定角が無い */
		"OpenFDTD 4 2\nplotfar1d = X 1\nend\n",           /* 分割数 < 2 */
		"OpenFDTD 4 2\nplotfar2d = 0 36\nend\n",          /* far2d の分割数 */
		"OpenFDTD 4 2\nplotnear1d = Exy X 0 0\nend\n",    /* 成分名が長い */
		"OpenFDTD 4 2\nplotnear1d = E W 0 0\nend\n",      /* 方向 */
		"OpenFDTD 4 2\nplotnear1d = E X 0\nend\n",        /* 位置が足りない */
		"OpenFDTD 4 2\nplotnear2d = E XY 0\nend\n",       /* 方向が 1 文字でない */
		"OpenFDTD 4 2\nplotnear2d = E Q 0\nend\n",        /* 方向 */
	};
	const int nbad = (int)(sizeof(bad) / sizeof(bad[0]));
	int caught = 0;
	for (int n = 0; n < nbad; n++) {
		const int r = run(bad[n]);
		if (r == 1) caught++;
		else printf("    post_data accepted invalid input #%d (ret %d)\n", n, r);
		reset();
	}
	CHECK(caught == nbad);

	/* 先頭行より前のコメント : 読み飛ばして先頭行を探す (キーの行の "# key = ..." は
	   コメントでなくても "#" というキーとして無視されるので、コメントの扱いが
	   効くのはここだけ) */
	CHECK(run("# generated by OpenFDTD-X\nOpenFDTD 4 2\nplotiter = 1\nend\n") == 0);
	CHECK(Piter == 1);
	reset();

	/* OpenTHFD の先頭行も受け付ける */
	CHECK(run("OpenTHFD 1 0\nend\n") == 0);
	reset();

	HTML = html0;
}
