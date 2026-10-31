/*
test_input.c

入力ファイルの解釈 input_data() (sol/input_data.c) のテスト。

期待値の出所 (入力形式の約束と .claude/rules/features.md の後方互換規則):
  - キーを省略したときは従来の既定値 (Mur、solver 3000/50/1e-3、HDF5 出力あり、
    内部抵抗 0、dt/Tw は自動) になること
  - メッシュ "x0 n1 x1 n2 x2 ..." は区間ごとの等分割で、節点はちょうど区切り点を通る
  - 周波数 "f0 f1 div" は div+1 点の等間隔。OpenTHFD の frequency は両方に入る
  - solver の反復数は nout の倍数に切り上げ、HDF5 の間隔も同様
  - 材料: type 1 (εr σ μr σm)、type 2 (分散性)、Ver. 2.1 以前の旧形式 (type 無し)
  - PML と PBC を同時に指定すると Mur に切り替わる (PBC と PML は併用できない)
  - CRLF の行末 (Windows で作ったファイル) を LF と同じに読むこと
  - 入力の誤りは 1 を返して止まり、落ちないこと。特に
      solver の nout = 0 (以前は 0 除算で SIGFPE)
      256 文字を超える title (以前は Title[256] があふれた)
      frequency2 の無い平面波 (以前は Freq2 の NULL 参照で SIGSEGV)
  - Ver. 2.1 以前の旧形式の material (6 トークン) を受け付けること
    (以前は一律に 7 トークンを要求していて旧形式の分岐に到達できなかった)
*/

#include "ofd.h"
#include "ofd_prototype.h"
#include "ofd_test.h"

#include <stdlib.h>
#include <string.h>

#define TMPNAME "ofd_test_input.ofd"

static void free_input(void)
{
	free(Material); Material = NULL; NMaterial = 0;
	free(Geometry); Geometry = NULL; NGeometry = 0;
	free(Feed); Feed = NULL; NFeed = 0;
	free(Point); Point = NULL; NPoint = 0;
	free(Inductor); Inductor = NULL; NInductor = 0;
	free(Tpa); Tpa = NULL; NTpa = 0;
	free(Freq1); Freq1 = NULL; NFreq1 = 0;
	free(Freq2); Freq2 = NULL; NFreq2 = 0;
	free(Xn); free(Yn); free(Zn); free(Xc); free(Yc); free(Zc);
	Xn = Yn = Zn = Xc = Yc = Zc = NULL;
	free(Gline); Gline = NULL;
	free(MGline); MGline = NULL;
	NGline = 0;
	IPlanewave = 0;
	WaveAmp = 0;
	PBCx = PBCy = PBCz = 0;
	iABC = 0;
}

/* text を一時ファイルに書いて input_data() に渡す */
static int parse(const char *text)
{
	FILE *fp = fopen(TMPNAME, "wb");
	if (fp == NULL) return -1;
	fputs(text, fp);
	fclose(fp);
	fp = fopen(TMPNAME, "r");
	if (fp == NULL) return -1;
	free_input();
	const int ret = input_data(fp);
	fclose(fp);
	remove(TMPNAME);
	return ret;
}

/* 最小の正しい入力 (メッシュ 2 軸は 1 区間、x は 2 区間) */
static const char *BASE =
	"OpenFDTD 4 2\n"
	"title = base\n"
	"xmesh = -0.05 10 -0.0125 6 0.0125\n"
	"ymesh = -0.04 8 0.04\n"
	"zmesh = 0 5 0.02\n"
	"feed = Z 0 0 0.01 1 0 50\n"
	"end\n";

/* BASE の "end" の前に行を足した入力 */
static int parse_with(const char *extra)
{
	static char buf[8192];
	const char *e = strstr(BASE, "end\n");
	const size_t head = (size_t)(e - BASE);
	memcpy(buf, BASE, head);
	buf[head] = '\0';
	strncat(buf, extra, sizeof(buf) - strlen(buf) - 8);
	strcat(buf, "end\n");
	return parse(buf);
}

static void test_defaults_and_mesh(void)
{
	CHECK(parse(BASE) == 0);

	/* 既定値 (キー省略時は従来と同じ) */
	CHECK(iABC == 0);
	CHECK(Solver.maxiter == 3000 && Solver.nout == 50);
	CHECK_NEAR(Solver.converg, 1e-3, 1e-15);
	CHECK(Hdf5Output == 1 && Hdf5Interval == 0);
	CHECK(rFeed == 0);
	CHECK(Dt == 0 && Tw == 0);
	CHECK(PBCx == 0 && PBCy == 0 && PBCz == 0);
	CHECK(NMaterial == 2);
	CHECK(NFreq1 == 0 && NFreq2 == 0);
	CHECK(strcmp(Title, "base") == 0);

	/* メッシュ : 区間ごとの等分割 */
	CHECK(Nx == 16 && Ny == 8 && Nz == 5);
	CHECK_NEAR(Xn[0], -0.05, 1e-15);
	CHECK_NEAR(Xn[10], -0.0125, 1e-15);
	CHECK_NEAR(Xn[16], 0.0125, 1e-15);
	CHECK_NEAR(Xn[5], -0.05 + 5 * 0.00375, 1e-15);      /* 区間 1 : 3.75 mm */
	CHECK_NEAR(Xn[13], -0.0125 + 3 * (0.025 / 6), 1e-15); /* 区間 2 : 4.17 mm */
	CHECK_NEAR(Yn[4], 0, 1e-15);
	CHECK_NEAR(Zc[2], 0.01, 1e-15);

	/* 給電 : Z 給電はセル中心の z に丸まる (z = 0.01 は Zc[2]) */
	CHECK(NFeed == 1);
	CHECK(Feed[0].dir == 'Z' && Feed[0].k == 2);
	CHECK_NEAR(Feed[0].volt, 1, 1e-15);
	CHECK_NEAR(Feed[0].z0, 50, 1e-15);

	/* CRLF の行末でも同じに読む */
	{
		char crlf[4096];
		size_t o = 0;
		for (const char *p = BASE; *p; p++) {
			if (*p == '\n') crlf[o++] = '\r';
			crlf[o++] = *p;
		}
		crlf[o] = '\0';
		CHECK(parse(crlf) == 0);
		CHECK(Nx == 16 && Ny == 8 && Nz == 5);
		CHECK(strcmp(Title, "base") == 0);
		CHECK(NFeed == 1 && Feed[0].k == 2);
	}

	/* コメント・空行は無視、end 以降は読まない */
	CHECK(parse_with("# comment = 1 2 3\n\nrfeed = 25\n") == 0);
	CHECK(rFeed == 25);
	{
		char buf[4096];
		snprintf(buf, sizeof(buf), "%srfeed = 99\n", BASE);   /* end の後 */
		CHECK(parse(buf) == 0);
		CHECK(rFeed == 0);
	}
}

static void test_keys(void)
{
	/* 周波数 : div+1 点の等間隔 */
	CHECK(parse_with("frequency1 = 1e9 3e9 4\nfrequency2 = 2e9 2e9 0\n") == 0);
	CHECK(NFreq1 == 5);
	CHECK_NEAR(Freq1[0], 1e9, 1);
	CHECK_NEAR(Freq1[2], 2e9, 1);
	CHECK_NEAR(Freq1[4], 3e9, 1);
	CHECK(NFreq2 == 1);
	CHECK_NEAR(Freq2[0], 2e9, 1);

	/* solver : maxiter は nout の倍数に切り上げ、hdf5 間隔も */
	CHECK(parse_with("solver = 1001 100 1e-4\nhdf5 = 1 150\n") == 0);
	CHECK(Solver.maxiter == 1100 && Solver.nout == 100);
	CHECK_NEAR(Solver.converg, 1e-4, 1e-18);
	CHECK(Hdf5Interval == 200);
	CHECK(parse_with("hdf5 = 0\n") == 0);
	CHECK(Hdf5Output == 0);

	/* 材料 : type 1 / type 2 */
	CHECK(parse_with("material = 1 2.5 0.1 1.5 20\nmaterial = 2 3 1e9 2e10 1e10\n") == 0);
	CHECK(NMaterial == 4);
	CHECK(Material[2].type == 1);
	CHECK_NEAR(Material[2].epsr, 2.5, 1e-15);
	CHECK_NEAR(Material[2].esgm, 0.1, 1e-15);
	CHECK_NEAR(Material[2].amur, 1.5, 1e-15);
	CHECK_NEAR(Material[2].msgm, 20, 1e-15);
	CHECK(Material[3].type == 2);
	CHECK_NEAR(Material[3].einf, 3, 1e-15);
	CHECK_NEAR(Material[3].ae, 1e9, 1);
	CHECK_NEAR(Material[3].be, 2e10, 1);
	CHECK_NEAR(Material[3].ce, 1e10, 1);

	/* 旧形式 (Ver. 2.1 : type の欄が無い) */
	{
		char buf[4096];
		snprintf(buf, sizeof(buf), "OpenFDTD 2 1\nxmesh = 0 4 1\nymesh = 0 4 1\nzmesh = 0 4 1\n"
			"material = 3.5 0.2 1.2 0\nfeed = Z 0.5 0.5 0.5 1 0 50\nend\n");
		CHECK(parse(buf) == 0);
		CHECK(NMaterial == 3 && Material[2].type == 1);
		CHECK_NEAR(Material[2].epsr, 3.5, 1e-15);
		CHECK_NEAR(Material[2].esgm, 0.2, 1e-15);
		CHECK_NEAR(Material[2].amur, 1.2, 1e-15);
	}

	/* 形状 : 直方体 (6 個) と円錐台 (8 個) のパラメータ */
	CHECK(parse_with("material = 1 2 0 1 0\n"
		"geometry = 2 1 -0.01 0.01 -0.02 0.02 0 0.01\n"
		"geometry = 1 53 0 0 0.005 0.002 0 0.01 0.003 0.001\n") == 0);
	CHECK(NGeometry == 2);
	CHECK(Geometry[0].m == 2 && Geometry[0].shape == 1);
	CHECK_NEAR(Geometry[0].g[3], 0.02, 1e-15);
	CHECK(Geometry[1].m == 1 && Geometry[1].shape == 53);
	CHECK_NEAR(Geometry[1].g[7], 0.001, 1e-15);

	/* 吸収境界 : PML、PBC と併用すると Mur に戻る */
	CHECK(parse_with("abc = 1 8 2 1e-6\n") == 0);
	CHECK(iABC == 1 && cPML.l == 8);
	CHECK_NEAR(cPML.m, 2, 1e-15);
	CHECK_NEAR(cPML.r0, 1e-6, 1e-20);
	CHECK(parse_with("abc = 1 8 2 1e-6\npbc = 1 1 0\n") == 0);
	CHECK(iABC == 0 && PBCx == 1 && PBCy == 1 && PBCz == 0);

	/* dt / Tw の指定 */
	CHECK(parse_with("timestep = 1.5e-12\npulsewidth = 2e-10\n") == 0);
	CHECK(Dt == 1.5e-12 && Tw == 2e-10);

	/* OpenTHFD の frequency は frequency1/2 の両方 */
	{
		char buf[4096];
		snprintf(buf, sizeof(buf), "OpenTHFD 1 0\nxmesh = 0 4 1\nymesh = 0 4 1\nzmesh = 0 4 1\n"
			"feed = Z 0.5 0.5 0.5 1 0 50\nfrequency = 1e12 2e12 2\nend\n");
		CHECK(parse(buf) == 0);
		CHECK(NFreq1 == 3 && NFreq2 == 3);
		CHECK_NEAR(Freq2[1], 1.5e12, 1);
	}

	/* 長い title は切り詰める (Title[256] をあふれさせない) */
	{
		char buf[8000], title[1000];
		memset(title, 'a', sizeof(title) - 1);
		title[sizeof(title) - 1] = '\0';
		snprintf(buf, sizeof(buf), "OpenFDTD 4 2\ntitle = %s\nxmesh = 0 4 1\nymesh = 0 4 1\nzmesh = 0 4 1\n"
			"feed = Z 0.5 0.5 0.5 1 0 50\nend\n", title);
		CHECK(parse(buf) == 0);
		CHECK(strlen(Title) == sizeof(Title) - 1);
		CHECK(Title[0] == 'a');
	}
}

static void test_errors(void)
{
	/* 見出し */
	CHECK(parse("OpenXYZ 1 0\nend\n") == 1);
	CHECK(parse("OpenFDTD\nend\n") == 1);

	/* メッシュ */
	CHECK(parse("OpenFDTD 4 2\nymesh = 0 4 1\nzmesh = 0 4 1\nfeed = Z 0 0 0 1 0 50\nend\n") == 1);
	CHECK(parse("OpenFDTD 4 2\nxmesh = 1 4 0\nymesh = 0 4 1\nzmesh = 0 4 1\nfeed = Z 0 0 0 1 0 50\nend\n") == 1);
	CHECK(parse("OpenFDTD 4 2\nxmesh = 0 0 1\nymesh = 0 4 1\nzmesh = 0 4 1\nfeed = Z 0 0 0 1 0 50\nend\n") == 1);
	CHECK(parse("OpenFDTD 4 2\nxmesh = 0 4\nymesh = 0 4 1\nzmesh = 0 4 1\nfeed = Z 0 0 0 1 0 50\nend\n") == 1);

	/* 波源 : 無し / 給電と平面波の両方 */
	CHECK(parse("OpenFDTD 4 2\nxmesh = 0 4 1\nymesh = 0 4 1\nzmesh = 0 4 1\nend\n") == 1);
	CHECK(parse_with("planewave = 90 0 1\n") == 1);

	/* 材料・形状 */
	CHECK(parse_with("material = 1 0 0 1 0\n") == 1);       /* εr <= 0 */
	CHECK(parse_with("material = 1 2 -1 1 0\n") == 1);      /* σ < 0 */
	CHECK(parse_with("material = 2 0 0 0 0\n") == 1);       /* ε∞ <= 0 */
	CHECK(parse_with("geometry = 5 1 0 1 0 1 0 1\n") == 1); /* 材料 id 範囲外 */
	CHECK(parse_with("geometry = 1 61 0 1 0 1 0 1\n") == 1); /* 未実装の形状 */
	CHECK(parse_with("geometry = 1 1 0 1 0\n") == 1);       /* パラメータ不足 */

	/* 給電・観測点・負荷の向き */
	CHECK(parse("OpenFDTD 4 2\nxmesh = 0 4 1\nymesh = 0 4 1\nzmesh = 0 4 1\nfeed = W 0 0 0 1 0 50\nend\n") == 1);
	CHECK(parse_with("point = Q 0 0 0.01 +X\n") == 1);
	CHECK(parse_with("load = X 0 0 0.01 Q 50\n") == 1);

	/* solver : nout = 0 は落ちずにエラー (以前は 0 除算) */
	CHECK(parse_with("solver = 100 0 1e-3\n") == 1);
	CHECK(parse_with("solver = 0 10 1e-3\n") == 1);

	/* 周波数の分割数が負 */
	CHECK(parse_with("frequency1 = 1e9 2e9 -1\n") == 1);

	/* 平面波は周波数が要る (無いとパルス幅を決められない)。
	   frequency1 だけでも動くこと (以前は Freq2 の NULL 参照で落ちた) */
	{
		char buf[4096];
		snprintf(buf, sizeof(buf), "OpenFDTD 4 2\nxmesh = 0 4 1\nymesh = 0 4 1\nzmesh = 0 4 1\n"
			"planewave = 90 0 1\nend\n");
		CHECK(parse(buf) == 1);
		snprintf(buf, sizeof(buf), "OpenFDTD 4 2\nxmesh = 0 4 1\nymesh = 0 4 1\nzmesh = 0 4 1\n"
			"planewave = 90 0 1\nfrequency1 = 1e8 3e8 2\nend\n");
		CHECK(parse(buf) == 0);
		/* パルス幅 a = 4 / (1.27 / f0)、f0 = frequency1 の中央 */
		CHECK_NEAR(Planewave.ai, 4 / (1.27 / 2e8), 1e-6 * Planewave.ai);
		snprintf(buf, sizeof(buf), "OpenFDTD 4 2\nxmesh = 0 4 1\nymesh = 0 4 1\nzmesh = 0 4 1\n"
			"planewave = 90 0 1\nfrequency1 = 1e8 3e8 2\nfrequency2 = 4e8 4e8 0\nend\n");
		CHECK(parse(buf) == 0);
		CHECK_NEAR(Planewave.ai, 4 / (1.27 / 4e8), 1e-6 * Planewave.ai);   /* frequency2 優先 */
	}

	/* waveamp は平面波専用、周波数が要る */
	CHECK(parse_with("waveamp = 1e6\n") == 1);
	{
		char buf[4096];
		snprintf(buf, sizeof(buf), "OpenFDTD 4 2\nxmesh = 0 4 1\nymesh = 0 4 1\nzmesh = 0 4 1\n"
			"planewave = 90 0 1\nwaveamp = 1e6\nend\n");
		CHECK(parse(buf) == 1);
		snprintf(buf, sizeof(buf), "OpenFDTD 4 2\nxmesh = 0 4 1\nymesh = 0 4 1\nzmesh = 0 4 1\n"
			"planewave = 90 0 1\nwaveamp = 1e6\nfrequency1 = 2e9 3e9 1\nend\n");
		CHECK(parse(buf) == 0);
		CHECK_NEAR(WaveOmega, 2 * PI * 2e9, 1e-3);
	}
}

void test_input(void)
{
	const double dt0 = Dt, tw0 = Tw;
	test_defaults_and_mesh();
	test_keys();
	test_errors();
	free_input();
	Dt = dt0;
	Tw = tw0;
}
