/*
test_setupid.c

sol/setupId.c (形状を Yee 格子へ焼き付ける処理) のユニットテスト。

setupId() は各 geometry について、6 成分それぞれの実位置

    Ex : (Xc[i], Yn[j], Zn[k])        Hx : (Xn[i], Yc[j], Zc[k])
    Ey : (Xn[i], Yc[j], Zn[k])        Hy : (Xc[i], Yn[j], Zc[k])
    Ez : (Xn[i], Yn[j], Zc[k])        Hz : (Xc[i], Yc[j], Zn[k])

を ingeometry() に渡して材料 ID (iEx..iHz) を焼き付ける。ここで節点
(Xn) とセル中心 (Xc) を取り違えると**材料が半セルずれる**が、実行は
正常に完走してしまい、結果の歪みでしか気付けない。このテストは
格子に対して半端な位置に直方体を置き、6 成分それぞれで「入るべき
格子点」の集合を固定する。

直方体は**軸ごとに寸法を変えてある** (x:[0.22,0.78], y:[0.42,0.68],
z:[0.12,0.48])。立方体だと座標の軸を取り違える変異 (y と z の
入れ替え等) が結果に現れないうえ、shape = 1 の高速経路が使う
getspan() は「境界をまたぐセルまで含めて広めに」範囲を返し、本体が
ingeometry() で再判定して削るため、**範囲の取り違えは結果の範囲が
狭くなる場合しか観測できない**。軸ごとに寸法を変え (x と y の
取り違えで範囲が狭くなる)、さらにドメイン端に接する直方体 (getspan の
探索上限の切り詰めが焼き付け漏れになる) を別に置いて、その方向の
間違いを捕まえる。境界 0.22 等はどのサンプル点 (節点 0.1k / 中心
0.05+0.1k) からも 0.02 以上離れていて、許容誤差 (~1.7e-6) の影響を
受けない。

さらに:
 - 同じ領域を直方体 1 個 (shape=1, getspan 経路) と三角柱 2 個
   (shape=33, 全域ループの汎用経路。対角線で 2 分割した和集合) の
   両方で焼き付け、6 配列が完全一致することで両経路の整合を見る。
 - setupId_surface() は曲面近似の補正で、E 辺の周囲 4 つの H 位置に
   PEC (または分散性材料) があれば E を昇格させる。PEC 直方体で
   昇格される集合を検証済みの H 集合から独立に導出して突き合わせる。
   通常の誘電体 (type 1) では何も起きないこと (no-op) も見る。

なお boundingbox() は全形状について計算されるが、setupId() が結果を
使うのは shape = 1 の getspan 経路だけで、他の形状では捨てられる
(ループは全域を回り ingeometry() が判定する)。そのため shape != 1 の
bounding box の間違いはこのテストでは (そして実行でも) 観測できない。
*/

#include "ofd.h"
#include "ofd_test.h"

#include <stdlib.h>
#include <string.h>

void setupId(void);
void setupId_surface(void);

#define N     (10)      /* Nx = Ny = Nz */
#define DELTA (0.1)

/* 直方体 x:[0.22,0.78] y:[0.42,0.68] z:[0.12,0.48] に入る格子点。
   inc* = セル中心 (0.05+0.1i)、inn* = 節点 (0.1i) */
static int incx(int i) { return (i >= 2) && (i <= 7); }   /* 0.25 .. 0.75 */
static int innx(int i) { return (i >= 3) && (i <= 7); }   /* 0.30 .. 0.70 */
static int incy(int j) { return (j >= 4) && (j <= 6); }   /* 0.45 .. 0.65 */
static int inny(int j) { return (j >= 5) && (j <= 6); }   /* 0.50 .. 0.60 */
static int incz(int k) { return (k >= 1) && (k <= 4); }   /* 0.15 .. 0.45 */
static int innz(int k) { return (k >= 2) && (k <= 4); }   /* 0.20 .. 0.40 */

/* PEC 表面補正で使う「セル中心の帯を 1 つ外側 (負方向の隣接参照の分)
   まで広げた」範囲 : [c1, c2+1] */
static int bandx(int i) { return (i >= 2) && (i <= 8); }
static int bandy(int j) { return (j >= 4) && (j <= 7); }
static int bandz(int k) { return (k >= 1) && (k <= 5); }

static id_t *ref[6];   /* 直方体で焼いた 6 配列の控え (経路比較用) */

static void setup_grid(void)
{
	Nx = Ny = Nz = N;
	iMin = jMin = kMin = 0;
	iMax = jMax = kMax = N;

	Xn = (double *)malloc((N + 1) * sizeof(double));
	Yn = (double *)malloc((N + 1) * sizeof(double));
	Zn = (double *)malloc((N + 1) * sizeof(double));
	Xc = (double *)malloc(N * sizeof(double));
	Yc = (double *)malloc(N * sizeof(double));
	Zc = (double *)malloc(N * sizeof(double));
	for (int i = 0; i <= N; i++) {
		Xn[i] = Yn[i] = Zn[i] = DELTA * i;
	}
	for (int i = 0; i < N; i++) {
		Xc[i] = Yc[i] = Zc[i] = DELTA * (i + 0.5);
	}

	/* 配列レイアウトは setupSize() と同じ式 (Mur 相当の余白 1 層)。
	   setupSize() 自体は setupABCsize() 経由で Mur/PML のセットアップ
	   関数群に依存していてテストバイナリに入れられないため、
	   式を直接使う。assert 条件 NA(iMin-1,...) == 0 も同じ。 */
	Nk = 1;
	Nj = (kMax - kMin) + 2 + 1;
	Ni = ((jMax - jMin) + 2 + 1) * Nj;
	N0 = -((iMin - 1) * Ni + (jMin - 1) * Nj + (kMin - 1) * Nk);
	NN = NA(iMax + 1, jMax + 1, kMax + 1) + 1;
	CHECK(NA(iMin - 1, jMin - 1, kMin - 1) == 0);

	iEx = (id_t *)malloc(NN * sizeof(id_t));
	iEy = (id_t *)malloc(NN * sizeof(id_t));
	iEz = (id_t *)malloc(NN * sizeof(id_t));
	iHx = (id_t *)malloc(NN * sizeof(id_t));
	iHy = (id_t *)malloc(NN * sizeof(id_t));
	iHz = (id_t *)malloc(NN * sizeof(id_t));
	for (int i = 0; i < 6; i++) {
		ref[i] = (id_t *)malloc(NN * sizeof(id_t));
	}

	/* 材料 : 0=真空, 1=PEC, 2=ユーザー材料 (誘電体 or 分散性) */
	NMaterial = 3;
	Material = (material_t *)calloc(NMaterial, sizeof(material_t));
	Material[0].type = 1; Material[0].epsr = 1; Material[0].amur = 1;
	Material[1].type = 1; Material[1].epsr = 1; Material[1].amur = 1;
	Material[2].type = 1; Material[2].epsr = 4; Material[2].amur = 1;

	Geometry = (geometry_t *)calloc(2, sizeof(geometry_t));
}

static void clear_id(void)
{
	memset(iEx, 0, NN * sizeof(id_t));
	memset(iEy, 0, NN * sizeof(id_t));
	memset(iEz, 0, NN * sizeof(id_t));
	memset(iHx, 0, NN * sizeof(id_t));
	memset(iHy, 0, NN * sizeof(id_t));
	memset(iHz, 0, NN * sizeof(id_t));
}

static void set_box(id_t m)
{
	NGeometry = 1;
	Geometry[0].m = m;
	Geometry[0].shape = 1;
	Geometry[0].g[0] = 0.22; Geometry[0].g[1] = 0.78;
	Geometry[0].g[2] = 0.42; Geometry[0].g[3] = 0.68;
	Geometry[0].g[4] = 0.12; Geometry[0].g[5] = 0.48;
}

/* 6 成分の期待集合との突き合わせ。
   ループ範囲は setupId() が書き込む格子範囲そのもの。
   不一致は 1 成分 1 カウントに集約する (半セルずれのバグは数百点が
   一斉にずれるので、CHECK を点ごとに出すとログが溢れる)。 */
static void check_box_ids(id_t m)
{
	int bad[6] = {0, 0, 0, 0, 0, 0};
	int i, j, k;

	for (i = 0; i < N; i++)  for (j = 0; j <= N; j++) for (k = 0; k <= N; k++)
		if (IEX(i, j, k) != ((incx(i) && inny(j) && innz(k)) ? m : 0)) bad[0]++;
	for (i = 0; i <= N; i++) for (j = 0; j < N; j++)  for (k = 0; k <= N; k++)
		if (IEY(i, j, k) != ((innx(i) && incy(j) && innz(k)) ? m : 0)) bad[1]++;
	for (i = 0; i <= N; i++) for (j = 0; j <= N; j++) for (k = 0; k < N; k++)
		if (IEZ(i, j, k) != ((innx(i) && inny(j) && incz(k)) ? m : 0)) bad[2]++;
	for (i = 0; i <= N; i++) for (j = 0; j < N; j++)  for (k = 0; k < N; k++)
		if (IHX(i, j, k) != ((innx(i) && incy(j) && incz(k)) ? m : 0)) bad[3]++;
	for (i = 0; i < N; i++)  for (j = 0; j <= N; j++) for (k = 0; k < N; k++)
		if (IHY(i, j, k) != ((incx(i) && inny(j) && incz(k)) ? m : 0)) bad[4]++;
	for (i = 0; i < N; i++)  for (j = 0; j < N; j++)  for (k = 0; k <= N; k++)
		if (IHZ(i, j, k) != ((incx(i) && incy(j) && innz(k)) ? m : 0)) bad[5]++;

	CHECK(bad[0] == 0);   /* Ex : (Xc, Yn, Zn) */
	CHECK(bad[1] == 0);   /* Ey : (Xn, Yc, Zn) */
	CHECK(bad[2] == 0);   /* Ez : (Xn, Yn, Zc) */
	CHECK(bad[3] == 0);   /* Hx : (Xn, Yc, Zc) */
	CHECK(bad[4] == 0);   /* Hy : (Xc, Yn, Zc) */
	CHECK(bad[5] == 0);   /* Hz : (Xc, Yc, Zn) */
}

/* 1. 直方体の 6 成分スタガリング */
static void test_box_staggering(void)
{
	set_box(2);
	clear_id();
	setupId();

	/* 代表点 : 節点と中心の食い違いが出る境界を明示しておく */
	CHECK(IEX(2, 5, 3) == 2);   /* Xc[2] = 0.25 は箱の中 */
	CHECK(IEX(1, 5, 3) == 0);   /* Xc[1] = 0.15 は箱の外 */
	CHECK(IEY(2, 5, 3) == 0);   /* Xn[2] = 0.20 は箱の外 (Ex とは i の意味が違う) */
	CHECK(IEY(3, 5, 3) == 2);   /* Xn[3] = 0.30 は箱の中 */
	CHECK(IHX(5, 4, 1) == 2);   /* Yc[4] = 0.45, Zc[1] = 0.15 はどちらも箱の中 */
	CHECK(IHZ(5, 4, 1) == 0);   /* 同じ添字でも Zn[1] = 0.10 は箱の外 */

	check_box_ids(2);

	/* 控えを取る (2. の経路比較に使う) */
	memcpy(ref[0], iEx, NN * sizeof(id_t));
	memcpy(ref[1], iEy, NN * sizeof(id_t));
	memcpy(ref[2], iEz, NN * sizeof(id_t));
	memcpy(ref[3], iHx, NN * sizeof(id_t));
	memcpy(ref[4], iHy, NN * sizeof(id_t));
	memcpy(ref[5], iHz, NN * sizeof(id_t));

	/* 通常の誘電体 (type 1) では表面補正は何もしないこと。
	   highest() が昇格させるのは PEC と分散性材料 (type 2) だけで、
	   誘電体の境界は元の焼き付けのまま残るのが仕様 */
	setupId_surface();
	check_box_ids(2);
}

/* 1b. ドメイン端 (+X 面) に接する直方体 :
       shape = 1 の getspan 経路はループ範囲を絞ってから焼くので、
       探索範囲の上限を切り詰める間違いは端に接する物体の
       焼き付け漏れとして現れる (内側の物体では観測できない)。 */
static void test_edge_box_staggering(void)
{
	NGeometry = 1;
	Geometry[0].m = 2;
	Geometry[0].shape = 1;
	Geometry[0].g[0] = 0.62; Geometry[0].g[1] = 1.0;    /* +X 面に接する */
	Geometry[0].g[2] = 0.42; Geometry[0].g[3] = 0.68;
	Geometry[0].g[4] = 0.12; Geometry[0].g[5] = 0.48;

	clear_id();
	setupId();

	/* x : セル中心 0.65..0.95 -> i=6..9、節点 0.70..1.00 -> i=7..10
	   (境界 1.0 = Xn[Nx] は「境界は内側」の判定で含まれる) */
	{
		int bad[6] = {0, 0, 0, 0, 0, 0};
		int i, j, k;
		for (i = 0; i < N; i++)  for (j = 0; j <= N; j++) for (k = 0; k <= N; k++)
			if (IEX(i, j, k) != (((i >= 6) && (inny(j) && innz(k))) ? 2 : 0)) bad[0]++;
		for (i = 0; i <= N; i++) for (j = 0; j < N; j++)  for (k = 0; k <= N; k++)
			if (IEY(i, j, k) != (((i >= 7) && (incy(j) && innz(k))) ? 2 : 0)) bad[1]++;
		for (i = 0; i <= N; i++) for (j = 0; j <= N; j++) for (k = 0; k < N; k++)
			if (IEZ(i, j, k) != (((i >= 7) && (inny(j) && incz(k))) ? 2 : 0)) bad[2]++;
		for (i = 0; i <= N; i++) for (j = 0; j < N; j++)  for (k = 0; k < N; k++)
			if (IHX(i, j, k) != (((i >= 7) && (incy(j) && incz(k))) ? 2 : 0)) bad[3]++;
		for (i = 0; i < N; i++)  for (j = 0; j <= N; j++) for (k = 0; k < N; k++)
			if (IHY(i, j, k) != (((i >= 6) && (inny(j) && incz(k))) ? 2 : 0)) bad[4]++;
		for (i = 0; i < N; i++)  for (j = 0; j < N; j++)  for (k = 0; k <= N; k++)
			if (IHZ(i, j, k) != (((i >= 6) && (incy(j) && innz(k))) ? 2 : 0)) bad[5]++;
		CHECK(bad[0] == 0);   /* Ex は端のセル i=9 まで焼かれること */
		CHECK(bad[1] == 0);   /* Ey は端の節点 i=10 まで焼かれること */
		CHECK(bad[2] == 0);
		CHECK(bad[3] == 0);
		CHECK(bad[4] == 0);
		CHECK(bad[5] == 0);
	}
}

/* 2. shape=1 の getspan 高速経路と汎用経路の整合 :
      同じ直方体を対角線で 2 分割した Z 三角柱 2 本 (shape=33, 汎用経路)
      で焼き、直方体 (shape=1, getspan 経路) の結果と完全一致すること */
static void test_pillar_union_equals_box(void)
{
	NGeometry = 2;
	/* 三角形 (0.22,0.42)-(0.78,0.42)-(0.22,0.68) : 斜辺は対角線 */
	Geometry[0].m = 2;
	Geometry[0].shape = 33;
	Geometry[0].g[0] = 0.12; Geometry[0].g[1] = 0.48;
	Geometry[0].g[2] = 0.22; Geometry[0].g[3] = 0.78; Geometry[0].g[4] = 0.22;
	Geometry[0].g[5] = 0.42; Geometry[0].g[6] = 0.42; Geometry[0].g[7] = 0.68;
	/* 三角形 (0.78,0.68)-(0.78,0.42)-(0.22,0.68) : 残り半分 */
	Geometry[1].m = 2;
	Geometry[1].shape = 33;
	Geometry[1].g[0] = 0.12; Geometry[1].g[1] = 0.48;
	Geometry[1].g[2] = 0.78; Geometry[1].g[3] = 0.78; Geometry[1].g[4] = 0.22;
	Geometry[1].g[5] = 0.68; Geometry[1].g[6] = 0.42; Geometry[1].g[7] = 0.68;

	clear_id();
	setupId();

	{
		int bad[6] = {0, 0, 0, 0, 0, 0};
		for (int64_t n = 0; n < NN; n++) {
			if (iEx[n] != ref[0][n]) bad[0]++;
			if (iEy[n] != ref[1][n]) bad[1]++;
			if (iEz[n] != ref[2][n]) bad[2]++;
			if (iHx[n] != ref[3][n]) bad[3]++;
			if (iHy[n] != ref[4][n]) bad[4]++;
			if (iHz[n] != ref[5][n]) bad[5]++;
		}
		CHECK(bad[0] == 0);
		CHECK(bad[1] == 0);
		CHECK(bad[2] == 0);
		CHECK(bad[3] == 0);
		CHECK(bad[4] == 0);
		CHECK(bad[5] == 0);
	}
}

/* 3. PEC 直方体の表面補正 :
      E 辺は自分の位置の判定に加え、周囲 4 つの H 位置のどれかが PEC なら
      PEC に昇格する。期待集合は 1. で検証済みの H 集合から導出した
      (実装と同じ式を呼び直すのではなく、集合演算で独立に求めている)。

      例えば Ex(i,j,k) の周囲は Hy(i,j,k), Hy(i,j,k-1), Hz(i,j,k), Hz(i,j-1,k)。
      Hy 集合 = incx&&inny&&incz を k-1 の分だけ広げると bandz、
      Hz 集合 = incx&&incy&&innz を j-1 の分だけ広げると bandy なので

          Ex が PEC  <=>  incx(i) && ((inny(j) && bandz(k)) || (bandy(j) && innz(k)))

      になる (元の Ex 集合はこの第 1 項に含まれる)。Ey, Ez は巡回対称。 */
static void test_pec_surface(void)
{
	set_box(PEC);
	clear_id();
	setupId();
	setupId_surface();

	/* 代表点 : 面のすぐ外は昇格、辺の斜め外と 2 つ外は昇格しない */
	CHECK(IEX(4, 4, 3) == PEC);   /* -Y 面のすぐ外 (Hz(4,4,3) が PEC) */
	CHECK(IEX(4, 7, 3) == PEC);   /* +Y 面のすぐ外 (Hz(4,6,3) が PEC) */
	CHECK(IEX(4, 5, 1) == PEC);   /* -Z 面のすぐ外 (Hy(4,5,1) が PEC) */
	CHECK(IEX(4, 5, 5) == PEC);   /* +Z 面のすぐ外 (Hy(4,5,4) が PEC) */
	CHECK(IEX(4, 4, 1) == 0);     /* 辺の斜め外 -> 昇格しない */
	CHECK(IEX(4, 3, 3) == 0);     /* 2 つ外 -> 昇格しない */

	{
		int bad[3] = {0, 0, 0};
		int i, j, k;
		for (i = 0; i < N; i++)  for (j = 0; j <= N; j++) for (k = 0; k <= N; k++) {
			const id_t e = (incx(i) && ((inny(j) && bandz(k)) || (bandy(j) && innz(k)))) ? PEC : 0;
			if (IEX(i, j, k) != e) bad[0]++;
		}
		for (i = 0; i <= N; i++) for (j = 0; j < N; j++)  for (k = 0; k <= N; k++) {
			const id_t e = (incy(j) && ((innz(k) && bandx(i)) || (bandz(k) && innx(i)))) ? PEC : 0;
			if (IEY(i, j, k) != e) bad[1]++;
		}
		for (i = 0; i <= N; i++) for (j = 0; j <= N; j++) for (k = 0; k < N; k++) {
			const id_t e = (incz(k) && ((innx(i) && bandy(j)) || (bandx(i) && inny(j)))) ? PEC : 0;
			if (IEZ(i, j, k) != e) bad[2]++;
		}
		CHECK(bad[0] == 0);
		CHECK(bad[1] == 0);
		CHECK(bad[2] == 0);
	}

	/* H は表面補正で変化しないこと */
	{
		int bad = 0;
		int i, j, k;
		for (i = 0; i <= N; i++) for (j = 0; j < N; j++) for (k = 0; k < N; k++)
			if (IHX(i, j, k) != ((innx(i) && incy(j) && incz(k)) ? PEC : 0)) bad++;
		for (i = 0; i < N; i++)  for (j = 0; j <= N; j++) for (k = 0; k < N; k++)
			if (IHY(i, j, k) != ((incx(i) && inny(j) && incz(k)) ? PEC : 0)) bad++;
		for (i = 0; i < N; i++)  for (j = 0; j < N; j++)  for (k = 0; k <= N; k++)
			if (IHZ(i, j, k) != ((incx(i) && incy(j) && innz(k)) ? PEC : 0)) bad++;
		CHECK(bad == 0);
	}
}

/* 4. 分散性材料 (type 2) も表面補正で昇格すること。
      highest() は 4 つの近傍を個別の引数で受けるので、4 近傍を
      **1 つずつ単独で** 分散性にできる点を選び、どの近傍の判定が
      欠けても検出できるようにする (Ex の近傍は順に
      Hy(i,j,k), Hy(i,j,k-1), Hz(i,j,k), Hz(i,j-1,k))。 */
static void test_dispersive_surface(void)
{
	Material[2].type = 2;   /* 同じ id 2 を分散性にする */

	set_box(2);
	clear_id();
	setupId();
	setupId_surface();

	/* 各点で分散性なのは指定した 1 近傍だけ (他の 3 つは真空) */
	CHECK(IEX(4, 5, 1) == 2);   /* Hy(i,j,k)   = Hy(4,5,1) 経由 */
	CHECK(IEX(4, 5, 5) == 2);   /* Hy(i,j,k-1) = Hy(4,5,4) 経由 */
	CHECK(IEX(4, 4, 3) == 2);   /* Hz(i,j,k)   = Hz(4,4,3) 経由 */
	CHECK(IEX(4, 7, 3) == 2);   /* Hz(i,j-1,k) = Hz(4,6,3) 経由 */

	CHECK(IEX(4, 3, 3) == 0);   /* 2 つ外 -> 昇格しない */
	CHECK(IEX(1, 5, 3) == 0);   /* 箱の外側の帯以外は変化しない */
	CHECK(IEX(4, 5, 3) == 2);   /* 内部はそのまま */

	Material[2].type = 1;   /* 後片付け */
}

void test_setupid(void)
{
	setup_grid();

	test_box_staggering();
	test_edge_box_staggering();
	test_pillar_union_equals_box();
	test_pec_surface();
	test_dispersive_surface();

	/* 後片付け : 他のテストがグローバルを再利用しても壊れないように */
	free(Xn); free(Yn); free(Zn);
	free(Xc); free(Yc); free(Zc);
	free(iEx); free(iEy); free(iEz);
	free(iHx); free(iHy); free(iHz);
	for (int i = 0; i < 6; i++) {
		free(ref[i]);
	}
	free(Material);
	free(Geometry);
	NGeometry = 0;
	NMaterial = 0;
}
