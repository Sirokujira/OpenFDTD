/*
test_vector.c

vector 経路の係数展開 (sol/setup_vector.c) が no-vector 経路と同じ係数を
作ることのテスト。

更新式は 2 系統ある (README「ベクトル化更新経路」):
  no-vector : 点ごとに材料 id を引き、C1[m] / C2[m] (H は D1 / D2) を掛ける
              (材料境界の id は setupId_surface() が決める)
  vector    : 全点の係数を K1 / K2 配列に展開しておく
              (境界の id は setupId_surface_vector() が別の実装で決める)
境界の決め方が 2 か所に別々に書かれているので、同じ形状に対して全点の
実効係数 (K1[n] と C1[id[n]] など) が一致することを確かめる。
data/sample/vector_check.sh は 2 経路の計算結果 (ログと HDF5) の一致を
通しで見るが、こちらはどの点で食い違ったかまで指す。

形状は「誘電体の箱に PEC の板と分散性の箱が接し、別の誘電体とも接する」
ものにして、PEC 優先・分散性優先・誘電体どうしの 3 種の境界規則を踏む。
分散性の点のリスト (setupDispersion) も両経路で同じ数になること。
*/

#include "ofd.h"
#include "ofd_prototype.h"
#include "ofd_test.h"

#include <stdlib.h>
#include <string.h>

#define NC 8

static void vc_setup(void)
{
	Nx = Ny = Nz = NC;
	iMin = jMin = kMin = 0;
	iMax = jMax = kMax = NC;
	double **nd[3] = {&Xn, &Yn, &Zn}, **ce[3] = {&Xc, &Yc, &Zc};
	for (int a = 0; a < 3; a++) {
		*nd[a] = (double *)malloc((NC + 1) * sizeof(double));
		*ce[a] = (double *)malloc(NC * sizeof(double));
		(*nd[a])[0] = 0;
		for (int i = 0; i < NC; i++) (*nd[a])[i + 1] = (*nd[a])[i] + 1e-3 * (1 + 0.1 * ((i + a) % 3));
		for (int i = 0; i < NC; i++) (*ce[a])[i] = ((*nd[a])[i] + (*nd[a])[i + 1]) / 2;
	}

	iABC = 0;
	NFeed = 1;
	IPlanewave = 0;
	NFreq1 = NFreq2 = 0;
	Dt = 1e-12;
	Tw = 1e-10;

	Nk = 1;
	Nj = (kMax - kMin) + 2 + 1;
	Ni = ((jMax - jMin) + 2 + 1) * Nj;
	N0 = -((iMin - 1) * Ni + (jMin - 1) * Nj + (kMin - 1) * Nk);
	NN = NA(iMax + 1, jMax + 1, kMax + 1) + 1;

	id_t **ids[6] = {&iEx, &iEy, &iEz, &iHx, &iHy, &iHz};
	real_t **k1[6] = {&K1Ex, &K1Ey, &K1Ez, &K1Hx, &K1Hy, &K1Hz};
	real_t **k2[6] = {&K2Ex, &K2Ey, &K2Ez, &K2Hx, &K2Hy, &K2Hz};
	for (int c = 0; c < 6; c++) {
		*ids[c] = (id_t *)calloc(NN, sizeof(id_t));
		*k1[c] = (real_t *)calloc(NN, sizeof(real_t));
		*k2[c] = (real_t *)calloc(NN, sizeof(real_t));
	}

	/* 0=真空, 1=PEC, 2=損失誘電体, 3=分散性, 4=磁性誘電体 */
	NMaterial = 5;
	Material = (material_t *)calloc(NMaterial, sizeof(material_t));
	for (int m = 0; m < NMaterial; m++) { Material[m].type = 1; Material[m].epsr = Material[m].amur = 1; }
	Material[2].epsr = 4;   Material[2].esgm = 0.3;
	Material[3].type = 2;   Material[3].einf = 2; Material[3].be = 4e10; Material[3].ce = 1e10;
	Material[4].epsr = 2.5; Material[4].amur = 1.8; Material[4].msgm = 500;
	C1 = (real_t *)calloc(NMaterial, sizeof(real_t));
	C2 = (real_t *)calloc(NMaterial, sizeof(real_t));
	D1 = (real_t *)calloc(NMaterial, sizeof(real_t));
	D2 = (real_t *)calloc(NMaterial, sizeof(real_t));
	real_t **rn[3] = {&RXn, &RYn, &RZn}, **rc[3] = {&RXc, &RYc, &RZc};
	for (int a = 0; a < 3; a++) {
		*rn[a] = (real_t *)calloc(NC + 1, sizeof(real_t));
		*rc[a] = (real_t *)calloc(NC, sizeof(real_t));
	}

	/* 誘電体の箱、その面に接する PEC の板、隣の分散性の箱、別の誘電体の箱
	   (座標は節点に合わせた値。境界は節点・セル中心の両方にまたがる) */
	/* 境界の点の材料は、形状を焼く順で後のものが勝つ。分散性優先の規則が
	   効くのは「面を共有する後の誘電体が面の E を上書きし、隣の H だけが
	   分散性」の点なので、分散性の箱を先に置き、誘電体の箱 (後) と各方向で
	   面を共有させる。周囲 4 つの H のどれが分散性かが方向ごとに変わるので、
	   規則の各分岐 (highest の t1..t4) を踏む */
	NGeometry = 10;
	Geometry = (geometry_t *)calloc(NGeometry, sizeof(geometry_t));
	const double g[10][6] = {
		{Xn[0], Xn[1], Yn[3], Yn[5], Zn[1], Zn[3]},   /* 1 : PEC の箱 (先に置く) */
		{Xn[5], Xn[7], Yn[3], Yn[5], Zn[1], Zn[4]},   /* 3 : 分散性 (誘電体の +x 側) */
		{Xn[3], Xn[4], Yn[1], Yn[2], Zn[2], Zn[3]},   /* 3 : 分散性 (-y 側) */
		{Xn[3], Xn[4], Yn[3], Yn[4], Zn[0], Zn[1]},   /* 3 : 分散性 (-z 側) */
		{Xn[3], Xn[4], Yn[6], Yn[7], Zn[2], Zn[3]},   /* 3 : 分散性 (+y 側) */
		{Xn[6], Xn[7], Yn[0], Yn[1], Zn[5], Zn[6]},   /* 3 : 分散性 (y = 0..1 : 添字 1 の境界) */
		{Xn[2], Xn[5], Yn[2], Yn[6], Zn[1], Zn[4]},   /* 2 : 誘電体 */
		{Xn[6], Xn[7], Yn[1], Yn[2], Zn[5], Zn[6]},   /* 4 : 磁性誘電体 (y = 1 で分散性と接する) */
		{Xn[1], Xn[2], Yn[2], Yn[6], Zn[1], Zn[6]},   /* 4 : 磁性誘電体 (PEC の箱と x = 1 で接する) */
		{Xn[2], Xn[5], Yn[2], Yn[6], Zn[4], Zn[4]},   /* 1 : PEC の板 (厚さ 0) */
	};
	/* PEC の箱も先に置き、後の誘電体が面の E を上書きした点で PEC 優先の
	   規則 (周囲の H に PEC があれば PEC) が効くようにする */
	const id_t gm[10] = {1, 3, 3, 3, 3, 3, 2, 4, 4, 1};
	for (int n = 0; n < NGeometry; n++) {
		Geometry[n].m = gm[n];
		Geometry[n].shape = 1;
		memcpy(Geometry[n].g, g[n], sizeof(g[n]));
	}
}

static void vc_teardown(void)
{
	free(Xn); free(Yn); free(Zn); free(Xc); free(Yc); free(Zc);
	Xn = Yn = Zn = Xc = Yc = Zc = NULL;
	id_t **ids[6] = {&iEx, &iEy, &iEz, &iHx, &iHy, &iHz};
	real_t **k1[6] = {&K1Ex, &K1Ey, &K1Ez, &K1Hx, &K1Hy, &K1Hz};
	real_t **k2[6] = {&K2Ex, &K2Ey, &K2Ez, &K2Hx, &K2Hy, &K2Hz};
	for (int c = 0; c < 6; c++) {
		free(*ids[c]); free(*k1[c]); free(*k2[c]);
		*ids[c] = NULL; *k1[c] = *k2[c] = NULL;
	}
	free(Material); free(C1); free(C2); free(D1); free(D2);
	Material = NULL; C1 = C2 = D1 = D2 = NULL;
	free(RXn); free(RYn); free(RZn); free(RXc); free(RYc); free(RZc);
	RXn = RYn = RZn = RXc = RYc = RZc = NULL;
	free(Geometry); Geometry = NULL; NGeometry = 0;
	NMaterial = 0;
	NFeed = 0;
	VECTOR = 0;
}

static void run_setup(int vector)
{
	VECTOR = vector;
	numMurHx = numMurHy = numMurHz = 0;
	setupMurHx(0); setupMurHy(0); setupMurHz(0);
	fMurHx = (mur_t *)calloc(numMurHx, sizeof(mur_t));
	fMurHy = (mur_t *)calloc(numMurHy, sizeof(mur_t));
	fMurHz = (mur_t *)calloc(numMurHz, sizeof(mur_t));
	setup();
	free(fMurHx); free(fMurHy); free(fMurHz);
	fMurHx = fMurHy = fMurHz = NULL;
}

static void free_dispersion(void)
{
	free(mDispersionEx); free(mDispersionEy); free(mDispersionEz);
	free(DispersionEx); free(DispersionEy); free(DispersionEz);
	mDispersionEx = mDispersionEy = mDispersionEz = NULL;
	DispersionEx = DispersionEy = DispersionEz = NULL;
}

/* 成分 c (0..2 = Ex..Ez, 3..5 = Hx..Hz) の更新範囲内か */
static int in_range(int c, int i, int j, int k)
{
	const int a = c % 3, isE = (c < 3);
	const int idx[3] = {i, j, k};
	for (int b = 0; b < 3; b++) {
		const int on_cell = isE ? (b == a) : (b != a);
		if ((idx[b] < 0) || (idx[b] > (on_cell ? NC - 1 : NC))) return 0;
	}
	return 1;
}

void test_vector(void)
{
	const double dt0 = Dt, tw0 = Tw;
	vc_setup();

	id_t **ids[6] = {&iEx, &iEy, &iEz, &iHx, &iHy, &iHz};
	real_t **k1[6] = {&K1Ex, &K1Ey, &K1Ez, &K1Hx, &K1Hy, &K1Hz};
	real_t **k2[6] = {&K2Ex, &K2Ey, &K2Ez, &K2Hx, &K2Hy, &K2Hz};

	/* no-vector : 実効係数を控える */
	run_setup(0);
	real_t *e1[6], *e2[6];
	int64_t nd0[3] = {numDispersionEx, numDispersionEy, numDispersionEz};
	int nmat[5] = {0, 0, 0, 0, 0};
	for (int c = 0; c < 6; c++) {
		e1[c] = (real_t *)calloc(NN, sizeof(real_t));
		e2[c] = (real_t *)calloc(NN, sizeof(real_t));
		const real_t *A1 = (c < 3) ? C1 : D1, *A2 = (c < 3) ? C2 : D2;
		for (int64_t q = 0; q < NN; q++) {
			const id_t m = (*ids[c])[q];
			e1[c][q] = A1[m];
			e2[c][q] = A2[m];
		}
	}
	for (int i = 0; i <= NC; i++) for (int j = 0; j <= NC; j++) for (int k = 0; k <= NC; k++) {
		if (in_range(0, i, j, k)) nmat[IEX(i, j, k)]++;
	}
	free_dispersion();

	/* 形状が境界規則を踏んでいること (全材料が実際に現れる) */
	CHECK(nmat[1] > 0 && nmat[2] > 0 && nmat[3] > 0 && nmat[4] > 0);
	CHECK(nd0[0] > 0 && nd0[1] > 0 && nd0[2] > 0);

	/* vector : K1/K2 が同じ値か (更新範囲の全点) */
	run_setup(1);
	CHECK(numDispersionEx == nd0[0] && numDispersionEy == nd0[1] && numDispersionEz == nd0[2]);
	for (int c = 0; c < 6; c++) {
		int bad = 0, first = -1;
		for (int i = -1; i <= NC + 1; i++) {
		for (int j = -1; j <= NC + 1; j++) {
		for (int k = -1; k <= NC + 1; k++) {
			if (!in_range(c, i, j, k)) continue;
			const int64_t q = NA(i, j, k);
			if (((*k1[c])[q] != e1[c][q]) || ((*k2[c])[q] != e2[c][q])) {
				if (first < 0) first = (int)q;
				bad++;
			}
		}
		}
		}
		CHECK(bad == 0);
		if (bad) printf("    component %d : %d points differ (first index %d)\n", c, bad, first);
		free(e1[c]); free(e2[c]);
	}
	free_dispersion();

	vc_teardown();
	Dt = dt0;
	Tw = tw0;
}
