/*
test_readhdf5.c

HDF5 の書き手 (sol/outputHdf5.c) と読み手 readhdf5() (post/readhdf5.c) の往復。

readhdf5() は ofd.out (readout) から HDF5 へ移行するための読み手で、既定では
呼ばれない (src/post_Main.c の NOTE。coupling / 断面積の展開が未実装)。
既定で呼ばれないため、書き手側でパスや型を変えても気付けない
(実際に /geometry や /convergence へ移したときに長く壊れていた、
data/sample/hdf5_layout_check.sh の冒頭)。ここでは書き手が出す項目のうち
readhdf5() が読む項目 — 格子の大きさと配列レイアウト、Dt、平面波の
パラメータ、格子座標、収束履歴、給電の時系列、周波数、形状の線分と材料 id、
等価面、Title — が、書いた値のまま戻ることを確かめる。
*/

#include "ofd.h"
#include "ofd_hdf5.h"
#include "ofd_test.h"
#include "hdf5.h"

#include <stdlib.h>
#include <string.h>

int readhdf5(void);

static double *dup(const double *a, int n)
{
	double *b = (double *)malloc((size_t)n * sizeof(double));
	memcpy(b, a, (size_t)n * sizeof(double));
	return b;
}

static int same(const double *a, const double *b, int n)
{
	if ((a == NULL) || (b == NULL)) return 0;
	for (int i = 0; i < n; i++) if (a[i] != b[i]) return 0;
	return 1;
}

void test_readhdf5(void)
{
	char title0[256];
	memcpy(title0, Title, sizeof(title0));
	const int out0 = Hdf5Output, int0 = Hdf5Interval;

	/* ---- 書く値 ---- */
	enum { NX = 3, NY = 2, NZ = 4, MAXIT = 6, NITER = 3, NGL = 2, NSURF = 2 };
	const double xn[NX + 1] = {0.0, 0.1, 0.25, 0.4}, xc[NX] = {0.05, 0.175, 0.325};
	const double yn[NY + 1] = {-0.2, 0.0, 0.3},      yc[NY] = {-0.1, 0.15};
	const double zn[NZ + 1] = {1, 2, 3, 4, 5},       zc[NZ] = {1.5, 2.5, 3.5, 4.5};
	const double eit[NITER] = {1.0, 0.5, 0.125}, hit[NITER] = {3.0, 2.0, 1.0};
	double vfeed[MAXIT + 1], ifeed[MAXIT + 1];
	for (int n = 0; n <= MAXIT; n++) { vfeed[n] = 0.5 * n - 1; ifeed[n] = -0.25 * n * n; }
	const double freq2[2] = {1.5e9, 2.5e9};
	const double gline[NGL][2][3] = {{{0.1, 0.2, 0.3}, {0.4, 0.5, 0.6}}, {{-1, -2, -3}, {7, 8, 9}}};
	const id_t mgline[NGL] = {PEC, 3};

	strcpy(Title, "readhdf5 round trip");
	Nx = NX; Ny = NY; Nz = NZ;
	iMin = jMin = kMin = 0; iMax = Nx; jMax = Ny; kMax = Nz;
	Nk = 1; Nj = (Nz + 3) * Nk; Ni = (Ny + 3) * Nj; N0 = Ni + Nj + Nk;
	NN = NA(Nx + 1, Ny + 1, Nz + 1) + 1;
	Xn = dup(xn, NX + 1); Yn = dup(yn, NY + 1); Zn = dup(zn, NZ + 1);
	Xc = dup(xc, NX);     Yc = dup(yc, NY);     Zc = dup(zc, NZ);
	NGline = NGL;
	Gline = (double (*)[2][3])malloc(sizeof(gline));
	memcpy(Gline, gline, sizeof(gline));
	MGline = (id_t *)malloc(sizeof(mgline));
	memcpy(MGline, mgline, sizeof(mgline));
	Solver.maxiter = MAXIT; Solver.nout = 2;
	Niter = NITER; Ntime = MAXIT;
	NFreq1 = 0; NFreq2 = 2; Freq2 = dup(freq2, 2); Freq1 = NULL;
	NFeed = 1; VFeed = dup(vfeed, MAXIT + 1); IFeed = dup(ifeed, MAXIT + 1);
	NPoint = 0; VPoint = NULL;
	IPlanewave = 1;
	Dt = 1.25e-12;
	Planewave.theta = 30; Planewave.phi = 45; Planewave.pol = 2;
	for (int m = 0; m < 3; m++) {
		Planewave.ei[m] = 0.1 * (m + 1); Planewave.hi[m] = -0.2 * (m + 1);
		Planewave.ri[m] = 0.3 * (m + 1); Planewave.r0[m] = 0.4 * (m + 1);
	}
	Planewave.ai = 7.5;
	NSurface = NSURF;
	Surface = (surface_t *)calloc(NSURF, sizeof(surface_t));
	Surface[0].nx = 1;  Surface[0].x = 0.4; Surface[0].y = 0.05; Surface[0].z = 2.5; Surface[0].ds = 0.02;
	Surface[1].nz = -1; Surface[1].x = 0.2; Surface[1].y = -0.1; Surface[1].z = 1.0; Surface[1].ds = 0.03;
	const surface_t s0 = Surface[0], s1 = Surface[1];
	const int64_t ni0 = Ni, nj0 = Nj, nk0 = Nk, n00 = N0, nn0 = NN;

	Hdf5Output = 1;
	Hdf5Interval = 0;
	remove(OFD_H5_FILE);
	CHECK(hdf5_open(1) == 1);
	hdf5_write_convergence(NITER, eit, hit);
	/* /post の近傍界 DFT : 余白込みの NN 全体 (余白にも値を入れる) */
	float *cf[12];
	for (int b = 0; b < 12; b++) {
		cf[b] = (float *)malloc((size_t)NFreq2 * NN * sizeof(float));
		for (int64_t n = 0; n < NFreq2 * NN; n++) cf[b][n] = (float)(b * 100000 + n % 99991) * 0.5f;
	}
	hdf5_write_freqdomain(cf[0], cf[1], cf[2], cf[3], cf[4], cf[5], cf[6], cf[7], cf[8], cf[9], cf[10], cf[11]);
	hdf5_close();
	/* /post の給電と入力電力 (outputChars の後に追記される) */
	Feed = (feed_t *)calloc(1, sizeof(feed_t));
	Feed[0].dir = 'Y'; Feed[0].i = 2; Feed[0].j = 1; Feed[0].k = 3;
	Feed[0].dx = 0.1; Feed[0].dy = 0.2; Feed[0].dz = 0.3;
	Feed[0].volt = 1.5; Feed[0].delay = 2e-10; Feed[0].z0 = 37.5;   /* 50 でも 73 でもない値 (読み手が定数を入れる変異を捕まえる) */
	const feed_t feed0 = Feed[0];
	const double pin0[2] = {0.75, 1.25}, pin1[2] = {0.8, 1.5};
	Pin[0] = dup(pin0, 2); Pin[1] = dup(pin1, 2);
	hdf5_write_post_feed();
	free(Feed); free(Pin[0]); free(Pin[1]);
	Feed = NULL; Pin[0] = Pin[1] = NULL;

	/* ---- 消す (読み手が確保し直すもの) ----
	   書いた配列はここでは解放しない : 解放すると readhdf5 の malloc が同じ領域を
	   (中身ごと) 再利用し、読み損ねても値が一致してしまう (等価面の ds を読まない
	   変異が実際にそれで通った)。比較の後で解放する */
	void *keep[12] = {Xn, Yn, Zn, Xc, Yc, Zc, Gline, MGline, Freq2, VFeed, IFeed, Surface};
	Xn = Yn = Zn = Xc = Yc = Zc = NULL;
	Gline = NULL; MGline = NULL; Freq2 = NULL; VFeed = IFeed = NULL; Surface = NULL;
	Eiter = Hiter = NULL;
	/* int64_t の項目は全ビットを立てておく : 4 バイト (int) で読むと上位が残る
	   (0 にしておくとリトルエンディアンでは下位 4 バイトだけで正しく見えてしまう) */
	Nx = Ny = Nz = 0; Ni = Nj = Nk = N0 = NN = -1;
	NGline = -1; NSurface = 0; NFreq2 = 0; NFeed = 0; Niter = 0; Ntime = 0;
	Solver.maxiter = Solver.nout = 0; IPlanewave = 0; Dt = 0;
	memset(&Planewave, 0, sizeof(Planewave));
	memset(Title, 0, sizeof(Title));

	/* ---- 読む ---- (書き手が出さない coupling 等の欠落で HDF5 がエラースタックを
	   表示するので、ここだけ止める) */
	H5E_auto2_t efunc;
	void *edata;
	H5Eget_auto2(H5E_DEFAULT, &efunc, &edata);
	H5Eset_auto2(H5E_DEFAULT, NULL, NULL);
	const int rret = readhdf5();
	H5Eset_auto2(H5E_DEFAULT, efunc, edata);

	CHECK((Nx == NX) && (Ny == NY) && (Nz == NZ));
	CHECK((Ni == ni0) && (Nj == nj0) && (Nk == nk0) && (N0 == n00) && (NN == nn0));
	CHECK((Solver.maxiter == MAXIT) && (Solver.nout == 2) && (Niter == NITER) && (Ntime == MAXIT));
	CHECK((NFreq2 == 2) && (NFeed == 1) && (IPlanewave == 1) && (NGline == NGL));
	CHECK(Dt == 1.25e-12);
	CHECK((Planewave.theta == 30) && (Planewave.phi == 45) && (Planewave.pol == 2) && (Planewave.ai == 7.5));
	CHECK((Planewave.ei[2] == 0.1 * 3) && (Planewave.hi[1] == -0.2 * 2) && (Planewave.ri[0] == 0.3) && (Planewave.r0[2] == 0.4 * 3));
	CHECK(same(Xn, xn, NX + 1) && same(Yn, yn, NY + 1) && same(Zn, zn, NZ + 1));
	CHECK(same(Xc, xc, NX) && same(Yc, yc, NY) && same(Zc, zc, NZ));
	CHECK(same(Eiter, eit, NITER) && same(Hiter, hit, NITER));
	CHECK(same(VFeed, vfeed, MAXIT + 1) && same(IFeed, ifeed, MAXIT + 1));
	CHECK(same(Freq2, freq2, 2));
	CHECK((Gline != NULL) && (memcmp(Gline, gline, sizeof(gline)) == 0));
	CHECK((MGline != NULL) && (MGline[0] == PEC) && (MGline[1] == 3));
	CHECK(NSurface == NSURF);
	CHECK((Surface != NULL) && (memcmp(&Surface[0], &s0, sizeof(surface_t)) == 0) && (memcmp(&Surface[1], &s1, sizeof(surface_t)) == 0));
	CHECK(strcmp(Title, "readhdf5 round trip") == 0);

	/* /post : 近傍界 DFT は余白込みで全要素、給電と入力電力 */
	CHECK(rret == 0);
	{
		float *dst[12] = {cEx_r, cEx_i, cEy_r, cEy_i, cEz_r, cEz_i, cHx_r, cHx_i, cHy_r, cHy_i, cHz_r, cHz_i};
		int ok = 1;
		for (int b = 0; b < 12; b++) {
			ok &= (dst[b] != NULL) && (memcmp(dst[b], cf[b], (size_t)NFreq2 * NN * sizeof(float)) == 0);
		}
		CHECK(ok);
	}
	CHECK((Feed != NULL) && (Feed[0].dir == 'Y') && (Feed[0].i == 2) && (Feed[0].j == 1) && (Feed[0].k == 3));
	CHECK((Feed != NULL) && (Feed[0].dz == feed0.dz) && (Feed[0].volt == feed0.volt) && (Feed[0].delay == feed0.delay) && (Feed[0].z0 == 37.5));
	CHECK(same(Pin[0], pin0, 2) && same(Pin[1], pin1, 2));
	/* 遠方界の等価面は ofd.out と同じ大きさで確保される (/metadata の NSurface は
	   ここでは 2 を書いたのでそれを読む。ソルバーが書く普通のファイルでは 0 で、
	   alloc_farfield の大きさになる : data/sample/hdf5_post_check.sh) */
	CHECK((SurfaceEx != NULL) && (SurfaceEx[0] != NULL) && (SurfaceHz[1] != NULL));
	for (int b = 0; b < 12; b++) free(cf[b]);
	if (strcmp(Title, "readhdf5 round trip") != 0) printf("    Title read back as \"%s\"\n", Title);

	for (int n = 0; n < 12; n++) free(keep[n]);

	/* ---- 後始末 (readhdf5 が確保したもの) ---- */
	free(Xn); free(Yn); free(Zn); free(Xc); free(Yc); free(Zc);
	free(Eiter); free(Hiter); free(VFeed); free(IFeed); free(Freq2);
	free(Gline); free(MGline); free(Surface);
	free(Feed); free(Zin); free(Ref); free(Pin[0]); free(Pin[1]); free(Spara);
	for (int f = 0; f < NFreq2; f++) {
		free(SurfaceEx[f]); free(SurfaceEy[f]); free(SurfaceEz[f]);
		free(SurfaceHx[f]); free(SurfaceHy[f]); free(SurfaceHz[f]);
	}
	free(SurfaceEx); free(SurfaceEy); free(SurfaceEz); free(SurfaceHx); free(SurfaceHy); free(SurfaceHz);
	free(cEx_r); free(cEx_i); free(cEy_r); free(cEy_i); free(cEz_r); free(cEz_i);
	free(cHx_r); free(cHx_i); free(cHy_r); free(cHy_i); free(cHz_r); free(cHz_i);
	Xn = Yn = Zn = Xc = Yc = Zc = NULL;
	Eiter = Hiter = VFeed = IFeed = Freq2 = NULL;
	Gline = NULL; MGline = NULL; Surface = NULL;
	Feed = NULL; Zin = NULL; Ref = NULL; Pin[0] = Pin[1] = NULL; Spara = NULL;
	SurfaceEx = SurfaceEy = SurfaceEz = SurfaceHx = SurfaceHy = SurfaceHz = NULL;
	cEx_r = cEx_i = cEy_r = cEy_i = cEz_r = cEz_i = NULL;
	cHx_r = cHx_i = cHy_r = cHy_i = cHz_r = cHz_i = NULL;
	remove(OFD_H5_FILE);
	Nx = Ny = Nz = 0; NN = 0; NGline = 0; NSurface = 0; NFreq2 = 0; NFeed = 0;
	Niter = 0; Ntime = 0; IPlanewave = 0;
	memset(&Planewave, 0, sizeof(Planewave));
	Hdf5Output = out0;
	Hdf5Interval = int0;
	memcpy(Title, title0, sizeof(title0));
}
