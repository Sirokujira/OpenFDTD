/*
test_outfile.c

ofd.out の書き手 writeout() (sol/writeout.c) と読み手 readout()
(post/readout.c, ofd_post) の往復テスト。

ofd.out はソルバーと後処理 ofd_post (別バイナリ) の間の契約で、書く順番・型・
個数が両側で一致していなければならない。片方にだけ項目を足すと以降が
すべてずれるが、ofd_post は "invalid file length" を stderr に出すだけで
作図を続ける。ここでは全項目に値を入れて書き、グローバルを消してから読み戻し、
すべて一致することを確かめる。

NGline / NSurface は int64_t の変数だが、ファイル上は 4 バイト (int)。以前は
変数を sizeof(int) で直接読み書きしていて、読み手では上位 4 バイトが前の値の
まま残った (リトルエンディアンで変数が 0 のときだけ正しい)。読む前に上位
ビットに値を入れておき、正しく置き換わることを確かめる。
*/

#include "ofd.h"
#include "ofd_prototype.h"
#include "ofd_test.h"

#include <stdlib.h>
#include <string.h>

#define TMPOUT "ofd_test_out.bin"

void readout(FILE *);   /* post/readout.c */

/* 値 v を持つ配列を確保 */
static double *dfill(size_t n, double v0)
{
	double *p = (double *)malloc((n + 1) * sizeof(double));
	for (size_t i = 0; i < n; i++) p[i] = v0 + 0.5 * (double)i;
	return p;
}
static d_complex_t *cfill(size_t n, double v0)
{
	d_complex_t *p = (d_complex_t *)malloc((n + 1) * sizeof(d_complex_t));
	for (size_t i = 0; i < n; i++) { p[i].r = v0 + i; p[i].i = -v0 - 0.25 * (double)i; }
	return p;
}
static float *ffill(size_t n, float v0)
{
	float *p = (float *)malloc((n + 1) * sizeof(float));
	for (size_t i = 0; i < n; i++) p[i] = v0 + 0.125f * (float)(i % 1000);
	return p;
}

static int deq(const double *a, const double *b, size_t n) { return !memcmp(a, b, n * sizeof(double)); }
static int ceq(const d_complex_t *a, const d_complex_t *b, size_t n) { return !memcmp(a, b, n * sizeof(d_complex_t)); }
static int feq(const float *a, const float *b, size_t n) { return !memcmp(a, b, n * sizeof(float)); }

void test_outfile(void)
{
	/* ---- 書く側の値 ---- */
	strcpy(Title, "round trip");
	Nx = 3; Ny = 4; Nz = 2;
	Ni = 35; Nj = 5; Nk = 1; N0 = 41;
	NN = 180;
	NFreq1 = 3; NFreq2 = 2;
	NFeed = 2; NPoint = 3;
	Niter = 7; Ntime = 123;
	Solver.maxiter = 20; Solver.nout = 5;
	Dt = 1.25e-12;
	NGline = 4;
	IPlanewave = 1;
	memset(&Planewave, 0, sizeof(Planewave));
	Planewave.theta = 30; Planewave.phi = 45; Planewave.pol = 2;
	Planewave.ei[0] = 0.1; Planewave.hi[2] = -0.7; Planewave.ri[1] = 0.6; Planewave.r0[2] = 3.5; Planewave.ai = 7e9;

	Xn = dfill(Nx + 1, 1); Yn = dfill(Ny + 1, 2); Zn = dfill(Nz + 1, 3);
	Xc = dfill(Nx, 4); Yc = dfill(Ny, 5); Zc = dfill(Nz, 6);
	Eiter = dfill(Niter, 7); Hiter = dfill(Niter, 8);
	const size_t nt = (size_t)(Solver.maxiter + 1);
	VFeed = dfill(NFeed * nt, 9); IFeed = dfill(NFeed * nt, 10); VPoint = dfill(NPoint * nt, 11);
	Freq1 = dfill(NFreq1, 1e9); Freq2 = dfill(NFreq2, 2e9);
	Feed = (feed_t *)calloc(NFeed, sizeof(feed_t));
	Feed[0].dir = 'Z'; Feed[0].i = 1; Feed[0].k = 2; Feed[0].volt = 1; Feed[0].z0 = 50;
	Feed[1].dir = 'X'; Feed[1].j = 3; Feed[1].delay = 1e-10; Feed[1].z0 = 75;
	Zin = cfill(NFeed * NFreq1, 12); Ref = dfill(NFeed * NFreq1, 13);
	Pin[0] = dfill(NFeed * NFreq2, 14); Pin[1] = dfill(NFeed * NFreq2, 15);
	Spara = cfill(NPoint * NFreq1, 16);
	Gline = (double (*)[2][3])dfill(NGline * 6, 17);
	NSurface = 5;
	Surface = (surface_t *)calloc(NSurface, sizeof(surface_t));
	for (int n = 0; n < NSurface; n++) { Surface[n].x = n; Surface[n].ds = 0.1 * n; }
	d_complex_t ***S[6] = {&SurfaceEx, &SurfaceEy, &SurfaceEz, &SurfaceHx, &SurfaceHy, &SurfaceHz};
	for (int c = 0; c < 6; c++) {
		*S[c] = (d_complex_t **)malloc(NFreq2 * sizeof(d_complex_t *));
		for (int f = 0; f < NFreq2; f++) (*S[c])[f] = cfill(NSurface, 100 * c + 10 * f);
	}
	float **cf[12] = {&cEx_r, &cEx_i, &cEy_r, &cEy_i, &cEz_r, &cEz_i, &cHx_r, &cHx_i, &cHy_r, &cHy_i, &cHz_r, &cHz_i};
	for (int c = 0; c < 12; c++) *cf[c] = ffill((size_t)NN * NFreq2, (float)(1000 * c));

	/* ---- 書く ---- */
	FILE *fp = fopen(TMPOUT, "wb");
	CHECK(fp != NULL);
	if (fp == NULL) return;
	writeout(fp);
	fclose(fp);

	/* ---- 控えて消す ---- */
	char title0[256];
	memcpy(title0, Title, sizeof(title0));
	const int sc[13] = {Nx, Ny, Nz, Ni, Nj, Nk, N0, NFreq1, NFreq2, NFeed, NPoint, Niter, Ntime};
	const int64_t nn0 = NN, ngl0 = NGline, nsf0 = NSurface;
	const int mi0 = Solver.maxiter, no0 = Solver.nout, ipw0 = IPlanewave;
	const double dt0 = Dt;
	planewave_t pw0 = Planewave;
	double *keep_d[17] = {Xn, Yn, Zn, Xc, Yc, Zc, Eiter, Hiter, VFeed, IFeed, VPoint, Freq1, Freq2, Ref, Pin[0], Pin[1], (double *)Gline};
	feed_t *feed0 = Feed;
	d_complex_t *zin0 = Zin, *sp0 = Spara;
	surface_t *surf0 = Surface;
	d_complex_t **s0[6];
	for (int c = 0; c < 6; c++) s0[c] = *S[c];
	float *cf0[12];
	for (int c = 0; c < 12; c++) cf0[c] = *cf[c];

	memset(Title, 0, sizeof(Title));
	Nx = Ny = Nz = Ni = Nj = Nk = N0 = 0;
	NFreq1 = NFreq2 = NFeed = NPoint = Niter = Ntime = 0;
	NN = 0;
	Solver.maxiter = Solver.nout = 0;
	Dt = 0;
	IPlanewave = 0;
	memset(&Planewave, 0, sizeof(Planewave));
	/* 上位 4 バイトにごみ : int64_t へ 4 バイトだけ読むと残る */
	NGline   = (int64_t)0x1234 << 32;
	NSurface = (int64_t)0x5678 << 32;

	/* ---- 読む (ofd_post の readout) ---- */
	fp = fopen(TMPOUT, "rb");
	CHECK(fp != NULL);
	if (fp == NULL) return;
	readout(fp);
	/* ファイルの最後まで読み切っていること */
	CHECK(fgetc(fp) == EOF);
	fclose(fp);
	remove(TMPOUT);

	/* ---- 比べる ---- */
	CHECK(!memcmp(Title, title0, sizeof(title0)));
	const int sr[13] = {Nx, Ny, Nz, Ni, Nj, Nk, N0, NFreq1, NFreq2, NFeed, NPoint, Niter, Ntime};
	CHECK(!memcmp(sc, sr, sizeof(sc)));
	CHECK(NN == nn0);
	CHECK(NGline == ngl0);
	CHECK(NSurface == nsf0);
	CHECK(Solver.maxiter == mi0 && Solver.nout == no0 && IPlanewave == ipw0);
	CHECK(Dt == dt0);
	CHECK(!memcmp(&Planewave, &pw0, sizeof(pw0)));

	const size_t len_d[17] = {(size_t)Nx + 1, (size_t)Ny + 1, (size_t)Nz + 1, (size_t)Nx, (size_t)Ny, (size_t)Nz,
		(size_t)Niter, (size_t)Niter, NFeed * nt, NFeed * nt, NPoint * nt, (size_t)NFreq1, (size_t)NFreq2,
		(size_t)NFeed * NFreq1, (size_t)NFeed * NFreq2, (size_t)NFeed * NFreq2, (size_t)NGline * 6};
	double *got_d[17] = {Xn, Yn, Zn, Xc, Yc, Zc, Eiter, Hiter, VFeed, IFeed, VPoint, Freq1, Freq2, Ref, Pin[0], Pin[1], (double *)Gline};
	int okd = 1;
	for (int a = 0; a < 17; a++) okd &= deq(got_d[a], keep_d[a], len_d[a]);
	CHECK(okd);
	CHECK(!memcmp(Feed, feed0, NFeed * sizeof(feed_t)));
	CHECK(ceq(Zin, zin0, (size_t)NFeed * NFreq1));
	CHECK(ceq(Spara, sp0, (size_t)NPoint * NFreq1));
	CHECK(!memcmp(Surface, surf0, NSurface * sizeof(surface_t)));
	int oks = 1;
	for (int c = 0; c < 6; c++) for (int f = 0; f < NFreq2; f++) oks &= ceq((*S[c])[f], s0[c][f], (size_t)NSurface);
	CHECK(oks);
	int okf = 1;
	for (int c = 0; c < 12; c++) okf &= feq(*cf[c], cf0[c], (size_t)NN * NFreq2);
	CHECK(okf);

	/* ---- 後片付け (読んだ側と控えの両方) ---- */
	for (int a = 0; a < 17; a++) { free(got_d[a]); free(keep_d[a]); }
	free(Feed); free(feed0); free(Zin); free(zin0); free(Spara); free(sp0);
	free(Surface); free(surf0);
	for (int c = 0; c < 6; c++) {
		for (int f = 0; f < NFreq2; f++) { free((*S[c])[f]); free(s0[c][f]); }
		free(*S[c]); free(s0[c]);
		*S[c] = NULL;
	}
	for (int c = 0; c < 12; c++) { free(*cf[c]); free(cf0[c]); *cf[c] = NULL; }
	Xn = Yn = Zn = Xc = Yc = Zc = NULL;
	Eiter = Hiter = VFeed = IFeed = VPoint = Freq1 = Freq2 = Ref = NULL;
	Pin[0] = Pin[1] = NULL;
	Gline = NULL;
	Feed = NULL; Zin = Spara = NULL; Surface = NULL;
	NFeed = NPoint = NFreq1 = NFreq2 = 0;
	NGline = NSurface = 0;
	IPlanewave = 0;
	memset(&Planewave, 0, sizeof(Planewave));
	Dt = 0;
}
