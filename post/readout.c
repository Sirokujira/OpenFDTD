/*
readout.c

read ofd.out
*/

#include "ofd.h"

// 個数 0 の配列は malloc(0) が NULL を返し得る。fread(NULL, size, 0, fp) は
// 読む量が 0 でも C の規格上は未定義動作なので、0 個なら呼ばない
// (sol/writeout.c の wr() と対)
static size_t rd(void *p, size_t size, size_t n, FILE *fp)
{
	return (n > 0) ? fread(p, size, n, fp) : 0;
}

void readout(FILE *fp)
{
	size_t num = 0;
	size_t size;
	// NGline / NSurface は ofd.out では 4 バイト (int)。int64_t の変数へ
	// sizeof(int) で直接読むと上位 4 バイトが前の値のまま残るので、int で
	// 受けてから代入する (sol/writeout.c と対)
	int ngline = 0, nsurface = 0;

	num += rd(Title,           sizeof(char),        256, fp);
	num += rd(&Nx,             sizeof(int),         1,   fp);
	num += rd(&Ny,             sizeof(int),         1,   fp);
	num += rd(&Nz,             sizeof(int),         1,   fp);
	num += rd(&Ni,             sizeof(int),         1,                             fp);
	num += rd(&Nj,             sizeof(int),         1,                             fp);
	num += rd(&Nk,             sizeof(int),         1,                             fp);
	num += rd(&N0,             sizeof(int),         1,                             fp);
	num += rd(&NN,             sizeof(int64_t),     1,   fp);
	num += rd(&NFreq1,         sizeof(int),         1,   fp);
	num += rd(&NFreq2,         sizeof(int),         1,   fp);
	num += rd(&NFeed,          sizeof(int),         1,   fp);
	num += rd(&NPoint,         sizeof(int),         1,   fp);
	num += rd(&Niter,          sizeof(int),         1,   fp);
	num += rd(&Ntime,          sizeof(int),         1,   fp);
	num += rd(&Solver.maxiter, sizeof(int),         1,   fp);
	num += rd(&Solver.nout,    sizeof(int),         1,   fp);
	num += rd(&Dt,             sizeof(double),      1,   fp);
	num += rd(&ngline,         sizeof(int),         1,   fp);
	NGline = ngline;
	num += rd(&IPlanewave,     sizeof(int),         1,   fp);
	num += rd(&Planewave,      sizeof(planewave_t), 1,   fp);

	Xn     =         (double *)malloc(sizeof(double)      * (Nx + 1));
	Yn     =         (double *)malloc(sizeof(double)      * (Ny + 1));
	Zn     =         (double *)malloc(sizeof(double)      * (Nz + 1));
	Xc     =         (double *)malloc(sizeof(double)      * Nx);
	Yc     =         (double *)malloc(sizeof(double)      * Ny);
	Zc     =         (double *)malloc(sizeof(double)      * Nz);
	Eiter  =         (double *)malloc(sizeof(double)      * Niter);
	Hiter  =         (double *)malloc(sizeof(double)      * Niter);
	VFeed  =         (double *)malloc(sizeof(double)      * NFeed  * (Solver.maxiter + 1));
	IFeed  =         (double *)malloc(sizeof(double)      * NFeed  * (Solver.maxiter + 1));
	VPoint =         (double *)malloc(sizeof(double)      * NPoint * (Solver.maxiter + 1));
	Freq1  =         (double *)malloc(sizeof(double)      * NFreq1);
	Freq2  =         (double *)malloc(sizeof(double)      * NFreq2);
	Feed   =         (feed_t *)malloc(sizeof(feed_t)      * NFeed);
	Zin    =    (d_complex_t *)malloc(sizeof(d_complex_t) * NFeed * NFreq1);
	Ref    =         (double *)malloc(sizeof(double)      * NFeed * NFreq1);
	Pin[0] =         (double *)malloc(sizeof(double)      * NFeed * NFreq2);
	Pin[1] =         (double *)malloc(sizeof(double)      * NFeed * NFreq2);
	Spara  =    (d_complex_t *)malloc(sizeof(d_complex_t) * NPoint * NFreq1);
	Gline  = (double (*)[2][3])malloc(sizeof(double)      * NGline * 2 * 3);

	num += rd(Xn,     sizeof(double),      Nx + 1,                        fp);
	num += rd(Yn,     sizeof(double),      Ny + 1,                        fp);
	num += rd(Zn,     sizeof(double),      Nz + 1,                        fp);
	num += rd(Xc,     sizeof(double),      Nx,                            fp);
	num += rd(Yc,     sizeof(double),      Ny,                            fp);
	num += rd(Zc,     sizeof(double),      Nz,                            fp);
	num += rd(Eiter,  sizeof(double),      Niter,                         fp);
	num += rd(Hiter,  sizeof(double),      Niter,                         fp);
	num += rd(VFeed,  sizeof(double),      NFeed  * (Solver.maxiter + 1), fp);
	num += rd(IFeed,  sizeof(double),      NFeed  * (Solver.maxiter + 1), fp);
	num += rd(VPoint, sizeof(double),      NPoint * (Solver.maxiter + 1), fp);
	num += rd(Freq1,  sizeof(double),      NFreq1,                        fp);
	num += rd(Freq2,  sizeof(double),      NFreq2,                        fp);
	num += rd(Feed,   sizeof(feed_t),      NFeed,                         fp);
	num += rd(Zin,    sizeof(d_complex_t), NFeed * NFreq1,                fp);
	num += rd(Ref,    sizeof(double),      NFeed * NFreq1,                fp);
	num += rd(Pin[0], sizeof(double),      NFeed * NFreq2,                fp);
	num += rd(Pin[1], sizeof(double),      NFeed * NFreq2,                fp);
	num += rd(Spara,  sizeof(d_complex_t), NPoint * NFreq1,               fp);
	num += rd(Gline,  sizeof(double),      NGline * 2 * 3,                fp);

	num += rd(&nsurface, sizeof(int),         1,        fp);
	NSurface = nsurface;

	size = NSurface * sizeof(surface_t);
	Surface = (surface_t *)malloc(size);
	num += rd(Surface,   sizeof(surface_t),   NSurface, fp);

	size = NFreq2 * sizeof(d_complex_t *);
	SurfaceEx = (d_complex_t **)malloc(size);
	SurfaceEy = (d_complex_t **)malloc(size);
	SurfaceEz = (d_complex_t **)malloc(size);
	SurfaceHx = (d_complex_t **)malloc(size);
	SurfaceHy = (d_complex_t **)malloc(size);
	SurfaceHz = (d_complex_t **)malloc(size);
	for (int ifreq = 0; ifreq < NFreq2; ifreq++) {
		size = NSurface * sizeof(d_complex_t);
		SurfaceEx[ifreq] = (d_complex_t *)malloc(size);
		SurfaceEy[ifreq] = (d_complex_t *)malloc(size);
		SurfaceEz[ifreq] = (d_complex_t *)malloc(size);
		SurfaceHx[ifreq] = (d_complex_t *)malloc(size);
		SurfaceHy[ifreq] = (d_complex_t *)malloc(size);
		SurfaceHz[ifreq] = (d_complex_t *)malloc(size);

		num += rd(SurfaceEx[ifreq], sizeof(d_complex_t), NSurface, fp);
		num += rd(SurfaceEy[ifreq], sizeof(d_complex_t), NSurface, fp);
		num += rd(SurfaceEz[ifreq], sizeof(d_complex_t), NSurface, fp);
		num += rd(SurfaceHx[ifreq], sizeof(d_complex_t), NSurface, fp);
		num += rd(SurfaceHy[ifreq], sizeof(d_complex_t), NSurface, fp);
		num += rd(SurfaceHz[ifreq], sizeof(d_complex_t), NSurface, fp);
	}

	size = NN * NFreq2 * sizeof(float);   // 要素は float (以前は sizeof(float *) で 2 倍確保していた)
	cEx_r = (float *)malloc(size);
	cEx_i = (float *)malloc(size);
	cEy_r = (float *)malloc(size);
	cEy_i = (float *)malloc(size);
	cEz_r = (float *)malloc(size);
	cEz_i = (float *)malloc(size);
	cHx_r = (float *)malloc(size);
	cHx_i = (float *)malloc(size);
	cHy_r = (float *)malloc(size);
	cHy_i = (float *)malloc(size);
	cHz_r = (float *)malloc(size);
	cHz_i = (float *)malloc(size);
	for (int ifreq = 0; ifreq < NFreq2; ifreq++) {
		int64_t n0 = ifreq * NN;
		num += rd(&cEx_r[n0], sizeof(float), NN, fp);
		num += rd(&cEx_i[n0], sizeof(float), NN, fp);
		num += rd(&cEy_r[n0], sizeof(float), NN, fp);
		num += rd(&cEy_i[n0], sizeof(float), NN, fp);
		num += rd(&cEz_r[n0], sizeof(float), NN, fp);
		num += rd(&cEz_i[n0], sizeof(float), NN, fp);
		num += rd(&cHx_r[n0], sizeof(float), NN, fp);
		num += rd(&cHx_i[n0], sizeof(float), NN, fp);
		num += rd(&cHy_r[n0], sizeof(float), NN, fp);
		num += rd(&cHy_i[n0], sizeof(float), NN, fp);
		num += rd(&cHz_r[n0], sizeof(float), NN, fp);
		num += rd(&cHz_i[n0], sizeof(float), NN, fp);
	}

	size_t num0;
	size = fread(&num0, sizeof(size_t), 1, fp);
	size = size;  // suppress gcc warning
	//printf("%zd %zd\n", num, num0);

	if (num != num0) {
		fprintf(stderr, "*** invalid file length : (%zd, %zd)\n", num0, num);
	}
}
