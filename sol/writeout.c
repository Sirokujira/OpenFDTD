/*
writeout.c

write ofd.out
*/

#include "ofd.h"

// 個数 0 の配列 (観測点が無いときの VPoint、形状が無いときの Gline など) は
// NULL のことがある。fwrite(NULL, size, 0, fp) は書く量が 0 でも C の規格上は
// 未定義動作 (引数は非 NULL と規定されている) なので、0 個なら呼ばない。
// 書く個数が 0 なのでファイルの中身は変わらない
static size_t wr(const void *p, size_t size, size_t n, FILE *fp)
{
	return (n > 0) ? fwrite(p, size, n, fp) : 0;
}

void writeout(FILE *fp)
{
	size_t num = 0;

	// NGline / NSurface は int64_t だが、ofd.out では従来どおり 4 バイト (int) で
	// 持つ (読み手 post/readout.c と OpenFDTD-X の互換のため)。以前は int64_t の
	// 変数を sizeof(int) で直接書いていて、リトルエンディアンでたまたま動いていた
	const int ngline   = (int)NGline;
	const int nsurface = (int)NSurface;

	num += wr(Title,           sizeof(char),        256,                           fp);
	num += wr(&Nx,             sizeof(int),         1,                             fp);
	num += wr(&Ny,             sizeof(int),         1,                             fp);
	num += wr(&Nz,             sizeof(int),         1,                             fp);
	num += wr(&Ni,             sizeof(int),         1,                             fp);
	num += wr(&Nj,             sizeof(int),         1,                             fp);
	num += wr(&Nk,             sizeof(int),         1,                             fp);
	num += wr(&N0,             sizeof(int),         1,                             fp);
	num += wr(&NN,             sizeof(int64_t),     1,                             fp);
	num += wr(&NFreq1,         sizeof(int),         1,                             fp);
	num += wr(&NFreq2,         sizeof(int),         1,                             fp);
	num += wr(&NFeed,          sizeof(int),         1,                             fp);
	num += wr(&NPoint,         sizeof(int),         1,                             fp);
	num += wr(&Niter,          sizeof(int),         1,                             fp);
	num += wr(&Ntime,          sizeof(int),         1,                             fp);
	num += wr(&Solver.maxiter, sizeof(int),         1,                             fp);
	num += wr(&Solver.nout,    sizeof(int),         1,                             fp);
	num += wr(&Dt,             sizeof(double),      1,                             fp);
	num += wr(&ngline,         sizeof(int),         1,                             fp);
	num += wr(&IPlanewave,     sizeof(int),         1,                             fp);
	num += wr(&Planewave,      sizeof(planewave_t), 1,                             fp);

	num += wr(Xn,              sizeof(double),      Nx + 1,                        fp);
	num += wr(Yn,              sizeof(double),      Ny + 1,                        fp);
	num += wr(Zn,              sizeof(double),      Nz + 1,                        fp);
	num += wr(Xc,              sizeof(double),      Nx,                            fp);
	num += wr(Yc,              sizeof(double),      Ny,                            fp);
	num += wr(Zc,              sizeof(double),      Nz,                            fp);
	num += wr(Eiter,           sizeof(double),      Niter,                         fp);
	num += wr(Hiter,           sizeof(double),      Niter,                         fp);
	num += wr(VFeed,           sizeof(double),      NFeed  * (Solver.maxiter + 1), fp);
	num += wr(IFeed,           sizeof(double),      NFeed  * (Solver.maxiter + 1), fp);
	num += wr(VPoint,          sizeof(double),      NPoint * (Solver.maxiter + 1), fp);
	num += wr(Freq1,           sizeof(double),      NFreq1,                        fp);
	num += wr(Freq2,           sizeof(double),      NFreq2,                        fp);
	num += wr(Feed,            sizeof(feed_t),      NFeed,                         fp);
	num += wr(Zin,             sizeof(d_complex_t), NFeed * NFreq1,                fp);
	num += wr(Ref,             sizeof(double),      NFeed * NFreq1,                fp);
	num += wr(Pin[0],          sizeof(double),      NFeed * NFreq2,                fp);
	num += wr(Pin[1],          sizeof(double),      NFeed * NFreq2,                fp);
	num += wr(Spara,           sizeof(d_complex_t), NPoint * NFreq1,               fp);
	num += wr(Gline,           sizeof(double),      NGline * 2 * 3,                fp);

	num += wr(&nsurface,       sizeof(int),         1,        fp);
	num += wr(Surface,         sizeof(surface_t),   NSurface, fp);
	for (int ifreq = 0; ifreq < NFreq2; ifreq++) {
		num += wr(SurfaceEx[ifreq], sizeof(d_complex_t), NSurface, fp);
		num += wr(SurfaceEy[ifreq], sizeof(d_complex_t), NSurface, fp);
		num += wr(SurfaceEz[ifreq], sizeof(d_complex_t), NSurface, fp);
		num += wr(SurfaceHx[ifreq], sizeof(d_complex_t), NSurface, fp);
		num += wr(SurfaceHy[ifreq], sizeof(d_complex_t), NSurface, fp);
		num += wr(SurfaceHz[ifreq], sizeof(d_complex_t), NSurface, fp);
	}

	for (int ifreq = 0; ifreq < NFreq2; ifreq++) {
		int64_t n0 = ifreq * NN;
		num += wr(&cEx_r[n0], sizeof(float), NN, fp);
		num += wr(&cEx_i[n0], sizeof(float), NN, fp);
		num += wr(&cEy_r[n0], sizeof(float), NN, fp);
		num += wr(&cEy_i[n0], sizeof(float), NN, fp);
		num += wr(&cEz_r[n0], sizeof(float), NN, fp);
		num += wr(&cEz_i[n0], sizeof(float), NN, fp);
		num += wr(&cHx_r[n0], sizeof(float), NN, fp);
		num += wr(&cHx_i[n0], sizeof(float), NN, fp);
		num += wr(&cHy_r[n0], sizeof(float), NN, fp);
		num += wr(&cHy_i[n0], sizeof(float), NN, fp);
		num += wr(&cHz_r[n0], sizeof(float), NN, fp);
		num += wr(&cHz_i[n0], sizeof(float), NN, fp);
	}

	fwrite(&num, sizeof(size_t), 1, fp);
}
