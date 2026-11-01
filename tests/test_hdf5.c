/*
test_hdf5.c

HDF5 出力 time_series_data.h5 (sol/outputHdf5.c) の値のテスト。

このファイルは表示側 (OpenFDTD-X) が読む契約で、構成 (どのパスに何が
あるか) は data/sample/hdf5_layout_check.sh が見ている。ここでは**値**、
特にソルバー内部の 1 次元配列 (余白込み、n = i*Ni + j*Nj + k*Nk + N0) から
(i,j,k,成分[,実部/虚部]) の自然な形状への詰め替えを確かめる。添字の順や
成分の取り違えは、構成の検査も "normal end" も通ってしまい、表示が
静かに転置・混色するだけになる。

期待値は include/ofd_hdf5.h のファイル構成の記述から作る:
  - 各成分・各節点に一意の整数値 (float で厳密) を入れ、余白には -1 を入れる。
    読み戻した配列の [..][i][j][k][c] がその値であり、-1 が 1 つも無いこと
  - /timeseries の time_H = time - Dt/2 (leapfrog)、itime はその順
  - hdf5 キーの interval で間引かれたステップは追記されない
  - /convergence/iter = n * nout
  - /metadata の Nx/Ny/Nz/Nsnapshot と Title
*/

#include "ofd.h"
#include "ofd_hdf5.h"
#include "ofd_test.h"
#include "hdf5.h"

#include <stdlib.h>
#include <string.h>

#define NXT (3)
#define NYT (2)
#define NZT (4)
#define NF  (2)

/* 一意な値 (float で厳密に表せる整数) */
static float val(int field, int comp, int i, int j, int k)
{
	return (float)((field * 10000) + (comp * 1000) + (i * 100) + (j * 10) + k + 1);
}

static float cval(int field, int ifreq, int comp, int reim, int i, int j, int k)
{
	return (float)((ifreq * 200000) + (reim * 100000) + val(field, comp, i, j, k));
}

/* データセットを丸ごと読む (要素数を返す。失敗時 -1) */
static int64_t read_all(hid_t file, const char *path, hid_t type, void *buf, size_t maxbytes, size_t elem)
{
	hid_t d = H5Dopen(file, path, H5P_DEFAULT);
	if (d < 0) return -1;
	hid_t s = H5Dget_space(d);
	const hssize_t n = H5Sget_simple_extent_npoints(s);
	int64_t ret = -1;
	if ((n > 0) && ((size_t)n * elem <= maxbytes) && (H5Dread(d, type, H5S_ALL, H5S_ALL, H5P_DEFAULT, buf) >= 0)) {
		ret = (int64_t)n;
	}
	H5Sclose(s);
	H5Dclose(d);
	return ret;
}

void test_hdf5(void)
{
	char title0[256];
	memcpy(title0, Title, sizeof(title0));
	strcpy(Title, "hdf5 test");

	/* 格子と配列レイアウト (余白 1 層。setupSize() と同じ式) */
	Nx = NXT; Ny = NYT; Nz = NZT;
	iMin = jMin = kMin = 0;
	iMax = Nx; jMax = Ny; kMax = Nz;
	Nk = 1;
	Nj = (Nz + 3) * Nk;
	Ni = (Ny + 3) * Nj;
	N0 = Ni + Nj + Nk;
	NN = NA(Nx + 1, Ny + 1, Nz + 1) + 1;
	Xn = (double *)malloc((Nx + 1) * sizeof(double)); Xc = (double *)malloc(Nx * sizeof(double));
	Yn = (double *)malloc((Ny + 1) * sizeof(double)); Yc = (double *)malloc(Ny * sizeof(double));
	Zn = (double *)malloc((Nz + 1) * sizeof(double)); Zc = (double *)malloc(Nz * sizeof(double));
	for (int i = 0; i <= Nx; i++) Xn[i] = 0.01 * i;
	for (int j = 0; j <= Ny; j++) Yn[j] = 0.02 * j;
	for (int k = 0; k <= Nz; k++) Zn[k] = 0.03 * k;
	for (int i = 0; i < Nx; i++) Xc[i] = 0.01 * (i + 0.5);
	for (int j = 0; j < Ny; j++) Yc[j] = 0.02 * (j + 0.5);
	for (int k = 0; k < Nz; k++) Zc[k] = 0.03 * (k + 0.5);
	NGline = 0; Gline = NULL; MGline = NULL;
	NFreq1 = 0; NFeed = 0; NPoint = 0; IPlanewave = 0;
	NSurface = 0; Surface = NULL;
	Dt = 1e-12;
	Solver.maxiter = 100; Solver.nout = 5;
	Ntime = 0; Niter = 0;

	/* 場 : 物理領域に一意の値、余白に -1 */
	real_t *f[6];
	for (int c = 0; c < 6; c++) {
		f[c] = (real_t *)malloc((size_t)NN * sizeof(real_t));
		for (int64_t n = 0; n < NN; n++) f[c][n] = -1;
	}
	for (int i = 0; i <= Nx; i++) {
	for (int j = 0; j <= Ny; j++) {
	for (int k = 0; k <= Nz; k++) {
		for (int c = 0; c < 6; c++) f[c][NA(i, j, k)] = val(c / 3, c % 3, i, j, k);
	}
	}
	}

	/* 近傍界 DFT : {NFreq2, NN} の 12 本 */
	NFreq2 = NF;
	Freq2 = (double *)malloc(NF * sizeof(double));
	Freq2[0] = 1e9; Freq2[1] = 2e9;
	float *cf[12];
	for (int b = 0; b < 12; b++) {
		cf[b] = (float *)malloc((size_t)NF * NN * sizeof(float));
		for (int64_t n = 0; n < NF * NN; n++) cf[b][n] = -1;
	}
	for (int ifreq = 0; ifreq < NF; ifreq++) {
	for (int i = 0; i <= Nx; i++) {
	for (int j = 0; j <= Ny; j++) {
	for (int k = 0; k <= Nz; k++) {
		const int64_t n = (int64_t)ifreq * NN + NA(i, j, k);
		/* b = field*6 + comp*2 + reim (cEx_r, cEx_i, cEy_r, ... の順) */
		for (int b = 0; b < 12; b++) cf[b][n] = cval(b / 6, ifreq, (b / 2) % 3, b % 2, i, j, k);
	}
	}
	}
	}

	/* 損失 {NFreq2, NN} (double) */
	double *ploss = (double *)malloc((size_t)NF * NN * sizeof(double));
	for (int64_t n = 0; n < NF * NN; n++) ploss[n] = -1;
	for (int ifreq = 0; ifreq < NF; ifreq++) {
	for (int i = 0; i <= Nx; i++) {
	for (int j = 0; j <= Ny; j++) {
	for (int k = 0; k <= Nz; k++) {
		ploss[(int64_t)ifreq * NN + NA(i, j, k)] = cval(0, ifreq, 0, 0, i, j, k) + 0.5;
	}
	}
	}
	}

	/* ---- 書く ---- */
	const int out0 = Hdf5Output, int0 = Hdf5Interval;
	Hdf5Output = 1;
	Hdf5Interval = 10;
	remove(OFD_H5_FILE);
	CHECK(hdf5_open(1) == 1);
	/* itime 10, 20 は書かれ、15 は間引かれる。2 枚目の E は値を 1 ずらす */
	hdf5_write_snapshot(10, 11 * Dt, f[0], f[1], f[2], f[3], f[4], f[5]);
	hdf5_write_snapshot(15, 16 * Dt, f[0], f[1], f[2], f[3], f[4], f[5]);
	for (int64_t n = 0; n < NN; n++) if (f[0][n] >= 0) f[0][n] += 0.5f;
	hdf5_write_snapshot(20, 21 * Dt, f[0], f[1], f[2], f[3], f[4], f[5]);
	for (int64_t n = 0; n < NN; n++) if (f[0][n] >= 0) f[0][n] -= 0.5f;
	hdf5_write_freqdomain(cf[0], cf[1], cf[2], cf[3], cf[4], cf[5],
	                      cf[6], cf[7], cf[8], cf[9], cf[10], cf[11]);
	hdf5_write_loss(ploss);
	const double eit[3] = {1.0, 0.5, 0.25}, hit[3] = {2.0, 1.0, 0.125};
	hdf5_write_convergence(3, eit, hit);
	hdf5_close();

	/* ---- 読む ---- */
	const int64_t nnode = (int64_t)(Nx + 1) * (Ny + 1) * (Nz + 1);
	const size_t big = (size_t)NF * nnode * 6 * 2 * sizeof(float) + 1024;
	float *buf = (float *)malloc(big);
	hid_t file = H5Fopen(OFD_H5_FILE, H5F_ACC_RDONLY, H5P_DEFAULT);
	CHECK(file >= 0);
	if ((file >= 0) && (buf != NULL)) {
		/* /timeseries/E, H : {nsnap, Nx+1, Ny+1, Nz+1, 3} */
		int ok_e = 1, ok_h = 1;
		CHECK(read_all(file, "/timeseries/E", H5T_NATIVE_FLOAT, buf, big, sizeof(float)) == 2 * nnode * 3);
		for (int s = 0; s < 2; s++) {
			int64_t m = (int64_t)s * nnode * 3;
			for (int i = 0; i <= Nx; i++) {
			for (int j = 0; j <= Ny; j++) {
			for (int k = 0; k <= Nz; k++) {
				for (int c = 0; c < 3; c++) {
					const float e = val(0, c, i, j, k) + (((s == 1) && (c == 0)) ? 0.5f : 0.0f);
					ok_e &= (buf[m++] == e);
				}
			}
			}
			}
		}
		CHECK(ok_e);
		CHECK(read_all(file, "/timeseries/H", H5T_NATIVE_FLOAT, buf, big, sizeof(float)) == 2 * nnode * 3);
		{
			int64_t m = nnode * 3;   /* 2 枚目 */
			for (int i = 0; i <= Nx; i++) {
			for (int j = 0; j <= Ny; j++) {
			for (int k = 0; k <= Nz; k++) {
				for (int c = 0; c < 3; c++) ok_h &= (buf[m++] == val(1, c, i, j, k));
			}
			}
			}
		}
		CHECK(ok_h);

		/* 時刻 : E は与えた t、H は半ステップ前。itime は間引き後の順 */
		double tv[8];
		int iv[8];
		CHECK(read_all(file, "/timeseries/itime", H5T_NATIVE_INT, iv, sizeof(iv), sizeof(int)) == 2);
		CHECK((iv[0] == 10) && (iv[1] == 20));
		CHECK(read_all(file, "/timeseries/time", H5T_NATIVE_DOUBLE, tv, sizeof(tv), sizeof(double)) == 2);
		CHECK_NEAR(tv[1], 21 * Dt, 1e-25);
		CHECK(read_all(file, "/timeseries/time_H", H5T_NATIVE_DOUBLE, tv, sizeof(tv), sizeof(double)) == 2);
		CHECK_NEAR(tv[0], 10.5 * Dt, 1e-25);
		CHECK_NEAR(tv[1], 20.5 * Dt, 1e-25);

		/* /freqdomain/E, H : {NFreq2, Nx+1, Ny+1, Nz+1, 3, 2} */
		for (int field = 0; field < 2; field++) {
			const char *path = (field == 0) ? "/freqdomain/E" : "/freqdomain/H";
			int ok = (read_all(file, path, H5T_NATIVE_FLOAT, buf, big, sizeof(float)) == NF * nnode * 6);
			int64_t m = 0;
			for (int ifreq = 0; ok && (ifreq < NF); ifreq++) {
			for (int i = 0; i <= Nx; i++) {
			for (int j = 0; j <= Ny; j++) {
			for (int k = 0; k <= Nz; k++) {
				for (int c = 0; c < 3; c++) {
				for (int r = 0; r < 2; r++) {
					ok &= (buf[m++] == cval(field, ifreq, c, r, i, j, k));
				}
				}
			}
			}
			}
			}
			CHECK(ok);
		}
		CHECK(read_all(file, "/freqdomain/freq", H5T_NATIVE_DOUBLE, tv, sizeof(tv), sizeof(double)) == NF);
		CHECK((tv[0] == 1e9) && (tv[1] == 2e9));

		/* /loss/P_loss : {NFreq2, Nx+1, Ny+1, Nz+1} */
		{
			int ok = (read_all(file, "/loss/P_loss", H5T_NATIVE_FLOAT, buf, big, sizeof(float)) == NF * nnode);
			int64_t m = 0;
			for (int ifreq = 0; ok && (ifreq < NF); ifreq++) {
			for (int i = 0; i <= Nx; i++) {
			for (int j = 0; j <= Ny; j++) {
			for (int k = 0; k <= Nz; k++) {
				ok &= (buf[m++] == cval(0, ifreq, 0, 0, i, j, k) + 0.5f);
			}
			}
			}
			}
			CHECK(ok);
		}

		/* /convergence : iter = n * nout */
		CHECK(read_all(file, "/convergence/iter", H5T_NATIVE_INT, iv, sizeof(iv), sizeof(int)) == 3);
		CHECK((iv[0] == 0) && (iv[1] == 5) && (iv[2] == 10));
		CHECK(read_all(file, "/convergence/H", H5T_NATIVE_DOUBLE, tv, sizeof(tv), sizeof(double)) == 3);
		CHECK((tv[0] == 2.0) && (tv[2] == 0.125));

		/* /metadata */
		CHECK(read_all(file, "/metadata/Ny", H5T_NATIVE_INT, iv, sizeof(iv), sizeof(int)) == 1);
		CHECK(iv[0] == NYT);
		CHECK(read_all(file, "/metadata/Nsnapshot", H5T_NATIVE_INT, iv, sizeof(iv), sizeof(int)) == 1);
		CHECK(iv[0] == 2);
		CHECK(read_all(file, "/geometry/Zn", H5T_NATIVE_DOUBLE, tv, sizeof(tv), sizeof(double)) == NZT + 1);
		CHECK_NEAR(tv[NZT], 0.03 * NZT, 1e-15);
		{
			char t[sizeof(Title)];
			memset(t, 0, sizeof(t));
			hid_t d = H5Dopen(file, "/metadata/Title", H5P_DEFAULT);
			CHECK(d >= 0);
			if (d >= 0) {
				hid_t ty = H5Dget_type(d);
				H5Dread(d, ty, H5S_ALL, H5S_ALL, H5P_DEFAULT, t);
				H5Tclose(ty);
				H5Dclose(d);
			}
			CHECK(strcmp(t, "hdf5 test") == 0);
		}
		H5Fclose(file);
	}
	remove(OFD_H5_FILE);

	/* hdf5 = 0 なら何も作らない */
	Hdf5Output = 0;
	CHECK(hdf5_open(1) == 0);
	{
		FILE *fp = fopen(OFD_H5_FILE, "rb");
		CHECK(fp == NULL);
		if (fp != NULL) fclose(fp);
	}

	Hdf5Output = out0;
	Hdf5Interval = int0;
	free(buf);
	free(ploss);
	for (int c = 0; c < 6; c++) free(f[c]);
	for (int b = 0; b < 12; b++) free(cf[b]);
	free(Freq2); Freq2 = NULL;
	free(Xn); free(Yn); free(Zn); free(Xc); free(Yc); free(Zc);
	Xn = Yn = Zn = Xc = Yc = Zc = NULL;
	NFreq2 = 0;
	Nx = Ny = Nz = 0; NN = 0;
	memcpy(Title, title0, sizeof(title0));
}
