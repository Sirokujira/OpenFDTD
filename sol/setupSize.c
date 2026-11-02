/*
setupSize.c

setup array size
*/

#include "ofd.h"
#include "ofd_prototype.h"

static void idminmax(int n, int _np, int np, int ip, int *min, int *max)
{
	if (_np == 1) {
		*min = 0;
		*max = n;
	}
	else {
		// MPI
		const int nc = MAX(n / np, 1);
		*min = (ip + 0) * nc;
		*max = (ip + 1) * nc;
		if (ip == np - 1) {
			*max = n;
		}
	}
}


// 戻り値 : 0 = 正常 / 1 = 分割数がセル数を超える (呼び出し側で止めること)
int setupSize(int npx, int npy, int npz, int comm_rank)
{
	// too many prosess (MPI)
	// 以前は表示するだけで続行していたが、分割数がセル数を超えると担当範囲が
	// 空や領域外 (例: Nz = 2 を 4 分割すると rank 2 は [2,3)) になり、結果が
	// 黙って変わる (1 プロセスと一致しない) のに normal end で終わっていた。
	// 全ランクが同じ値で判定するので、呼び出し側は comm_check() で止める。
	if ((npx > Nx) || (npy > Ny) || (npz > Nz)) {
		if (comm_rank == 0) {
			fprintf(stderr, "*** too many processes = %dx%dx%d (limit = %dx%dx%d)\n", npx, npy, npz, Nx, Ny, Nz);
			fflush(stderr);
		}
		return 1;
	}

	// get block numbers : (Ipx, Ipy, Ipz) (MPI)
	int ip = 0;
	for (int i = 0; i < npx; i++) {
	for (int j = 0; j < npy; j++) {
	for (int k = 0; k < npz; k++) {
		if (comm_rank == ip) {
			Ipx = i;
			Ipy = j;
			Ipz = k;
		}
		ip++;
	}
	}
	}
	//printf("%d %d %d %d\n", comm_rank, Ipx, Ipy, Ipz); fflush(stdout);
	//assert(comm_rank == (Ipx * npy * npz) + (Ipy * npz) + Ipz);

	// min, max
	idminmax(Nx, npx, Npx, Ipx, &iMin, &iMax);
	idminmax(Ny, npy, Npy, Ipy, &jMin, &jMax);
	idminmax(Nz, npz, Npz, Ipz, &kMin, &kMax);
	//printf("%d %d %d %d %d %d %d\n", comm_rank, iMin, iMax, jMin, jMax, kMin, kMax);

	// array index

	const int lx = (iABC == 0) ? 1 : (iABC == 1) ? cPML.l : 0;
	const int ly = (iABC == 0) ? 1 : (iABC == 1) ? cPML.l : 0;
	const int lz = (iABC == 0) ? 1 : (iABC == 1) ? cPML.l : 0;

	Nk = 1;
	Nj = ((kMax - kMin) + (2 * lz) + 1);
	Ni = ((jMax - jMin) + (2 * ly) + 1) * Nj;
	N0 = -((iMin - lx) * Ni + (jMin - ly) * Nj + (kMin - lz) * Nk);
	NN = NA(iMax + lx, jMax + ly, kMax + lz) + 1;
	//printf("%d %d %zd %zd %zd %zd %zd\n", comm_size, comm_rank, Ni, Nj, Nk, N0, NN);
	assert(NA(iMin - lx, jMin - ly, kMin - lz) == 0);
	assert(NA(iMax + lx, jMax + ly, kMax + lz) == NN - 1);

	return 0;
}


// ABC (array size)
void setupABCsize(void)
{
	numMurHx = numMurHy = numMurHz = 0;
	numPmlEx = numPmlEy = numPmlEz = numPmlHx = numPmlHy = numPmlHz = 0;

	if      (iABC == 0) {
		setupMurHx(0);
		setupMurHy(0);
		setupMurHz(0);
		//printf("%zd %zd %zd\n", numMurHx, numMurHy, numMurHz);
	}
	else if (iABC == 1) {
		setupPmlEx(0);
		setupPmlEy(0);
		setupPmlEz(0);
		setupPmlHx(0);
		setupPmlHy(0);
		setupPmlHz(0);
		//printf("%zd %zd %zd %zd %zd %zd\n", numPmlEx, numPmlEy, numPmlEz, numPmlHx, numPmlHy, numPmlHz);
	}
}
