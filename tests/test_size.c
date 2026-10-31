/*
test_size.c

領域分割と配列レイアウト setupSize() (sol/setupSize.c) のテスト。

期待値の出所 (分割と索引の定義):
  - 全ランクの担当範囲 [iMin, iMax) を集めると [0, Nx) をちょうど覆う
    (隙間・重なりが無い)。各ランクは空でない。最後のランクが端数を持つ
  - ランク番号は (Ipx, Ipy, Ipz) に Ipz が最も速く回る順で対応する
  - 配列は担当範囲に吸収境界の余白 (Mur 1 層 / PML L 層) を足した直方体で、
    NA() は余白の最小隅で 0、最大隅で NN - 1、その間で 1 対 1
  - 分割数がセル数を超えると 1 を返す (以前は表示だけして続行し、担当範囲が
    空や領域外になって結果が黙って変わっていた)
*/

#include "ofd.h"
#include "ofd_prototype.h"
#include "ofd_test.h"

#include <stdlib.h>
#include <string.h>

static void check_partition(int nx, int ny, int nz, int px, int py, int pz, int abc, int l)
{
	Nx = nx; Ny = ny; Nz = nz;
	/* setupSize() は分割数を引数でも受け取るが、ブロック幅の計算は
	   グローバルの Npx/Npy/Npz で行う (実行時は同じ値)。両方を揃える */
	Npx = px; Npy = py; Npz = pz;
	iABC = abc;
	cPML.l = l;
	const int np = px * py * pz;
	int *covx = (int *)calloc(nx, sizeof(int));
	int *covy = (int *)calloc(ny, sizeof(int));
	int *covz = (int *)calloc(nz, sizeof(int));
	int ok_rank = 1, ok_layout = 1, ok_nonempty = 1;
	const int lm = (abc == 0) ? 1 : l;

	for (int rank = 0; rank < np; rank++) {
		CHECK(setupSize(px, py, pz, rank) == 0);
		/* ランク番号 = Ipx * py * pz + Ipy * pz + Ipz */
		ok_rank &= (rank == (Ipx * py * pz) + (Ipy * pz) + Ipz);
		ok_nonempty &= (iMin < iMax) && (jMin < jMax) && (kMin < kMax);
		/* 1 方向に 1 回だけ数える (他の 2 方向の添字が 0 のランクで) */
		if ((Ipy == 0) && (Ipz == 0)) for (int i = iMin; i < iMax; i++) covx[i]++;
		if ((Ipx == 0) && (Ipz == 0)) for (int j = jMin; j < jMax; j++) covy[j]++;
		if ((Ipx == 0) && (Ipy == 0)) for (int k = kMin; k < kMax; k++) covz[k]++;
		/* 最後のブロックは端まで */
		if (Ipx == px - 1) ok_rank &= (iMax == nx);
		if (Ipy == py - 1) ok_rank &= (jMax == ny);
		if (Ipz == pz - 1) ok_rank &= (kMax == nz);

		/* レイアウト : 余白の隅で 0 / NN-1、全点が相異なる */
		ok_layout &= (NA(iMin - lm, jMin - lm, kMin - lm) == 0);
		ok_layout &= (NA(iMax + lm, jMax + lm, kMax + lm) == NN - 1);
		const int64_t count = (int64_t)(iMax - iMin + 2 * lm + 1) * (jMax - jMin + 2 * lm + 1) * (kMax - kMin + 2 * lm + 1);
		ok_layout &= (count == NN);
		char *seen = (char *)calloc(NN, 1);
		for (int i = iMin - lm; i <= iMax + lm; i++) {
		for (int j = jMin - lm; j <= jMax + lm; j++) {
		for (int k = kMin - lm; k <= kMax + lm; k++) {
			const int64_t q = NA(i, j, k);
			if ((q < 0) || (q >= NN) || seen[q]) { ok_layout = 0; continue; }
			seen[q] = 1;
		}
		}
		}
		free(seen);
	}
	int ok_cover = 1;
	for (int i = 0; i < nx; i++) ok_cover &= (covx[i] == 1);
	for (int j = 0; j < ny; j++) ok_cover &= (covy[j] == 1);
	for (int k = 0; k < nz; k++) ok_cover &= (covz[k] == 1);
	CHECK(ok_rank);
	CHECK(ok_nonempty);
	CHECK(ok_cover);
	CHECK(ok_layout);
	if (!(ok_rank && ok_nonempty && ok_cover && ok_layout)) {
		printf("    %dx%dx%d / %dx%dx%d abc=%d : rank=%d nonempty=%d cover=%d layout=%d\n",
			nx, ny, nz, px, py, pz, abc, ok_rank, ok_nonempty, ok_cover, ok_layout);
	}
	free(covx); free(covy); free(covz);
}

void test_size(void)
{
	const int abc0 = iABC;

	check_partition(7, 5, 9, 1, 1, 1, 0, 0);
	check_partition(7, 5, 9, 2, 1, 3, 0, 0);    /* 割り切れない (端数は最後) */
	check_partition(7, 5, 9, 3, 2, 2, 1, 4);    /* PML の余白 */
	check_partition(4, 4, 4, 4, 1, 1, 0, 0);    /* ちょうどセル数 */
	check_partition(10, 3, 6, 1, 3, 6, 1, 2);

	/* 分割数がセル数を超える : 1 を返す (rank 0 だけが表示する) */
	Nx = 10; Ny = 10; Nz = 2;
	Npx = Npy = Npz = 1;
	iABC = 0;
	CHECK(setupSize(1, 1, 3, 1) == 1);
	CHECK(setupSize(1, 1, 4, 1) == 1);
	CHECK(setupSize(11, 1, 1, 1) == 1);
	CHECK(setupSize(1, 1, 2, 0) == 0);

	iABC = abc0;
	Nx = Ny = Nz = 0;
	Npx = Npy = Npz = 1;
}
