/*
comm_E.c (MPI)

領域境界の E ハロー交換

既存の comm_X/Y/Z は H しか交換していない。この FDTD の領域分割では
E 更新が読むのは H の近傍だけ (sol/updateEx.c 参照)、分散性材料や
インダクタ (sol/dispersionEx.c, sol/eload.c) はセル局所なので、
E のハローは必要なかった。

TPA (sol/updateTpa.c) は |E|^2 の colocated 近似のため隣接セルの
E 成分を読む初めての処理で、領域分割時にはハローが要る。必要なのは
各方向につき 1 成分・1 面だけ:

  updateTpaEx : Ey[j-1], Ez[k-1]
  updateTpaEy : Ex[i-1], Ez[k-1]
  updateTpaEz : Ex[i-1], Ey[j-1]

すなわち X 方向は Ex、Y 方向は Ey、Z 方向は Ez の -側ガード面。
comm_X/Y/Z と同じ Sendrecv 構造にして +側も同時に埋める (+側の値も
隣プロセスが正しく計算した値なので有害ではない)。

バッファはこのファイル内に閉じて確保する (ofd.h のグローバルを増やさない)。
*/

#ifdef _MPI
#include <mpi.h>
#endif

#include "ofd.h"

#ifdef _MPI

static real_t *sbuf = NULL, *rbuf = NULL;
static int64_t bufsize = 0;

// 必要なバッファを確保する (面の最大要素数)
static void alloc_buffer(int64_t num)
{
	if (num <= bufsize) return;
	free(sbuf);
	free(rbuf);
	sbuf = (real_t *)malloc(num * sizeof(real_t));
	rbuf = (real_t *)malloc(num * sizeof(real_t));
	bufsize = num;
}

/*
軸 idir (0=X, 1=Y, 2=Z) に垂直な面上の成分 e を交換する。
p1range/p2range は面内 2 方向の添字範囲 (X 面なら j, k)。
send/recv は交換する面の添字 (軸方向)。
*/
static void pack(int idir, int ip, const real_t e[], real_t buf[],
	const int p1[], const int p2[])
{
	const int64_t n2 = p2[1] - p2[0] + 1;
	for (int a = p1[0]; a <= p1[1]; a++) {
	for (int b = p2[0]; b <= p2[1]; b++) {
		const int64_t m = ((int64_t)(a - p1[0]) * n2) + (b - p2[0]);
		buf[m] = (idir == 0) ? e[NA(ip, a, b)]
		       : (idir == 1) ? e[NA(a, ip, b)]
		                     : e[NA(a, b, ip)];
	}
	}
}


static void unpack(int idir, int ip, const real_t buf[], real_t e[],
	const int p1[], const int p2[])
{
	const int64_t n2 = p2[1] - p2[0] + 1;
	for (int a = p1[0]; a <= p1[1]; a++) {
	for (int b = p2[0]; b <= p2[1]; b++) {
		const int64_t m = ((int64_t)(a - p1[0]) * n2) + (b - p2[0]);
		if      (idir == 0) e[NA(ip, a, b)] = buf[m];
		else if (idir == 1) e[NA(a, ip, b)] = buf[m];
		else                e[NA(a, b, ip)] = buf[m];
	}
	}
}


// idir 方向のハロー交換の共通処理
static void comm_E(int idir, real_t e[], int nproc, int iproc,
	int pmin, int pmax, const int p1[], const int p2[])
{
	if (nproc <= 1) return;

	MPI_Status status;
	const int tag = 0;
	const int64_t num = (int64_t)(p1[1] - p1[0] + 1) * (p2[1] - p2[0] + 1);
	alloc_buffer(num);
	if ((sbuf == NULL) || (rbuf == NULL)) return;

	// side 0 : -側の隣へ送り、-側ガード面 (pmin - 1) に受け取る
	// side 1 : +側の隣へ送り、+側ガード面 (pmax    ) に受け取る
	const int active[2] = {iproc > 0, iproc < nproc - 1};
	const int nbr[2]    = {iproc - 1, iproc + 1};
	const int isend[2]  = {pmin,      pmax - 1};
	const int irecv[2]  = {pmin - 1,  pmax};

	for (int side = 0; side < 2; side++) {
		if (!active[side]) continue;

		pack(idir, isend[side], e, sbuf, p1, p2);

		const int ipx = (idir == 0) ? nbr[side] : Ipx;
		const int ipy = (idir == 1) ? nbr[side] : Ipy;
		const int ipz = (idir == 2) ? nbr[side] : Ipz;
		const int dst = (ipx * Npy * Npz) + (ipy * Npz) + ipz;

		MPI_Sendrecv(sbuf, (int)num, MPI_REAL_T, dst, tag,
		             rbuf, (int)num, MPI_REAL_T, dst, tag, MPI_COMM_WORLD, &status);

		unpack(idir, irecv[side], rbuf, e, p1, p2);
	}
}

#endif  // _MPI


// X 面の Ex を交換する (updateTpaEy / updateTpaEz が Ex[i-1] を読む)
void comm_E_X(void)
{
#ifdef _MPI
	const int jr[2] = {jMin, jMax};
	const int kr[2] = {kMin, kMax};
	comm_E(0, Ex, Npx, Ipx, iMin, iMax, jr, kr);
#endif
}


// Y 面の Ey を交換する (updateTpaEx / updateTpaEz が Ey[j-1] を読む)
void comm_E_Y(void)
{
#ifdef _MPI
	const int ir[2] = {iMin, iMax};
	const int kr[2] = {kMin, kMax};
	comm_E(1, Ey, Npy, Ipy, jMin, jMax, ir, kr);
#endif
}


// Z 面の Ez を交換する (updateTpaEx / updateTpaEy が Ez[k-1] を読む)
void comm_E_Z(void)
{
#ifdef _MPI
	const int ir[2] = {iMin, iMax};
	const int jr[2] = {jMin, jMax};
	comm_E(2, Ez, Npz, Ipz, kMin, kMax, ir, jr);
#endif
}
