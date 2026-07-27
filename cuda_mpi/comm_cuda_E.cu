/*
comm_cuda_E.cu (CUDA + MPI)

領域境界の E ハロー交換 (GPU 版)

CPU+MPI 版 mpi/comm_E.c の GPU 対応版。交換する成分・面・向きは
まったく同じで、違いは E が device メモリ上にあり得ることだけ:
パック/アンパックをカーネルで行い、MPI には host バッファを渡す
(comm_cuda_X.cu の d2h / h2d と同じ構造)。

背景 (詳細は mpi/comm_E.c 冒頭):
  既存の comm_cuda_X/Y/Z は H 面しか交換していない。E 更新が読むのは
  H の近傍だけなので、これまで E のハローは不要だった。
  TPA (cuda/updateTpa.cu) は |E|^2 の colocated 近似のため隣接セルの
  E 成分を読む初めての処理で、領域分割時にはハローが要る。

  updateTpaEx : Ey[j-1], Ez[k-1]
  updateTpaEy : Ex[i-1], Ez[k-1]
  updateTpaEz : Ex[i-1], Ey[j-1]

すなわち X 方向は Ex、Y 方向は Ey、Z 方向は Ez の -側ガード面。
comm_X/Y/Z と同じ Sendrecv 構造にして +側も同時に埋める。

バッファはこのファイル内に閉じて確保する (ofd.h のグローバルを増やさない)。
面の最大要素数まで遅延確保し、以降は再利用する。
*/

#ifdef _MPI
#include <mpi.h>
#endif

#include "ofd.h"
#include "ofd_cuda.h"
#include "ofd_prototype.h"

#ifdef _MPI

// 面内座標 (a, b) と軸方向添字 ip から配列添字を作る
// idir = 0 (X 面) : i = ip, j = a,  k = b
// idir = 1 (Y 面) : i = a,  j = ip, k = b
// idir = 2 (Z 面) : i = a,  j = b,  k = ip
__host__ __device__ __inline__
static int64_t eaddr(int idir, int ip, int a, int b, const param_t *p)
{
	const int i = (idir == 0) ? ip : a;
	const int j = (idir == 0) ? a : ((idir == 1) ? ip : b);
	const int k = (idir == 2) ? ip : b;

	return (p->Ni * i) + (p->Nj * j) + (p->Nk * k) + p->N0;
}


// E to buffer (1 点)
__host__ __device__ __inline__
static void _pack(int idir, int ip, int a, int b, const real_t *e, real_t *buf,
	int a0, int b0, int b1, const param_t *p)
{
	const int64_t m = ((int64_t)(a - a0) * (b1 - b0 + 1)) + (b - b0);
	buf[m] = e[eaddr(idir, ip, a, b, p)];
}

// buffer to E (1 点)
__host__ __device__ __inline__
static void _unpack(int idir, int ip, int a, int b, real_t *e, const real_t *buf,
	int a0, int b0, int b1, const param_t *p)
{
	const int64_t m = ((int64_t)(a - a0) * (b1 - b0 + 1)) + (b - b0);
	e[eaddr(idir, ip, a, b, p)] = buf[m];
}


__global__
static void pack_gpu(int idir, int ip, const real_t *e, real_t *buf,
	int a0, int a1, int b0, int b1)
{
	const int a = a0 + threadIdx.x + (blockIdx.x * blockDim.x);
	const int b = b0 + threadIdx.y + (blockIdx.y * blockDim.y);
	if ((a <= a1) &&
	    (b <= b1)) {
		_pack(idir, ip, a, b, e, buf, a0, b0, b1, &d_Param);
	}
}


__global__
static void unpack_gpu(int idir, int ip, real_t *e, const real_t *buf,
	int a0, int a1, int b0, int b1)
{
	const int a = a0 + threadIdx.x + (blockIdx.x * blockDim.x);
	const int b = b0 + threadIdx.y + (blockIdx.y * blockDim.y);
	if ((a <= a1) &&
	    (b <= b1)) {
		_unpack(idir, ip, a, b, e, buf, a0, b0, b1, &d_Param);
	}
}


// バッファ (host は MPI 用、device はパック/アンパック用)
static real_t *sbuf = NULL, *rbuf = NULL;          // host
static real_t *d_sbuf = NULL, *d_rbuf = NULL;      // device (GPU 実行時のみ)
static int64_t bufnum = 0;

static void alloc_buffer(int64_t num)
{
	if (num <= bufnum) return;

	free(sbuf);
	free(rbuf);
	sbuf = (real_t *)malloc(num * sizeof(real_t));
	rbuf = (real_t *)malloc(num * sizeof(real_t));

	if (GPU) {
		if (d_sbuf != NULL) cuda_free(GPU, d_sbuf);
		if (d_rbuf != NULL) cuda_free(GPU, d_rbuf);
		cuda_malloc(GPU, UM, (void **)&d_sbuf, num * sizeof(real_t));
		cuda_malloc(GPU, UM, (void **)&d_rbuf, num * sizeof(real_t));
	}

	bufnum = num;
}


// E -> host バッファ
static void pack(int idir, int ip, const real_t *e, int a0, int a1, int b0, int b1, size_t size)
{
	if (GPU) {
		cudaMemcpyToSymbol(d_Param, &h_Param, sizeof(param_t));

		dim3 grid(CEIL(a1 - a0 + 1, bufBlock.x),
		          CEIL(b1 - b0 + 1, bufBlock.y));
		pack_gpu<<<grid, bufBlock>>>(idir, ip, e, d_sbuf, a0, a1, b0, b1);

		cuda_memcpy(GPU, sbuf, d_sbuf, size, cudaMemcpyDeviceToHost);

		if (UM) cudaDeviceSynchronize();
	}
	else {
		for (int a = a0; a <= a1; a++) {
		for (int b = b0; b <= b1; b++) {
			_pack(idir, ip, a, b, e, sbuf, a0, b0, b1, &h_Param);
		}
		}
	}
}


// host バッファ -> E
static void unpack(int idir, int ip, real_t *e, int a0, int a1, int b0, int b1, size_t size)
{
	if (GPU) {
		cudaMemcpyToSymbol(d_Param, &h_Param, sizeof(param_t));

		cuda_memcpy(GPU, d_rbuf, rbuf, size, cudaMemcpyHostToDevice);

		dim3 grid(CEIL(a1 - a0 + 1, bufBlock.x),
		          CEIL(b1 - b0 + 1, bufBlock.y));
		unpack_gpu<<<grid, bufBlock>>>(idir, ip, e, d_rbuf, a0, a1, b0, b1);

		if (UM) cudaDeviceSynchronize();
	}
	else {
		for (int a = a0; a <= a1; a++) {
		for (int b = b0; b <= b1; b++) {
			_unpack(idir, ip, a, b, e, rbuf, a0, b0, b1, &h_Param);
		}
		}
	}
}


// idir 方向のハロー交換 (mpi/comm_E.c の comm_E と同じ手順)
static void comm_cuda_E(int idir, real_t *e, int nproc, int iproc,
	int pmin, int pmax, int a0, int a1, int b0, int b1)
{
	if (nproc <= 1) return;

	MPI_Status status;
	const int tag = 0;
	const int64_t num = (int64_t)(a1 - a0 + 1) * (b1 - b0 + 1);
	alloc_buffer(num);
	if ((sbuf == NULL) || (rbuf == NULL)) return;
	const size_t size = num * sizeof(real_t);

	// side 0 : -側の隣へ送り、-側ガード面 (pmin - 1) に受け取る
	// side 1 : +側の隣へ送り、+側ガード面 (pmax    ) に受け取る
	const int active[2] = {iproc > 0, iproc < nproc - 1};
	const int nbr[2]    = {iproc - 1, iproc + 1};
	const int isend[2]  = {pmin,      pmax - 1};
	const int irecv[2]  = {pmin - 1,  pmax};

	for (int side = 0; side < 2; side++) {
		if (!active[side]) continue;

		pack(idir, isend[side], e, a0, a1, b0, b1, size);

		const int ipx = (idir == 0) ? nbr[side] : Ipx;
		const int ipy = (idir == 1) ? nbr[side] : Ipy;
		const int ipz = (idir == 2) ? nbr[side] : Ipz;
		const int dst = (ipx * Npy * Npz) + (ipy * Npz) + ipz;

		MPI_Sendrecv(sbuf, (int)num, MPI_REAL_T, dst, tag,
		             rbuf, (int)num, MPI_REAL_T, dst, tag, MPI_COMM_WORLD, &status);

		unpack(idir, irecv[side], e, a0, a1, b0, b1, size);
	}
}

#endif  // _MPI


// X 面の Ex を交換する (updateTpaEy / updateTpaEz が Ex[i-1] を読む)
void comm_cuda_E_X(void)
{
#ifdef _MPI
	comm_cuda_E(0, Ex, Npx, Ipx, iMin, iMax, jMin, jMax, kMin, kMax);
#endif
}


// Y 面の Ey を交換する (updateTpaEx / updateTpaEz が Ey[j-1] を読む)
void comm_cuda_E_Y(void)
{
#ifdef _MPI
	comm_cuda_E(1, Ey, Npy, Ipy, jMin, jMax, iMin, iMax, kMin, kMax);
#endif
}


// Z 面の Ez を交換する (updateTpaEx / updateTpaEy が Ez[k-1] を読む)
void comm_cuda_E_Z(void)
{
#ifdef _MPI
	comm_cuda_E(2, Ez, Npz, Ipz, kMin, kMax, iMin, iMax, jMin, jMax);
#endif
}
