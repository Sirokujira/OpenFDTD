/*
updateTpa.cu

二光子吸収 (TPA) 非線形減衰 (GPU)

物理モデルと係数の導出は sol/updateTpa.c 冒頭のコメントを参照 (同一の式)。
  E *= 1 / (1 + γ),  γ = (2/3) β ε0 c^2 Δt |E|^2
|E|^2 は Yee 格子の colocated 近似 (成分位置で他成分を 4 点平均) で評価する。

CPU 版 (sol/updateTpa.c) と対になるファイルで、CUDA ビルドではこちらが
リンクされる (CMakeLists の SOURCES2 に sol/updateTpa.c を入れないこと)。
セル毎の演算は __host__ __device__ の共通関数にしてあるので、GPU 実行時も
CPU 実行時 (GPU = 0) も同じ式を通る。
*/

#include "ofd.h"
#include "ofd_cuda.h"
#include "ofd_prototype.h"   // cuda_malloc / cuda_memcpy
#include "finc_cuda.h"

// material id 毎の β [m/W] テーブル (デバイス側)
static double *d_TpaBeta = NULL;


// 入射平面波の波形 (振幅 1)
__host__ __device__ __forceinline__
static real_t tpa_finc(real_t x, real_t y, real_t z, real_t t, const param_t *p)
{
	real_t fi = 0, dfi = 0;
	finc_cuda_cw(x, y, z, t, p->r0, p->ri, (real_t)1, p->ai, p->dt, p->waveAmp, p->waveOmega, &fi, &dfi);

	return fi;
}


// Ex 位置 1 セル分
__host__ __device__
static void updateTpaEx_(
	int i, int j, int k, real_t ex[], const real_t ey[], const real_t ez[],
	const id_t iex[], const double beta[], double cf, real_t t,
	const real_t xc[], const real_t yn[], const real_t zn[], const param_t *p)
{
	const int64_t n = LA(p, i, j, k);
	const double b = beta[iex[n]];
	if (b <= 0) return;

	double vx = ex[n];
	double vy = 0.25 * (ey[n] + ey[n + p->Ni] + ey[n - p->Nj] + ey[n + p->Ni - p->Nj]);
	double vz = 0.25 * (ez[n] + ez[n + p->Ni] + ez[n - p->Nk] + ez[n + p->Ni - p->Nk]);
	double f = 0;
	if (p->IPlanewave) {
		// 散乱界 -> 全電界
		f = tpa_finc(xc[i], yn[j], zn[k], t, p);
		vx += f * p->ei[0];
		vy += f * p->ei[1];
		vz += f * p->ei[2];
	}
	const double e2 = (vx * vx) + (vy * vy) + (vz * vz);
	const double gamma = cf * b * e2;
	ex[n] = (real_t)((vx / (1 + gamma)) - (f * p->ei[0]));
}


// Ey 位置 1 セル分
__host__ __device__
static void updateTpaEy_(
	int i, int j, int k, real_t ey[], const real_t ex[], const real_t ez[],
	const id_t iey[], const double beta[], double cf, real_t t,
	const real_t xn[], const real_t yc[], const real_t zn[], const param_t *p)
{
	const int64_t n = LA(p, i, j, k);
	const double b = beta[iey[n]];
	if (b <= 0) return;

	double vy = ey[n];
	double vx = 0.25 * (ex[n] + ex[n - p->Ni] + ex[n + p->Nj] + ex[n - p->Ni + p->Nj]);
	double vz = 0.25 * (ez[n] + ez[n + p->Nj] + ez[n - p->Nk] + ez[n + p->Nj - p->Nk]);
	double f = 0;
	if (p->IPlanewave) {
		f = tpa_finc(xn[i], yc[j], zn[k], t, p);
		vx += f * p->ei[0];
		vy += f * p->ei[1];
		vz += f * p->ei[2];
	}
	const double e2 = (vx * vx) + (vy * vy) + (vz * vz);
	const double gamma = cf * b * e2;
	ey[n] = (real_t)((vy / (1 + gamma)) - (f * p->ei[1]));
}


// Ez 位置 1 セル分
__host__ __device__
static void updateTpaEz_(
	int i, int j, int k, real_t ez[], const real_t ex[], const real_t ey[],
	const id_t iez[], const double beta[], double cf, real_t t,
	const real_t xn[], const real_t yn[], const real_t zc[], const param_t *p)
{
	const int64_t n = LA(p, i, j, k);
	const double b = beta[iez[n]];
	if (b <= 0) return;

	double vz = ez[n];
	double vx = 0.25 * (ex[n] + ex[n - p->Ni] + ex[n + p->Nk] + ex[n - p->Ni + p->Nk]);
	double vy = 0.25 * (ey[n] + ey[n - p->Nj] + ey[n + p->Nk] + ey[n - p->Nj + p->Nk]);
	double f = 0;
	if (p->IPlanewave) {
		f = tpa_finc(xn[i], yn[j], zc[k], t, p);
		vx += f * p->ei[0];
		vy += f * p->ei[1];
		vz += f * p->ei[2];
	}
	const double e2 = (vx * vx) + (vy * vy) + (vz * vz);
	const double gamma = cf * b * e2;
	ez[n] = (real_t)((vz / (1 + gamma)) - (f * p->ei[2]));
}


// gpu
__global__
static void updateTpaEx_gpu(
	real_t ex[], const real_t ey[], const real_t ez[], const id_t iex[],
	const double beta[], double cf, real_t t,
	const real_t xc[], const real_t yn[], const real_t zn[])
{
	const int i = d_Param.iMin + threadIdx.z + (blockIdx.z * blockDim.z);
	const int j = d_Param.jMin + threadIdx.y + (blockIdx.y * blockDim.y);
	const int k = d_Param.kMin + threadIdx.x + (blockIdx.x * blockDim.x);
	if ((i <  d_Param.iMax) &&
	    (j <= d_Param.jMax) &&
	    (k <= d_Param.kMax)) {
		updateTpaEx_(i, j, k, ex, ey, ez, iex, beta, cf, t, xc, yn, zn, &d_Param);
	}
}


__global__
static void updateTpaEy_gpu(
	real_t ey[], const real_t ex[], const real_t ez[], const id_t iey[],
	const double beta[], double cf, real_t t,
	const real_t xn[], const real_t yc[], const real_t zn[])
{
	const int i = d_Param.iMin + threadIdx.z + (blockIdx.z * blockDim.z);
	const int j = d_Param.jMin + threadIdx.y + (blockIdx.y * blockDim.y);
	const int k = d_Param.kMin + threadIdx.x + (blockIdx.x * blockDim.x);
	if ((i <= d_Param.iMax) &&
	    (j <  d_Param.jMax) &&
	    (k <= d_Param.kMax)) {
		updateTpaEy_(i, j, k, ey, ex, ez, iey, beta, cf, t, xn, yc, zn, &d_Param);
	}
}


__global__
static void updateTpaEz_gpu(
	real_t ez[], const real_t ex[], const real_t ey[], const id_t iez[],
	const double beta[], double cf, real_t t,
	const real_t xn[], const real_t yn[], const real_t zc[])
{
	const int i = d_Param.iMin + threadIdx.z + (blockIdx.z * blockDim.z);
	const int j = d_Param.jMin + threadIdx.y + (blockIdx.y * blockDim.y);
	const int k = d_Param.kMin + threadIdx.x + (blockIdx.x * blockDim.x);
	if ((i <= d_Param.iMax) &&
	    (j <= d_Param.jMax) &&
	    (k <  d_Param.kMax)) {
		updateTpaEz_(i, j, k, ez, ex, ey, iez, beta, cf, t, xn, yn, zc, &d_Param);
	}
}


// cpu (GPU = 0 のとき)
static void updateTpa_cpu(double cf, real_t t)
{
	int i;
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for (    i = iMin; i <  iMax; i++) {
	for (int j = jMin; j <= jMax; j++) {
	for (int k = kMin; k <= kMax; k++) {
		updateTpaEx_(i, j, k, Ex, Ey, Ez, iEx, TpaBeta, cf, t, h_Xc, h_Yn, h_Zn, &h_Param);
	}
	}
	}

#ifdef _OPENMP
#pragma omp parallel for
#endif
	for (    i = iMin; i <= iMax; i++) {
	for (int j = jMin; j <  jMax; j++) {
	for (int k = kMin; k <= kMax; k++) {
		updateTpaEy_(i, j, k, Ey, Ex, Ez, iEy, TpaBeta, cf, t, h_Xn, h_Yc, h_Zn, &h_Param);
	}
	}
	}

#ifdef _OPENMP
#pragma omp parallel for
#endif
	for (    i = iMin; i <= iMax; i++) {
	for (int j = jMin; j <= jMax; j++) {
	for (int k = kMin; k <  kMax; k++) {
		updateTpaEz_(i, j, k, Ez, Ex, Ey, iEz, TpaBeta, cf, t, h_Xn, h_Yn, h_Zc, &h_Param);
	}
	}
	}
}


// material id 毎の β テーブルを作り、GPU 実行時はデバイスへ転送する
void setupTpa(void)
{
	if (NTpa <= 0) return;

	free(TpaBeta);
	TpaBeta = (double *)malloc(NMaterial * sizeof(double));
	memset(TpaBeta, 0, NMaterial * sizeof(double));
	for (int n = 0; n < NTpa; n++) {
		TpaBeta[Tpa[n].m] = Tpa[n].beta;
	}

	if (GPU) {
		const size_t size = NMaterial * sizeof(double);
		cuda_malloc(GPU, UM, (void **)&d_TpaBeta, size);
		cuda_memcpy(GPU, d_TpaBeta, TpaBeta, size, cudaMemcpyHostToDevice);
	}
}


// E 更新後に呼ぶ (t は E の時刻 = (itime+1)*Dt)
void updateTpa(double t)
{
	if ((NTpa <= 0) || (TpaBeta == NULL)) return;

	// γ = (2/3) β ε0 c^2 Δt |E|^2 の |E|^2 を除く係数 (導出は sol/updateTpa.c)
	const double cf = (2.0 / 3.0) * EPS0 * C * C * Dt;

	if (GPU) {
		cudaMemcpyToSymbol(d_Param, &h_Param, sizeof(param_t));

		dim3 gridEx(
			CEIL(kMax - kMin + 1, updateBlock.x),
			CEIL(jMax - jMin + 1, updateBlock.y),
			CEIL(iMax - iMin + 0, updateBlock.z));
		updateTpaEx_gpu<<<gridEx, updateBlock>>>(
			Ex, Ey, Ez, d_iEx, d_TpaBeta, cf, (real_t)t, d_Xc, d_Yn, d_Zn);

		dim3 gridEy(
			CEIL(kMax - kMin + 1, updateBlock.x),
			CEIL(jMax - jMin + 0, updateBlock.y),
			CEIL(iMax - iMin + 1, updateBlock.z));
		updateTpaEy_gpu<<<gridEy, updateBlock>>>(
			Ey, Ex, Ez, d_iEy, d_TpaBeta, cf, (real_t)t, d_Xn, d_Yc, d_Zn);

		dim3 gridEz(
			CEIL(kMax - kMin + 0, updateBlock.x),
			CEIL(jMax - jMin + 1, updateBlock.y),
			CEIL(iMax - iMin + 1, updateBlock.z));
		updateTpaEz_gpu<<<gridEz, updateBlock>>>(
			Ez, Ex, Ey, d_iEz, d_TpaBeta, cf, (real_t)t, d_Xn, d_Yn, d_Zc);

		if (UM) cudaDeviceSynchronize();
	}
	else {
		updateTpa_cpu(cf, (real_t)t);
	}
}
