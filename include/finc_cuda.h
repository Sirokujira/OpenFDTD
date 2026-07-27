/*
finc_cuda.h (CUDA)

incidence function
波形 : 既定はガウス微分パルス。waveamp キー指定時 (WaveAmp > 0) は
CW 正弦波 (振幅 WaveAmp [V/m], 角周波数 WaveOmega) に切り替わる。
CPU 版 include/finc.h と同じ式を実装する (TPA 等の非線形解析では
振幅の絶対値が意味を持つため、両者がずれてはいけない)。

CW のパラメータはホスト側グローバル (WaveAmp/WaveOmega) をデバイスから
参照できないので、param_t の waveAmp/waveOmega 経由で渡す。
*/

// CW 立ち上げ周期数 (raised-cosine ランプ : 立ち上がり過渡の帯域を抑える)
// CPU 版 include/finc.h の FINC_CW_RAMP_CYCLES と同じ値にすること
#define FINC_CUDA_CW_RAMP_CYCLES ((real_t)5.0)

__host__ __device__ __forceinline__
static void finc_cuda_cw(
	real_t x, real_t y, real_t z, real_t t,
	const real_t r0[], const real_t ri[], real_t fc, real_t ai, real_t dt,
	real_t wamp, real_t womega,
	real_t *fi, real_t *dfi)
{
	const real_t c = (real_t)2.99792458e8;
	const real_t pi = (real_t)3.14159265358979323846;

	t -= ((x - r0[0]) * ri[0]
	    + (y - r0[1]) * ri[1]
	    + (z - r0[2]) * ri[2]) / c;

	// CW 正弦波源 (waveamp キー指定時)
	// f(t) = E0 * s(t) * sin(ωt), s(t) : raised-cosine ランプ (0 <= t <= tr)
	if (wamp > 0) {
		if (t <= 0) {
			*fi  = 0;
			*dfi = 0;
			return;
		}
		const real_t tr = FINC_CUDA_CW_RAMP_CYCLES * (2 * pi / womega);
		real_t s, ds;
		if (t < tr) {
			s  = (real_t)0.5 * (1 - (real_t)cos(pi * t / tr));
			ds = (real_t)0.5 * (pi / tr) * (real_t)sin(pi * t / tr);
		}
		else {
			s  = 1;
			ds = 0;
		}
		const real_t sn = (real_t)sin(womega * t);
		const real_t cs = (real_t)cos(womega * t);
		*fi  = wamp * fc * s * sn;
		*dfi = dt * wamp * fc * ((ds * sn) + (s * womega * cs));
		return;
	}

	const real_t at = ai * t;
	const real_t ex = (at * at < 16) ? (real_t)exp(-at * at) : 0;
	*fi = at * ex * fc;
	*dfi = dt * ai * (1 - 2 * at * at) * ex * fc;
}


// 既存の呼び出し互換 : param_t を持たない箇所向け (CW 無効 = 従来のパルス)
__host__ __device__ __forceinline__
static void finc_cuda(
	real_t x, real_t y, real_t z, real_t t,
	const real_t r0[], const real_t ri[], real_t fc, real_t ai, real_t dt,
	real_t *fi, real_t *dfi)
{
	finc_cuda_cw(x, y, z, t, r0, ri, fc, ai, dt, 0, 0, fi, dfi);
}
