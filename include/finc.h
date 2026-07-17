/*
finc.h

incidence function (gauss derivative)
波形 : 既定はガウス微分パルス。waveamp キー指定時 (WaveAmp > 0) は
CW 正弦波 (振幅 WaveAmp [V/m], 角周波数 WaveOmega) に切り替わる。
TPA 等の非線形解析では振幅の絶対値が意味を持つため CW 波源を用いる。
注意 : 本ヘッダは ofd.h (WaveAmp/WaveOmega の宣言) の後に include すること。
*/

#include <math.h>

// CW 立ち上げ周期数 (raised-cosine ランプ : 立ち上がり過渡の帯域を抑える)
#define FINC_CW_RAMP_CYCLES (5.0)

static inline void finc(
	double x, double y, double z, double t,
	const double r0[], const double ri[], double fc, double ai, double dt,
	real_t *fi, real_t *dfi)
{
	const double c = 2.99792458e8;

	t -= ((x - r0[0]) * ri[0]
	    + (y - r0[1]) * ri[1]
	    + (z - r0[2]) * ri[2]) / c;

	// CW 正弦波源 (waveamp キー指定時)
	// f(t) = E0 * s(t) * sin(ωt), s(t) : raised-cosine ランプ (0 <= t <= tr)
	// dfi は既存のパルス波源と同じく dt * df/dt を返す
	if (WaveAmp > 0) {
		if (t <= 0) {
			*fi  = 0;
			*dfi = 0;
			return;
		}
		const double tr = FINC_CW_RAMP_CYCLES * (2 * PI / WaveOmega);
		double s, ds;
		if (t < tr) {
			s  = 0.5 * (1 - cos(PI * t / tr));
			ds = 0.5 * (PI / tr) * sin(PI * t / tr);
		}
		else {
			s  = 1;
			ds = 0;
		}
		const double sn = sin(WaveOmega * t);
		const double cs = cos(WaveOmega * t);
		*fi  = (real_t)(WaveAmp * fc * s * sn);
		*dfi = (real_t)(dt * WaveAmp * fc * ((ds * sn) + (s * WaveOmega * cs)));
		return;
	}

	const double at = ai * t;
	const double ex = (at * at < 16) ? exp(-at * at) : 0;
	//const double ex = exp(-at * at);
	*fi = (real_t)(at * ex * fc);
	*dfi = (real_t)(dt * ai * (1 - 2 * at * at) * ex * fc);
}
