/*
updateTpa.c

二光子吸収 (TPA) 非線形減衰
出典 : Honda, Shoji, Amemiya, "Proposal and analysis of an optical activation
       function based on metamaterial-loaded Si waveguides for nanophotonic
       neural networks", Opt. Lett. 49, 5811 (2024) (Si: β = 424 cm/GW)

物理モデルと等価非線形導電率 σ_NL の導出
----------------------------------------
TPA は強度依存吸収 dI/dz = -β I^2 (α_TPA = β I) を与える。
これを FDTD では、E 更新完了後に各 E 成分へ一次減衰

    E *= 1 / (1 + γ),   γ = Δt σ_NL / (2 ε0 εr)

を乗ずる陽的・無条件安定な方式で表す (γ >= 0 なのでエネルギーは単調減少)。
この操作は連続極限で ε ∂E/∂t = ∇×H - (σ_NL/2) E に相当する。すなわち
「実効導電率は σ_NL/2」であり、標準的な半陰的伝導項 (1-a)/(1+a) の σ とは
係数が 2 倍異なる点に注意 (下の導出はこの 2 倍を含めて厳密に行う)。

σ_NL(t) = A |E(t)|^2 (瞬時値) とおき、CW 定常でサイクル平均が
dI/dz = -β I^2 を再現するよう A を決める:

 (1) 屈折率 n = √εr の媒質中の進行波 E(t) = E0 cos(ωt-kz) の
     サイクル平均強度は I = (1/2) n ε0 c E0^2。
 (2) 減衰項 -(σ_NL/2) E / ε = -(A E^2 / 2ε) E の cos^3 を調和平衡すると
     基本波成分は cos^3 = (3/4)cos + (1/4)cos3ω より 3/4 倍。
     ゆえに実効 (線形換算) 導電率は σ_eff = (1/2)(3/4) A E0^2 = (3/8) A E0^2。
 (3) 低損失媒質の強度減衰係数は α_I = σ_eff η0 / n  (η0 = 1/(ε0 c))。
     TPA の要求 α_I = β I を課すと
        (3/8) A E0^2 / (ε0 c n) = β (1/2) n ε0 c E0^2
     →  A = (4/3) β εr ε0^2 c^2          (εr = n^2)
 (4) したがって
        γ = Δt σ_NL / (2 ε0 εr) = (2/3) β ε0 c^2 Δt |E|^2
     となり、γ は εr に依存しない (媒質係数は σ_NL と分母で相殺)。

|E|^2 は Yee 格子の各 E 成分位置での colocated 近似
(例: Ex 位置では Ex^2 + [4点平均した Ey]^2 + [4点平均した Ez]^2) を用いる。
1D 検証 (単一偏波の平面波) では単一成分のみ非零なので厳密。

平面波入射 (散乱界定式化) では格納配列は散乱界なので、全電界
E_tot = E_scat + E_inc に減衰を適用し、E_scat = E_tot/(1+γ) - E_inc を
書き戻す。給電 (feed) 時は格納配列が全電界なのでそのまま減衰する。
検証は data/sample/tpa_slab.ofd (解析解 T = 1/(1 + β I0 L) との比較) を参照。
*/

#include "ofd.h"
#include "finc.h"

// material id 毎の β [m/W] テーブルを作成する (β=0 : TPA なし)
void setupTpa(void)
{
	if (NTpa <= 0) return;

	free(TpaBeta);
	TpaBeta = (double *)malloc(NMaterial * sizeof(double));
	memset(TpaBeta, 0, NMaterial * sizeof(double));
	for (int n = 0; n < NTpa; n++) {
		TpaBeta[Tpa[n].m] = Tpa[n].beta;
	}
}


// 入射平面波の波形 (振幅 1) : 各成分は ei[] を乗じて得る
static inline double finc_waveform(double x, double y, double z, double t)
{
	real_t fi, dfi;
	finc(x, y, z, t, Planewave.r0, Planewave.ri, 1, Planewave.ai, Dt, &fi, &dfi);
	return (double)fi;
}


static void updateTpaEx(double t, double cf)
{
	int i;
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for (    i = iMin; i <  iMax; i++) {
	for (int j = jMin; j <= jMax; j++) {
		int64_t n = NA(i, j, kMin);
		for (int k = kMin; k <= kMax; k++, n++) {
			const double b = TpaBeta[iEx[n]];
			if (b <= 0) continue;

			// colocated 近似 : Ex 位置 (Xc[i], Yn[j], Zn[k]) での |E|^2
			double ex = Ex[n];
			double ey = 0.25 * (Ey[n] + Ey[n + Ni] + Ey[n - Nj] + Ey[n + Ni - Nj]);
			double ez = 0.25 * (Ez[n] + Ez[n + Ni] + Ez[n - Nk] + Ez[n + Ni - Nk]);
			double f = 0;
			if (IPlanewave) {
				// 散乱界 -> 全電界 (入射波形は colocated 近似で同一点評価)
				f = finc_waveform(Xc[i], Yn[j], Zn[k], t);
				ex += f * Planewave.ei[0];
				ey += f * Planewave.ei[1];
				ez += f * Planewave.ei[2];
			}
			const double e2 = (ex * ex) + (ey * ey) + (ez * ez);
			const double gamma = cf * b * e2;  // = Δt σ_NL / (2 ε0 εr)
			Ex[n] = (real_t)((ex / (1 + gamma)) - (f * Planewave.ei[0]));
		}
	}
	}
}


static void updateTpaEy(double t, double cf)
{
	int i;
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for (    i = iMin; i <= iMax; i++) {
	for (int j = jMin; j <  jMax; j++) {
		int64_t n = NA(i, j, kMin);
		for (int k = kMin; k <= kMax; k++, n++) {
			const double b = TpaBeta[iEy[n]];
			if (b <= 0) continue;

			// colocated 近似 : Ey 位置 (Xn[i], Yc[j], Zn[k]) での |E|^2
			double ey = Ey[n];
			double ex = 0.25 * (Ex[n] + Ex[n - Ni] + Ex[n + Nj] + Ex[n - Ni + Nj]);
			double ez = 0.25 * (Ez[n] + Ez[n + Nj] + Ez[n - Nk] + Ez[n + Nj - Nk]);
			double f = 0;
			if (IPlanewave) {
				f = finc_waveform(Xn[i], Yc[j], Zn[k], t);
				ex += f * Planewave.ei[0];
				ey += f * Planewave.ei[1];
				ez += f * Planewave.ei[2];
			}
			const double e2 = (ex * ex) + (ey * ey) + (ez * ez);
			const double gamma = cf * b * e2;
			Ey[n] = (real_t)((ey / (1 + gamma)) - (f * Planewave.ei[1]));
		}
	}
	}
}


static void updateTpaEz(double t, double cf)
{
	int i;
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for (    i = iMin; i <= iMax; i++) {
	for (int j = jMin; j <= jMax; j++) {
		int64_t n = NA(i, j, kMin);
		for (int k = kMin; k <  kMax; k++, n++) {
			const double b = TpaBeta[iEz[n]];
			if (b <= 0) continue;

			// colocated 近似 : Ez 位置 (Xn[i], Yn[j], Zc[k]) での |E|^2
			double ez = Ez[n];
			double ex = 0.25 * (Ex[n] + Ex[n - Ni] + Ex[n + Nk] + Ex[n - Ni + Nk]);
			double ey = 0.25 * (Ey[n] + Ey[n - Nj] + Ey[n + Nk] + Ey[n - Nj + Nk]);
			double f = 0;
			if (IPlanewave) {
				f = finc_waveform(Xn[i], Yn[j], Zc[k], t);
				ex += f * Planewave.ei[0];
				ey += f * Planewave.ei[1];
				ez += f * Planewave.ei[2];
			}
			const double e2 = (ex * ex) + (ey * ey) + (ez * ez);
			const double gamma = cf * b * e2;
			Ez[n] = (real_t)((ez / (1 + gamma)) - (f * Planewave.ei[2]));
		}
	}
	}
}


// E 更新後に呼ぶ (t は E の時刻 = (itime+1)*Dt)
void updateTpa(double t)
{
	if ((NTpa <= 0) || (TpaBeta == NULL)) return;

	// γ = (2/3) β ε0 c^2 Δt |E|^2 の |E|^2 を除く係数 (導出は冒頭コメント)
	const double cf = (2.0 / 3.0) * EPS0 * C * C * Dt;

	updateTpaEx(t, cf);
	updateTpaEy(t, cf);
	updateTpaEz(t, cf);
}
