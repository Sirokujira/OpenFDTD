/*
test_pml.c

PML 吸収境界 (abc = 1 L M R0) のテスト。

対象:
  setupPml()        sol/setupPml.c       損失分布 gPml[XYZ][nc] と材料係数 rPml*
  setupPml[EH][xyz] sol/setupPml[EH]?.c  PML 内の格子点リストと境界材料
  pml[EH][xyz]()    sol/pml[EH]?.c       分割界 (split field) の更新

物理での裏付けは data/sample/pml_reflection_check.sh (実際に反射が小さい
こと) が受け持つ。あちらは符号の誤りや勾配の欠落のような大きな壊れ方しか
捕まえられない (正常な PML でも角と近傍界で 1% 程度の差が残り、
係数のわずかな誤りが埋もれる) ので、ここで式の単位で押さえる。

期待値の出所:
  - 損失分布: 多項式分布の PML の理論反射率
        R(0) = exp(-2 ∫[0,L] σ(x)/(ε0 c) dx),  σ(x) = σmax (x/L)^M
    を、分布の標本から σmax を取り出して計算し、入力の R0 に戻ること。
    更新式 X/(1 + g) の g は σ dt/ε0 なので、セル単位では
    σ/(ε0 c) [1/セル] = g d / (c dt)。
    標本の位置は「格子点の深さより半セル内側」(E 節点の深さ p に
    σ(p - 1/2)、H セル中心の深さ p + 1/2 に σ(p))。これは既存の設計で、
    E と H を深さ順に並べると狭義単調増加になることで位置ずれを検出する。
  - PML 内のセル幅は境界のセル幅 (Xn[1]-Xn[0] / Xn[Nx]-Xn[Nx-1]) を延長
    したもの。不等間隔メッシュにして、低い側と高い側で別の幅にする。
  - 格子点の範囲: 外周の PEC 壁 (深さ L) の接線 E は更新しない。
    Ex なら x はセル i ∈ [-L, Nx+L-1]、y,z は節点 [-L+1, N+L-1]。
  - 境界材料: PML 内の点は最も近い内部の点の材料を引き継ぐ。
  - 更新式: Maxwell の回転
        ∂Ea/∂t =  (1/ε)(∂Hc/∂b - ∂Hb/∂c),  ∂Ha/∂t = -(1/μ)(∂Ec/∂b - ∂Eb/∂c)
    ((a,b,c) は (x,y,z) の巡回) を b 方向と c 方向の分割界に分け、
    それぞれ自分の方向の損失 g/√(εμ) で割る (E と H が同じ減衰率になる
    整合条件)。E は H の後退差分と節点の双対セル幅、H は E の前進差分と
    セル幅を使う。
  - 6 成分が同じ式であることを、Nx≠Ny≠Nz・不等間隔・非真空の材料
    (εr≠1, μr≠1) で確かめる。対称な設定では軸や係数の取り違えが隠れる。
*/

#include "ofd.h"
#include "ofd_prototype.h"
#include "ofd_test.h"

#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#define L     4
#define ORDER 2.0
#define R0    1e-4

/* 各軸のセル幅 [m] (不等間隔、境界セルは両端で別の幅) */
static const double DX[] = {1.0e-3, 1.3e-3, 0.8e-3, 1.1e-3, 1.6e-3};
static const double DY[] = {0.7e-3, 0.9e-3, 1.2e-3, 1.0e-3, 0.8e-3, 1.2e-3};
static const double DZ[] = {1.5e-3, 1.0e-3, 0.9e-3, 1.1e-3, 1.3e-3, 0.95e-3, 0.6e-3};

/* 軸ごとの情報 (0/1/2 = x/y/z) */
static int     n_[3];
static double *node_[3], *cell_[3];
static real_t *gn_[3], *gc_[3];

static real_t **efield[3], **hfield[3];

static double cdt;

static void make_axis(int a, const double *d, int n, double **nodes, double **cells)
{
	double *xn = (double *)malloc((n + 1) * sizeof(double));
	double *xc = (double *)malloc(n * sizeof(double));
	xn[0] = -0.5e-3 * (a + 1);
	for (int i = 0; i < n; i++) {
		xn[i + 1] = xn[i] + d[i];
		xc[i] = (xn[i] + xn[i + 1]) / 2;
	}
	*nodes = xn;
	*cells = xc;
}

/* 境界セル幅を延長した、節点 (双対セル) の幅 */
static double dual_width(int a, int i)
{
	const int n = n_[a];
	const double *x = node_[a];
	if (i <= 0) return x[1] - x[0];
	if (i >= n) return x[n] - x[n - 1];
	return (x[i + 1] - x[i - 1]) / 2;
}

/* 境界セル幅を延長した、セル (主セル) の幅 */
static double cell_width(int a, int i)
{
	const int n = n_[a];
	const double *x = node_[a];
	if (i < 0) return x[1] - x[0];
	if (i >= n) return x[n] - x[n - 1];
	return x[i + 1] - x[i];
}

static uint32_t rng = 12345;
static double urand(void)
{
	rng = rng * 1664525u + 1013904223u;
	return ((rng >> 8) & 0xffff) / 65536.0 - 0.5;
}

static void pml_setup(void)
{
	n_[0] = Nx = (int)(sizeof(DX) / sizeof(DX[0]));
	n_[1] = Ny = (int)(sizeof(DY) / sizeof(DY[0]));
	n_[2] = Nz = (int)(sizeof(DZ) / sizeof(DZ[0]));
	iMin = jMin = kMin = 0;
	iMax = Nx; jMax = Ny; kMax = Nz;

	make_axis(0, DX, Nx, &Xn, &Xc);
	make_axis(1, DY, Ny, &Yn, &Yc);
	make_axis(2, DZ, Nz, &Zn, &Zc);
	node_[0] = Xn; node_[1] = Yn; node_[2] = Zn;
	cell_[0] = Xc; cell_[1] = Yc; cell_[2] = Zc;

	iABC = 1;
	cPML.l = L;
	cPML.m = ORDER;
	cPML.r0 = R0;

	Dt = 1.0e-12;
	cdt = C * Dt;

	/* setupSize() と同じ配列レイアウト (PML 分の余白 L 層) */
	Nk = 1;
	Nj = (kMax - kMin) + (2 * L) + 1;
	Ni = ((jMax - jMin) + (2 * L) + 1) * Nj;
	N0 = -((iMin - L) * Ni + (jMin - L) * Nj + (kMin - L) * Nk);
	NN = NA(iMax + L, jMax + L, kMax + L) + 1;

	iEx = (id_t *)calloc(NN, sizeof(id_t));
	iEy = (id_t *)calloc(NN, sizeof(id_t));
	iEz = (id_t *)calloc(NN, sizeof(id_t));
	iHx = (id_t *)calloc(NN, sizeof(id_t));
	iHy = (id_t *)calloc(NN, sizeof(id_t));
	iHz = (id_t *)calloc(NN, sizeof(id_t));
	Ex = (real_t *)calloc(NN, sizeof(real_t));
	Ey = (real_t *)calloc(NN, sizeof(real_t));
	Ez = (real_t *)calloc(NN, sizeof(real_t));
	Hx = (real_t *)calloc(NN, sizeof(real_t));
	Hy = (real_t *)calloc(NN, sizeof(real_t));
	Hz = (real_t *)calloc(NN, sizeof(real_t));
	efield[0] = &Ex; efield[1] = &Ey; efield[2] = &Ez;
	hfield[0] = &Hx; hfield[1] = &Hy; hfield[2] = &Hz;

	/* 0=真空, 1=PEC, 2,3=磁性誘電体 (εr, μr ≠ 1) */
	NMaterial = 4;
	Material = (material_t *)calloc(NMaterial, sizeof(material_t));
	for (int m = 0; m < NMaterial; m++) {
		Material[m].type = 1;
		Material[m].epsr = Material[m].amur = 1;
	}
	Material[2].epsr = 4;   Material[2].amur = 2.25;   /* 1/√(εμ) = 1/3 */
	Material[3].epsr = 2;   Material[3].amur = 1.5;

	/* 内部の材料は位置で 2/3 を混ぜる (どの内部点を引き継いだか分かるように) */
	for (int i = 0; i <= Nx; i++) {
	for (int j = 0; j <= Ny; j++) {
	for (int k = 0; k <= Nz; k++) {
		const id_t m = (id_t)(2 + ((i + 2 * j + 3 * k) % 2));
		const id_t mh = (id_t)(2 + ((2 * i + j + 5 * k) % 2));
		iEx[NA(i, j, k)] = iEy[NA(i, j, k)] = iEz[NA(i, j, k)] = m;
		iHx[NA(i, j, k)] = iHy[NA(i, j, k)] = iHz[NA(i, j, k)] = mh;
	}
	}
	}

	/* 節点 (RXn) は双対セル幅、セル (RXc) はセル幅。端の節点は境界セル幅
	   (sol/setup.c の setup_mesh と同じ定義を、ここでは幅の関数から作る) */
	real_t **rn[3] = {&RXn, &RYn, &RZn};
	real_t **rc[3] = {&RXc, &RYc, &RZc};
	real_t **gn[3] = {&gPmlXn, &gPmlYn, &gPmlZn};
	real_t **gc[3] = {&gPmlXc, &gPmlYc, &gPmlZc};
	for (int a = 0; a < 3; a++) {
		const int n = n_[a];
		*rn[a] = (real_t *)malloc((n + 1) * sizeof(real_t));
		*rc[a] = (real_t *)malloc(n * sizeof(real_t));
		for (int i = 0; i <= n; i++) (*rn[a])[i] = (real_t)(cdt / dual_width(a, i));
		for (int i = 0; i <  n; i++) (*rc[a])[i] = (real_t)(cdt / cell_width(a, i));
		*gn[a] = (real_t *)calloc(n + 2 * L + 1, sizeof(real_t));
		*gc[a] = (real_t *)calloc(n + 2 * L + 1, sizeof(real_t));
		gn_[a] = *gn[a];
		gc_[a] = *gc[a];
	}

	rPmlE = (real_t *)calloc(NMaterial, sizeof(real_t));
	rPmlH = (real_t *)calloc(NMaterial, sizeof(real_t));
	rPml  = (real_t *)calloc(NMaterial, sizeof(real_t));

	setupPml();
}

static void pml_teardown(void)
{
	free(Xn); free(Yn); free(Zn); free(Xc); free(Yc); free(Zc);
	free(iEx); free(iEy); free(iEz); free(iHx); free(iHy); free(iHz);
	free(Ex); free(Ey); free(Ez); free(Hx); free(Hy); free(Hz);
	free(RXn); free(RYn); free(RZn); free(RXc); free(RYc); free(RZc);
	free(gPmlXn); free(gPmlYn); free(gPmlZn); free(gPmlXc); free(gPmlYc); free(gPmlZc);
	free(rPmlE); free(rPmlH); free(rPml);
	free(fPmlEx); free(fPmlEy); free(fPmlEz); free(fPmlHx); free(fPmlHy); free(fPmlHz);
	free(Exy); free(Exz); free(Eyz); free(Eyx); free(Ezx); free(Ezy);
	free(Hxy); free(Hxz); free(Hyz); free(Hyx); free(Hzx); free(Hzy);
	free(Material);

	Xn = Yn = Zn = Xc = Yc = Zc = NULL;
	iEx = iEy = iEz = iHx = iHy = iHz = NULL;
	Ex = Ey = Ez = Hx = Hy = Hz = NULL;
	RXn = RYn = RZn = RXc = RYc = RZc = NULL;
	gPmlXn = gPmlYn = gPmlZn = gPmlXc = gPmlYc = gPmlZc = NULL;
	rPmlE = rPmlH = rPml = NULL;
	fPmlEx = fPmlEy = fPmlEz = fPmlHx = fPmlHy = fPmlHz = NULL;
	Exy = Exz = Eyz = Eyx = Ezx = Ezy = NULL;
	Hxy = Hxz = Hyz = Hyx = Hzx = Hzy = NULL;
	numPmlEx = numPmlEy = numPmlEz = numPmlHx = numPmlHy = numPmlHz = 0;
	Material = NULL;
	NMaterial = 0;
	iABC = 0;
	memset(&cPML, 0, sizeof(cPML));
}

/* ---- 損失分布 ---- */

/* 正規化した損失 σ/(ε0 c) [1/セル] : g d / (c dt)。
   深さ 2*h 半セル (h = 1, 2, ...) の標本を返す。side 0/1 = 低い側/高い側 */
static double sigma_at(int a, int side, int h)
{
	const int n = n_[a];
	const double *x = node_[a];
	const double d = side ? (x[n] - x[n - 1]) : (x[1] - x[0]);
	real_t g;
	if (h % 2) {
		/* 奇数 : H (セル中心), 深さ h/2 = p + 1/2 → セル i = -p-1 / n+p */
		const int p = h / 2;
		const int i = side ? (n + p) : (-p - 1);
		g = gc_[a][i + L];
	}
	else {
		/* 偶数 : E (節点), 深さ p = h/2 → 節点 i = -p / n+p */
		const int p = h / 2;
		const int i = side ? (n + p) : (-p);
		g = gn_[a][i + L];
	}
	return g * d / cdt;
}

static void test_profile(void)
{
	for (int a = 0; a < 3; a++) {
		const int n = n_[a];

		/* 内部 (境界の節点を含む) は損失 0 */
		int zero = 1;
		for (int i = 0; i <= n; i++) zero &= (gn_[a][i + L] == 0);
		for (int i = 0; i <  n; i++) zero &= (gc_[a][i + L] == 0);
		CHECK(zero);

		/* 深さ 1/2 .. L-1/2 を半セル刻みで並べる (h = 1 .. 2L-1)。
		   深さ 1/2 の H は 0、その先は狭義単調増加、両側で一致 */
		for (int side = 0; side < 2; side++) {
			CHECK(sigma_at(a, side, 1) == 0);
			int inc = 1;
			for (int h = 2; h <= 2 * L - 1; h++) {
				inc &= (sigma_at(a, side, h) > sigma_at(a, side, h - 1));
			}
			CHECK(inc);
		}
		double sym = 0;
		for (int h = 1; h <= 2 * L - 1; h++) {
			const double lo = sigma_at(a, 0, h), hi = sigma_at(a, 1, h);
			sym = fmax(sym, fabs(lo - hi) / fmax(fabs(lo), 1e-30));
		}
		CHECK(sym < 1e-5);

		/* 多項式分布 σ(x) = σmax (x/L)^M を半セル内側で標本化したもの */
		const double smax_sample = sigma_at(a, 0, 2 * L - 1);   /* x = L - 1 */
		const double smax = smax_sample / pow((L - 1.0) / L, ORDER);
		double shape = 0;
		for (int h = 2; h <= 2 * L - 1; h++) {
			const double x = (h - 1) / 2.0;
			const double expect = smax * pow(x / L, ORDER);
			shape = fmax(shape, fabs(sigma_at(a, 0, h) - expect) / expect);
		}
		CHECK(shape < 1e-5);

		/* 理論反射率 R(0) = exp(-2 ∫[0,L] σ dx) = exp(-2 σmax L/(M+1)) が R0 */
		const double r = exp(-2 * smax * L / (ORDER + 1));
		CHECK_NEAR(log(r), log(R0), 1e-5 * fabs(log(R0)));
	}

	/* 軸の取り違え : 各軸の分布は自軸の境界セル幅で正規化されている */
	CHECK(fabs(sigma_at(0, 0, 4) - sigma_at(1, 1, 4)) < 1e-5 * sigma_at(0, 0, 4));
	CHECK(gPmlXn[-2 + L] != gPmlYn[-2 + L]);   /* 境界セル幅 1.0 mm と 0.7 mm */

	/* 材料係数 : 1/εr, 1/μr, 1/√(εr μr)。PEC は 0 */
	CHECK_NEAR(rPmlE[0], 1, 1e-7);
	CHECK_NEAR(rPmlH[0], 1, 1e-7);
	CHECK_NEAR(rPml[0], 1, 1e-7);
	CHECK(rPmlE[PEC] == 0 && rPmlH[PEC] == 0 && rPml[PEC] == 0);
	CHECK_NEAR(rPmlE[2], 0.25, 1e-7);
	CHECK_NEAR(rPmlH[2], 1 / 2.25, 1e-7);
	CHECK_NEAR(rPml[2], 1 / 3.0, 1e-7);
	CHECK_NEAR(rPml[3], 1 / sqrt(3.0), 1e-7);
}

/* ---- 格子点リスト ---- */

typedef struct {
	const char *name;
	void (*setup)(int);
	void (*update)(void);
	pml_t **list;
	int64_t *num;
	int isE;
	int a;              /* 成分の向き */
	real_t **split_b;   /* a の次の軸 b 方向の分割界 */
	real_t **split_c;   /* c 方向の分割界 */
	id_t **id;
} comp_t;

static comp_t comps[6];

static void make_comps(void)
{
	comp_t c[6] = {
		{"Ex", setupPmlEx, pmlEx, &fPmlEx, &numPmlEx, 1, 0, &Exy, &Exz, &iEx},
		{"Ey", setupPmlEy, pmlEy, &fPmlEy, &numPmlEy, 1, 1, &Eyz, &Eyx, &iEy},
		{"Ez", setupPmlEz, pmlEz, &fPmlEz, &numPmlEz, 1, 2, &Ezx, &Ezy, &iEz},
		{"Hx", setupPmlHx, pmlHx, &fPmlHx, &numPmlHx, 0, 0, &Hxy, &Hxz, &iHx},
		{"Hy", setupPmlHy, pmlHy, &fPmlHy, &numPmlHy, 0, 1, &Hyz, &Hyx, &iHy},
		{"Hz", setupPmlHz, pmlHz, &fPmlHz, &numPmlHz, 0, 2, &Hzx, &Hzy, &iHz},
	};
	memcpy(comps, c, sizeof(c));
}

/* 成分の格子点の範囲 (各軸) : 自分の向きは E ならセル・H なら節点、
   他の 2 軸はその逆。PML を含む範囲と内部の範囲 */
static void index_range(const comp_t *c, int b, int *lo, int *hi, int *ilo, int *ihi)
{
	const int n = n_[b];
	const int on_cell = c->isE ? (b == c->a) : (b != c->a);
	if (on_cell) { *lo = -L;     *hi = n + L - 1; *ilo = 0; *ihi = n - 1; }
	else         { *lo = -L + 1; *hi = n + L - 1; *ilo = 0; *ihi = n;     }
}

static int clampi(int v, int lo, int hi) { return (v < lo) ? lo : ((v > hi) ? hi : v); }

static void test_lists(void)
{
	for (int t = 0; t < 6; t++) {
		comp_t *c = &comps[t];
		int lo[3], hi[3], ilo[3], ihi[3];
		for (int b = 0; b < 3; b++) index_range(c, b, &lo[b], &hi[b], &ilo[b], &ihi[b]);

		c->setup(0);
		const int64_t all = (int64_t)(hi[0] - lo[0] + 1) * (hi[1] - lo[1] + 1) * (hi[2] - lo[2] + 1);
		const int64_t inner = (int64_t)(ihi[0] - ilo[0] + 1) * (ihi[1] - ilo[1] + 1) * (ihi[2] - ilo[2] + 1);
		const int64_t num = *c->num;
		CHECK(num == all - inner);
		if (num != all - inner) {
			printf("    %s : num = %lld, expected %lld\n", c->name, (long long)num, (long long)(all - inner));
		}

		*c->list = (pml_t *)malloc(num * sizeof(pml_t));
		c->setup(1);

		/* 範囲内・内部の外・重複なし・材料は最も近い内部の点 */
		char *seen = (char *)calloc(NN, 1);
		int inrange = 1, outside = 1, unique = 1, material = 1;
		for (int64_t n = 0; n < num; n++) {
			const pml_t *p = &(*c->list)[n];
			const int idx[3] = {p->i, p->j, p->k};
			int inside = 1;
			for (int b = 0; b < 3; b++) {
				inrange &= (idx[b] >= lo[b]) && (idx[b] <= hi[b]);
				inside &= (idx[b] >= ilo[b]) && (idx[b] <= ihi[b]);
			}
			outside &= !inside;
			if (!inrange) continue;
			const int64_t q = NA(p->i, p->j, p->k);
			unique &= !seen[q];
			seen[q] = 1;
			const int ci = clampi(p->i, ilo[0], ihi[0]);
			const int cj = clampi(p->j, ilo[1], ihi[1]);
			const int ck = clampi(p->k, ilo[2], ihi[2]);
			material &= (p->m == (*c->id)[NA(ci, cj, ck)]);
		}
		free(seen);
		CHECK(inrange);
		CHECK(outside);
		CHECK(unique);
		CHECK(material);
		if (!(inrange && outside && unique && material)) {
			printf("    %s : inrange=%d outside=%d unique=%d material=%d\n", c->name, inrange, outside, unique, material);
		}

		*c->split_b = (real_t *)malloc(num * sizeof(real_t));
		*c->split_c = (real_t *)malloc(num * sizeof(real_t));
	}
}

/* ---- 更新式 ---- */

static real_t fval(real_t **f, int i, int j, int k)
{
	return (*f)[NA(i, j, k)];
}

/* 成分 c の点 p での、軸 b 方向の分割界の増分と損失 g。
   E : +(1/ε) ∂Hc/∂b (分割 b), -(1/ε) ∂Hb/∂c (分割 c)  … H の後退差分
   H : -(1/μ) ∂Ec/∂b (分割 b), +(1/μ) ∂Eb/∂c (分割 c)  … E の前進差分 */
static void expect_split(const comp_t *c, const pml_t *p, int which, double *delta, double *g)
{
	const int a = c->a;
	const int b = (a + 1) % 3, cc = (a + 2) % 3;
	const int dir = which ? cc : b;       /* 差分をとる軸 */
	const int other = which ? b : cc;     /* 差分される成分の向き */
	const int idx[3] = {p->i, p->j, p->k};
	const material_t *m = &Material[p->m];
	int nb[3] = {idx[0], idx[1], idx[2]};
	double diff, width, loss;

	if (c->isE) {
		nb[dir] -= 1;
		diff = fval(hfield[other], idx[0], idx[1], idx[2]) - fval(hfield[other], nb[0], nb[1], nb[2]);
		width = dual_width(dir, idx[dir]);
		loss = gn_[dir][idx[dir] + L];
		*delta = (which ? -1 : +1) * (cdt / width) * diff / m->epsr;
	}
	else {
		nb[dir] += 1;
		diff = fval(efield[other], nb[0], nb[1], nb[2]) - fval(efield[other], idx[0], idx[1], idx[2]);
		width = cell_width(dir, idx[dir]);
		loss = gc_[dir][idx[dir] + L];
		*delta = (which ? +1 : -1) * (cdt / width) * diff / m->amur;
	}
	*g = loss / sqrt(m->epsr * m->amur);
}

static void test_update(void)
{
	/* 全域を乱数で埋める */
	for (int64_t q = 0; q < NN; q++) {
		Ex[q] = (real_t)urand(); Ey[q] = (real_t)urand(); Ez[q] = (real_t)urand();
		Hx[q] = (real_t)urand(); Hy[q] = (real_t)urand(); Hz[q] = (real_t)urand();
	}

	/* 内部の点は PML の更新で変わらないこと (番兵) */
	const real_t ex0 = EX(2, 3, 3), hz0 = HZ(2, 3, 3);

	for (int t = 0; t < 6; t++) {
		const comp_t *c = &comps[t];
		const int64_t num = *c->num;
		real_t *sb = *c->split_b, *sc = *c->split_c;
		real_t **f = c->isE ? efield[c->a] : hfield[c->a];

		for (int64_t n = 0; n < num; n++) {
			sb[n] = (real_t)urand();
			sc[n] = (real_t)urand();
		}
		real_t *b0 = (real_t *)malloc(num * sizeof(real_t));
		real_t *c0 = (real_t *)malloc(num * sizeof(real_t));
		memcpy(b0, sb, num * sizeof(real_t));
		memcpy(c0, sc, num * sizeof(real_t));

		/* E は H だけ、H は E だけを読むので、更新前の場をそのまま使える */
		c->update();

		double err = 0, scale = 0;
		int lossy = 0;
		for (int64_t n = 0; n < num; n++) {
			const pml_t *p = &(*c->list)[n];
			double db, gb, dc, gc;
			expect_split(c, p, 0, &db, &gb);
			expect_split(c, p, 1, &dc, &gc);
			const double eb = (b0[n] + db) / (1 + gb);
			const double ec = (c0[n] + dc) / (1 + gc);
			err = fmax(err, fabs(sb[n] - eb));
			err = fmax(err, fabs(sc[n] - ec));
			err = fmax(err, fabs(fval(f, p->i, p->j, p->k) - (eb + ec)));
			scale = fmax(scale, fabs(eb) + fabs(ec));
			lossy += (gb > 0) && (gc > 0);
		}
		CHECK(err < 1e-5 * scale);
		CHECK(lossy > 0);   /* 両方向に損失のある (角の) 点も検査している */
		if (err >= 1e-5 * scale) {
			printf("    %s : max error %.3e (scale %.3e)\n", c->name, err, scale);
		}

		free(b0);
		free(c0);
	}

	CHECK(EX(2, 3, 3) == ex0);
	CHECK(HZ(2, 3, 3) == hz0);
}

void test_pml(void)
{
	const double dt0 = Dt;
	const int abc0 = iABC;

	pml_setup();
	make_comps();
	test_profile();
	test_lists();
	test_update();
	pml_teardown();

	Dt = dt0;
	iABC = abc0;
}
