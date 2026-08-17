/*
test_farfield.c

sol/farfield.c (遠方界変換) と sol/nearfield_c.c (節点補間) のユニットテスト。

遠方界は等価電磁流の面積分

    Z0*N = ∮ (Z0*J) e^{jk r̂·r'} dS,   Z0*J = n̂ x (Z0*H)
    L    = ∮ M       e^{jk r̂·r'} dS,   M    = -n̂ x E
    Eθ = ffctr (Z0*Nθ + Lφ),  Eφ = ffctr (Z0*Nφ - Lθ)

で計算される。既存の検証は `data/sample/sphere_rcs_check.sh` (完全導体球の
Mie 厳密解と比較) だけで、これは端から端まで通した判定なので**どの段で
壊れたかは分からない**うえ許容幅も 0.625〜1.6 倍と広い。ここでは各段を
単体で、解析的に分かっている値と突き合わせる。

期待値の出所 (`.claude/rules/features.md` の方針どおり、実装式の写しではなく
独立な根拠を使う):

 1. farfield() … 面を微小パッチに退化させ、**E と H を直接指定**して
    教科書的なパターンと比べる。等価電磁流の向きは外積から手で求めて
    コメントに書いてある (実装の符号規約をテスト側で再利用しない —
    そうすると符号を反転させる変異が検出できなくなる)。
      n̂ = x̂, Z0*Hy = h  ->  Z0*J = x̂ x ŷh = +ẑh   (Hertz ダイポール)
      n̂ = x̂, Z0*Hz = h  ->  Z0*J = x̂ x ẑh = -ŷh
      n̂ = x̂, Ez = e     ->  M = -x̂ x ẑe = +ŷe
      n̂ = x̂, Ey = e     ->  M = -x̂ x ŷe = -ẑe
    球面座標への射影は t̂ = (cosθcosφ, cosθsinφ, -sinθ),
    p̂ = (-sinφ, cosφ, 0) から求める。位相因子は**軸から外した 1 枚の
    パッチ**で実部と虚部の両方を見る (対称配置だと cos が偶関数なので
    符号反転が観測できない)。
 2. farComponent() … 偏波の分解。直線偏波なら短軸 0・RHCP=LHCP、円偏波なら
    長軸=短軸で片方の円偏波成分が 0、さらにどの偏波でも
    |E|² = 長軸² + 短軸² = RHCP² + LHCP² が恒等的に成り立つ。
 3. NodeE_c / NodeH_c … **一次関数の場**を入れる。等間隔格子なら内点の
    2 点平均も境界の 1.5c1-0.5c2 も一次関数を厳密に再現するので、節点での
    値が解析値と一致することを見る (定数場だとどちらも通ってしまい、
    境界外挿を単純な代入に置き換える変異が検出できない)。
 4. setup_farfield() … 面の幾何 (法線・位置・面素) と 4 点平均。面素の総和が
    箱の表面積 2(ab+bc+ca) と一致すること、面に垂直な成分が 0 に落とされて
    いること、PBC 指定でその向きの面素が 0 になることを見る。
 5. farfactor() … 周波数に比例し給電電力の平方根に反比例すること、
    給電が平面波より優先されること、MatchingLoss が Pin の行を替えること。
*/

#include "ofd.h"
#include "complex.h"
#include "ofd_test.h"

#include <stdlib.h>
#include <string.h>

void   farfield(int ifreq, double theta, double phi, double ffctr, d_complex_t *etheta, d_complex_t *ephi);
void   farComponent(d_complex_t etheta, d_complex_t ephi, double e[]);
double farfactor(int ifreq);
void   alloc_farfield(void);
void   setup_farfield(void);
void   NodeE_c(int ifreq, int i, int j, int k, d_complex_t *cex, d_complex_t *cey, d_complex_t *cez);
void   NodeH_c(int ifreq, int i, int j, int k, d_complex_t *chx, d_complex_t *chy, d_complex_t *chz);

#define TOL (1e-12)

/* ---------------------------------------------------------------- */
/* 1. farfield() : 微小ダイポールの放射パターン                      */
/* ---------------------------------------------------------------- */

/* 微小パッチを 1 枚置く。指定するのは**物理量である E と H** で、
   等価電磁流は farfield() の中で外積から作られる。 */
static void set_patch(int n, double nx, double ny, double nz,
                      double x, double y, double z, double ds,
                      double ex, double ey, double ez,
                      double hx, double hy, double hz)
{
	Surface[n].nx = nx;
	Surface[n].ny = ny;
	Surface[n].nz = nz;
	Surface[n].x = x;
	Surface[n].y = y;
	Surface[n].z = z;
	Surface[n].ds = ds;

	SurfaceEx[0][n] = d_complex(ex, 0);
	SurfaceEy[0][n] = d_complex(ey, 0);
	SurfaceEz[0][n] = d_complex(ez, 0);
	SurfaceHx[0][n] = d_complex(hx, 0);
	SurfaceHy[0][n] = d_complex(hy, 0);
	SurfaceHz[0][n] = d_complex(hz, 0);
}

static void alloc_patches(int n)
{
	NFreq2 = 1;
	Freq2 = (double *)malloc(sizeof(double));
	Freq2[0] = C;            /* k = 2π/λ で λ = 1 [m] になる */
	NSurface = n;
	Surface = (surface_t *)calloc(n, sizeof(surface_t));
	SurfaceEx = (d_complex_t **)malloc(sizeof(d_complex_t *));
	SurfaceEy = (d_complex_t **)malloc(sizeof(d_complex_t *));
	SurfaceEz = (d_complex_t **)malloc(sizeof(d_complex_t *));
	SurfaceHx = (d_complex_t **)malloc(sizeof(d_complex_t *));
	SurfaceHy = (d_complex_t **)malloc(sizeof(d_complex_t *));
	SurfaceHz = (d_complex_t **)malloc(sizeof(d_complex_t *));
	SurfaceEx[0] = (d_complex_t *)calloc(n, sizeof(d_complex_t));
	SurfaceEy[0] = (d_complex_t *)calloc(n, sizeof(d_complex_t));
	SurfaceEz[0] = (d_complex_t *)calloc(n, sizeof(d_complex_t));
	SurfaceHx[0] = (d_complex_t *)calloc(n, sizeof(d_complex_t));
	SurfaceHy[0] = (d_complex_t *)calloc(n, sizeof(d_complex_t));
	SurfaceHz[0] = (d_complex_t *)calloc(n, sizeof(d_complex_t));
}

static void free_patches(void)
{
	free(SurfaceEx[0]); free(SurfaceEy[0]); free(SurfaceEz[0]);
	free(SurfaceHx[0]); free(SurfaceHy[0]); free(SurfaceHz[0]);
	free(SurfaceEx); free(SurfaceEy); free(SurfaceEz);
	free(SurfaceHx); free(SurfaceHy); free(SurfaceHz);
	free(Surface);
	free(Freq2);
}

static void test_dipole_pattern(void)
{
	alloc_patches(2);
	NSurface = 1;

	d_complex_t et, ep;
	const double th = 50.0, ph = 35.0;
	const double ct = cos(th * DTOR), st = sin(th * DTOR);
	const double cp = cos(ph * DTOR), sp = sin(ph * DTOR);

	/* (a) n̂ = x̂, Z0*Hy = 1  ->  Z0*J = +ẑ (Hertz ダイポール)
	   Nθ = t̂·ẑ = -sinθ, Nφ = p̂·ẑ = 0  ->  Eθ = -sinθ, Eφ = 0 */
	set_patch(0, 1, 0, 0, 0, 0, 0, 1.0,  0, 0, 0,  0, 1, 0);
	for (int it = 0; it <= 6; it++) {
		const double t = it * 30.0;
		farfield(0, t, ph, 1.0, &et, &ep);
		CHECK_NEAR(et.r, -sin(t * DTOR), TOL);
		CHECK_NEAR(et.i, 0.0, TOL);
		CHECK_NEAR(ep.r, 0.0, TOL);
		CHECK_NEAR(ep.i, 0.0, TOL);
	}
	/* パターンの要点 : 軸上 (θ=0) で 0、赤道 (θ=90) で最大、φ に依らない */
	{
		double e0[7], e90[7], e90b[7];
		farfield(0, 0.0, 0.0, 1.0, &et, &ep);    farComponent(et, ep, e0);
		farfield(0, 90.0, 0.0, 1.0, &et, &ep);   farComponent(et, ep, e90);
		farfield(0, 90.0, 123.0, 1.0, &et, &ep); farComponent(et, ep, e90b);
		CHECK_NEAR(e0[0], 0.0, TOL);
		CHECK_NEAR(e90[0], 1.0, TOL);
		CHECK_NEAR(e90[0], e90b[0], TOL);
	}

	/* (b) n̂ = x̂, Z0*Hz = 1  ->  Z0*J = x̂ x ẑ = -ŷ
	   Nθ = t̂·(-ŷ) = -cosθ sinφ, Nφ = p̂·(-ŷ) = -cosφ */
	set_patch(0, 1, 0, 0, 0, 0, 0, 1.0,  0, 0, 0,  0, 0, 1);
	farfield(0, th, ph, 1.0, &et, &ep);
	CHECK_NEAR(et.r, -ct * sp, TOL);
	CHECK_NEAR(ep.r, -cp, TOL);
	CHECK_NEAR(et.i, 0.0, TOL);

	/* (c) n̂ = x̂, Ez = 1  ->  M = -x̂ x ẑ = +ŷ
	   Lθ = t̂·ŷ = cosθ sinφ, Lφ = p̂·ŷ = cosφ
	   Eθ = +Lφ = cosφ,  Eφ = -Lθ = -cosθ sinφ
	   (Lφ が Eθ に足される経路はここだけで踏まれる) */
	set_patch(0, 1, 0, 0, 0, 0, 0, 1.0,  0, 0, 1,  0, 0, 0);
	farfield(0, th, ph, 1.0, &et, &ep);
	CHECK_NEAR(et.r, +cp, TOL);
	CHECK_NEAR(ep.r, -ct * sp, TOL);

	/* (d) n̂ = x̂, Ey = 1  ->  M = -x̂ x ŷ = -ẑ (磁気ダイポール)
	   Lθ = t̂·(-ẑ) = +sinθ, Lφ = 0  ->  Eθ = 0, Eφ = -Lθ = -sinθ
	   電気ダイポール (a) と双対な φ 偏波の sinθ パターン */
	set_patch(0, 1, 0, 0, 0, 0, 0, 1.0,  0, 1, 0,  0, 0, 0);
	for (int it = 1; it <= 5; it++) {
		const double t = it * 30.0;
		farfield(0, t, ph, 1.0, &et, &ep);
		CHECK_NEAR(et.r, 0.0, TOL);
		CHECK_NEAR(ep.r, -sin(t * DTOR), TOL);
		CHECK_NEAR(ep.i, 0.0, TOL);
	}

	/* (e) Huygens 源 : +x 方向に進む平面波の等価流 (Ey = 1, Z0*Hz = 1)。
	   J と M が同時に効く。(b) と (d) の重ね合わせなので
	   Eθ = -cosθ sinφ + 0 = -cosθ sinφ,  Eφ = -cosφ - sinθ
	   +x 方向 (θ=90, φ=0) では Eφ = -1-1 = -2 (前方に集まる)、
	   -x 方向 (θ=90, φ=180) では Eφ = +1-1 = 0 (後方に放射しない) */
	set_patch(0, 1, 0, 0, 0, 0, 0, 1.0,  0, 1, 0,  0, 0, 1);
	farfield(0, th, ph, 1.0, &et, &ep);
	CHECK_NEAR(et.r, -ct * sp, TOL);
	CHECK_NEAR(ep.r, -cp - st, TOL);
	{
		double ef[7], eb[7];
		farfield(0, 90.0, 0.0, 1.0, &et, &ep);   farComponent(et, ep, ef);
		CHECK_NEAR(ep.r, -2.0, TOL);
		farfield(0, 90.0, 180.0, 1.0, &et, &ep); farComponent(et, ep, eb);
		CHECK_NEAR(eb[0], 0.0, TOL);             /* 後方の打ち消し */
		CHECK_NEAR(ef[0], 2.0, TOL);
	}

	/* (e2) 外積の 3 成分すべてを踏む。法線を x̂ に固定すると
	   n̂ x H の x 成分が常に 0 になり、その成分だけを壊す変異が
	   検出できないので、法線を 3 軸すべてに振って総当たりする。
	   期待値の等価流は外積を手で計算した値 (下の表の jx..mz):
	     Z0*J = n̂ x (Z0*H) :  x̂xŷ=+ẑ  x̂xẑ=-ŷ  ŷxẑ=+x̂  ŷxx̂=-ẑ  ẑxx̂=+ŷ  ẑxŷ=-x̂
	     M    = -n̂ x E     : 上の符号を反転
	   遠方界は Eθ = t̂·J + p̂·M,  Eφ = p̂·J - t̂·M。 */
	{
		/* nx,ny,nz, ex,ey,ez, hx,hy,hz,  jx,jy,jz, mx,my,mz */
		static const double tbl[12][15] = {
			/* n=x̂ に H を当てる */
			{1,0,0,  0,0,0,  0,1,0,   0,0, 1,  0,0,0},   /* x̂xŷ = +ẑ */
			{1,0,0,  0,0,0,  0,0,1,   0,-1,0,  0,0,0},   /* x̂xẑ = -ŷ */
			/* n=ŷ に H を当てる (J の x 成分と z 成分) */
			{0,1,0,  0,0,0,  0,0,1,   1,0,0,   0,0,0},   /* ŷxẑ = +x̂ */
			{0,1,0,  0,0,0,  1,0,0,   0,0,-1,  0,0,0},   /* ŷxx̂ = -ẑ */
			/* n=ẑ に H を当てる */
			{0,0,1,  0,0,0,  1,0,0,   0,1,0,   0,0,0},   /* ẑxx̂ = +ŷ */
			{0,0,1,  0,0,0,  0,1,0,   -1,0,0,  0,0,0},   /* ẑxŷ = -x̂ */
			/* n=x̂ に E を当てる (M は符号反転) */
			{1,0,0,  0,1,0,  0,0,0,   0,0,0,   0,0,-1},
			{1,0,0,  0,0,1,  0,0,0,   0,0,0,   0,1,0},
			/* n=ŷ に E を当てる (M の x 成分と z 成分) */
			{0,1,0,  0,0,1,  0,0,0,   0,0,0,   -1,0,0},
			{0,1,0,  1,0,0,  0,0,0,   0,0,0,   0,0,1},
			/* n=ẑ に E を当てる */
			{0,0,1,  1,0,0,  0,0,0,   0,0,0,   0,-1,0},
			{0,0,1,  0,1,0,  0,0,0,   0,0,0,   1,0,0},
		};
		const double t1[3] = {ct * cp, ct * sp, -st};
		const double p1[3] = {-sp, cp, 0};
		NSurface = 1;
		for (int m = 0; m < 12; m++) {
			const double *g = tbl[m];
			set_patch(0, g[0], g[1], g[2], 0, 0, 0, 1.0,
			          g[3], g[4], g[5], g[6], g[7], g[8]);
			const double jx = g[9],  jy = g[10], jz = g[11];
			const double mx = g[12], my = g[13], mz = g[14];
			const double nt = (t1[0] * jx) + (t1[1] * jy) + (t1[2] * jz);
			const double np = (p1[0] * jx) + (p1[1] * jy) + (p1[2] * jz);
			const double lt = (t1[0] * mx) + (t1[1] * my) + (t1[2] * mz);
			const double lp = (p1[0] * mx) + (p1[1] * my) + (p1[2] * mz);
			farfield(0, th, ph, 1.0, &et, &ep);
			CHECK_NEAR(et.r, nt + lp, TOL);
			CHECK_NEAR(ep.r, np - lt, TOL);
		}
	}

	/* (f) 位相因子 exp(+jk r̂·r')。軸から外した 1 枚のパッチで
	   実部と虚部の両方を見る (r̂ の 3 成分すべてが効く位置・方向を選ぶ) */
	{
		const double x0 = 0.1, y0 = 0.2, z0 = 0.3;   /* λ = 1 */
		const double kw = 2 * PI;
		const double rr = (st * cp * x0) + (st * sp * y0) + (ct * z0);
		set_patch(0, 1, 0, 0, x0, y0, z0, 1.0,  0, 0, 0,  0, 1, 0);
		farfield(0, th, ph, 1.0, &et, &ep);
		CHECK_NEAR(et.r, -st * cos(kw * rr), TOL);
		CHECK_NEAR(et.i, -st * sin(kw * rr), TOL);   /* 虚部の符号 = 位相の符号 */
		CHECK(fabs(et.i) > 0.1);                     /* 実際に位相が乗っていること */
		/* 面素は積分の重み (振幅に比例する) */
		set_patch(0, 1, 0, 0, x0, y0, z0, 2.5,  0, 0, 0,  0, 1, 0);
		farfield(0, th, ph, 1.0, &et, &ep);
		CHECK_NEAR(et.r, 2.5 * (-st * cos(kw * rr)), TOL);
		CHECK_NEAR(et.i, 2.5 * (-st * sin(kw * rr)), TOL);
	}

	/* (g) アレイファクタ : ẑ 向き電流を x = ±λ/4 に 2 枚 (間隔 λ/2) 置くと
	   x 方向で逆位相に打ち消し、y 方向では同位相で 2 倍になる */
	NSurface = 2;
	set_patch(0, 1, 0, 0, -0.25, 0, 0, 1.0,  0, 0, 0,  0, 1, 0);
	set_patch(1, 1, 0, 0, +0.25, 0, 0, 1.0,  0, 0, 0,  0, 1, 0);
	{
		double e[7];
		farfield(0, 90.0, 0.0, 1.0, &et, &ep);   farComponent(et, ep, e);
		CHECK_NEAR(e[0], 0.0, TOL);
		farfield(0, 90.0, 90.0, 1.0, &et, &ep);  farComponent(et, ep, e);
		CHECK_NEAR(e[0], 2.0, TOL);
		const double kw = 2 * PI;
		farfield(0, 90.0, 60.0, 1.0, &et, &ep);
		CHECK_NEAR(et.r, -1.0 * 2 * cos(kw * 0.25 * cos(60.0 * DTOR)), TOL);
	}

	/* (h) ffctr は全体に掛かる定数 (線形性) */
	{
		d_complex_t a_t, a_p, b_t, b_p;
		farfield(0, 55.0, 20.0, 1.0, &a_t, &a_p);
		farfield(0, 55.0, 20.0, 3.5, &b_t, &b_p);
		CHECK_NEAR(b_t.r, 3.5 * a_t.r, TOL);
		CHECK_NEAR(b_t.i, 3.5 * a_t.i, TOL);
		CHECK_NEAR(b_p.r, 3.5 * a_p.r, TOL);
	}

	/* (i) 面素 ds = 0 の面は寄与しないこと (PBC で使う経路) */
	{
		double e1[7], e2[7];
		Surface[1].ds = 0;
		farfield(0, 90.0, 90.0, 1.0, &et, &ep);  farComponent(et, ep, e1);
		NSurface = 1;
		farfield(0, 90.0, 90.0, 1.0, &et, &ep);  farComponent(et, ep, e2);
		CHECK_NEAR(e1[0], e2[0], TOL);
	}

	free_patches();
}

/* ---------------------------------------------------------------- */
/* 2. farComponent() : 偏波の分解                                    */
/* ---------------------------------------------------------------- */
static void test_far_component(void)
{
	double e[7];

	/* (a) θ 向きの直線偏波 : 長軸 = |E|、短軸 = 0、円偏波成分は等しい */
	farComponent(d_complex(2, 0), d_complex(0, 0), e);
	CHECK_NEAR(e[0], 2.0, TOL);            /* |E| */
	CHECK_NEAR(e[1], 2.0, TOL);            /* |Eθ| */
	CHECK_NEAR(e[2], 0.0, TOL);            /* |Eφ| */
	CHECK_NEAR(e[3], 2.0, TOL);            /* 長軸 */
	CHECK_NEAR(e[4], 0.0, TOL);            /* 短軸 */
	CHECK_NEAR(e[5], 2.0 / sqrt(2.0), TOL);
	CHECK_NEAR(e[6], 2.0 / sqrt(2.0), TOL);

	/* (b) 45° 傾いた直線偏波 : 成分は等しく、短軸はやはり 0 */
	farComponent(d_complex(1, 0), d_complex(1, 0), e);
	CHECK_NEAR(e[0], sqrt(2.0), TOL);
	CHECK_NEAR(e[3], sqrt(2.0), TOL);
	CHECK_NEAR(e[4], 0.0, TOL);

	/* (c) 円偏波 (Eφ = +j Eθ) : 長軸 = 短軸、片方の円偏波成分が 0 */
	farComponent(d_complex(1, 0), d_complex(0, 1), e);
	CHECK_NEAR(e[0], sqrt(2.0), TOL);
	CHECK_NEAR(e[3], 1.0, TOL);
	CHECK_NEAR(e[4], 1.0, TOL);
	CHECK_NEAR(e[5], 0.0, TOL);
	CHECK_NEAR(e[6], sqrt(2.0), TOL);

	/* 逆旋回では 2 つが入れ替わること */
	farComponent(d_complex(1, 0), d_complex(0, -1), e);
	CHECK_NEAR(e[5], sqrt(2.0), TOL);
	CHECK_NEAR(e[6], 0.0, TOL);

	/* (d) 楕円偏波 : Eθ = a, Eφ = j b は長軸 max(a,b)・短軸 min(a,b) */
	farComponent(d_complex(3, 0), d_complex(0, 1), e);
	CHECK_NEAR(e[3], 3.0, TOL);
	CHECK_NEAR(e[4], 1.0, TOL);

	/* (e) どの偏波でも成り立つ恒等式 : |E|² = 長軸² + 短軸² = RHCP² + LHCP² */
	{
		const d_complex_t ts[5][2] = {
			{{1, 0},     {0, 0}},
			{{0.3, 0.7}, {-0.2, 0.5}},
			{{1, 0},     {0, 1}},
			{{2, -1},    {0.5, 0.5}},
			{{0, 0},     {1.5, -0.25}},
		};
		for (int i = 0; i < 5; i++) {
			farComponent(ts[i][0], ts[i][1], e);
			CHECK_NEAR(e[0] * e[0], (e[3] * e[3]) + (e[4] * e[4]), 1e-12);
			CHECK_NEAR(e[0] * e[0], (e[5] * e[5]) + (e[6] * e[6]), 1e-12);
			CHECK_NEAR(e[0] * e[0], (e[1] * e[1]) + (e[2] * e[2]), 1e-12);
			CHECK(e[3] >= e[4] - 1e-15);   /* 長軸 >= 短軸 */
		}
	}
}

/* ---------------------------------------------------------------- */
/* 3./4. 節点補間と面の構成                                          */
/* ---------------------------------------------------------------- */
#define NG (4)   /* Nx = Ny = Nz */

/* 一様格子と DFT 配列を用意する。dnu != 0 なら不等間隔にする */
static double gxn[NG + 1], gyn[NG + 1], gzn[NG + 1];
static float *dftbuf[12];

static void alloc_grid(int nonuniform)
{
	const double ux[NG + 1] = {0.00, 0.25, 0.50, 0.75, 1.00};
	const double nx[NG + 1] = {0.00, 0.10, 0.30, 0.60, 1.00};
	const double ny[NG + 1] = {0.00, 0.20, 0.50, 0.70, 0.80};
	const double nz[NG + 1] = {0.00, 0.40, 0.50, 0.90, 1.20};

	Nx = Ny = Nz = NG;
	iMin = jMin = kMin = 0;
	iMax = jMax = kMax = NG;
	NFreq2 = 1;
	Freq2 = (double *)malloc(sizeof(double));
	Freq2[0] = C;
	PBCx = PBCy = PBCz = 0;

	Xn = (double *)malloc((NG + 1) * sizeof(double));
	Yn = (double *)malloc((NG + 1) * sizeof(double));
	Zn = (double *)malloc((NG + 1) * sizeof(double));
	Xc = (double *)malloc(NG * sizeof(double));
	Yc = (double *)malloc(NG * sizeof(double));
	Zc = (double *)malloc(NG * sizeof(double));
	for (int i = 0; i <= NG; i++) {
		gxn[i] = nonuniform ? nx[i] : ux[i];
		gyn[i] = nonuniform ? ny[i] : ux[i];
		gzn[i] = nonuniform ? nz[i] : ux[i];
		Xn[i] = gxn[i]; Yn[i] = gyn[i]; Zn[i] = gzn[i];
	}
	for (int i = 0; i < NG; i++) {
		Xc[i] = (gxn[i] + gxn[i + 1]) / 2;
		Yc[i] = (gyn[i] + gyn[i + 1]) / 2;
		Zc[i] = (gzn[i] + gzn[i + 1]) / 2;
	}

	/* 配列レイアウト (Mur 相当の余白 1 層。setupSize() と同じ式) */
	Nk = 1;
	Nj = NG + 2 + 1;
	Ni = (NG + 2 + 1) * Nj;
	N0 = -((-1) * Ni + (-1) * Nj + (-1) * Nk);
	NN = NA(NG + 1, NG + 1, NG + 1) + 1;

	for (int b = 0; b < 12; b++) {
		dftbuf[b] = (float *)calloc((size_t)NN, sizeof(float));
	}
	cEx_r = dftbuf[0]; cEx_i = dftbuf[1];
	cEy_r = dftbuf[2]; cEy_i = dftbuf[3];
	cEz_r = dftbuf[4]; cEz_i = dftbuf[5];
	cHx_r = dftbuf[6]; cHx_i = dftbuf[7];
	cHy_r = dftbuf[8]; cHy_i = dftbuf[9];
	cHz_r = dftbuf[10]; cHz_i = dftbuf[11];
}

static void free_grid(void)
{
	for (int b = 0; b < 12; b++) free(dftbuf[b]);
	free(Xn); free(Yn); free(Zn);
	free(Xc); free(Yc); free(Zc);
	free(Freq2);
}

/* 3. NodeE_c / NodeH_c : 一次関数の場を厳密に再現すること。
      等間隔格子では
        内点   : (c[i-1] + c[i]) / 2      -> 節点の値 (一次関数なら厳密)
        境界   : 1.5 c[0] - 0.5 c[1]      -> 節点への一次外挿 (同じく厳密)
      なので、期待値は「その節点における一次関数の値」で書ける。
      定数場だとどちらの式も通ってしまい、境界外挿を単純な代入に
      置き換える変異が検出できない。 */
static void test_node_interp(void)
{
	alloc_grid(0);   /* 等間隔 (一次外挿が厳密になる条件) */
	const double d = 0.25;

	/* Ex は (Xc[i], Yn[j], Zn[k]) にある。x 方向に一次な場を入れる。
	   添字は -1..NG+1 まで回るので、セル中心の位置を d*(i+0.5) で
	   延長して埋める (場の一次性を格子の外まで保つ) */
	const double a0 = 1.0, a1 = 2.0;      /* f(x) = a0 + a1 x */
	const double b0 = 3.0, b1 = 4.0, b2 = 5.0;  /* g(y,z) = b0 + b1 y + b2 z */
	for (int i = -1; i <= NG + 1; i++) {
	for (int j = -1; j <= NG + 1; j++) {
	for (int k = -1; k <= NG + 1; k++) {
		const double xc = d * (i + 0.5);
		const double yc = d * (j + 0.5);
		const double zc = d * (k + 0.5);
		cEx_r[NA(i, j, k)] = (float)(a0 + (a1 * xc));
		/* Hx は (Xn[i], Yc[j], Zc[k]) にあり、y と z の 4 点平均で節点へ運ぶ */
		cHx_r[NA(i, j, k)] = (float)(b0 + (b1 * yc) + (b2 * zc));
	}
	}
	}

	d_complex_t cex, cey, cez, chx, chy, chz;
	int bad_e = 0, bad_h = 0;
	for (int i = 0; i <= NG; i++) {
	for (int j = 1; j <= NG - 1; j++) {
	for (int k = 1; k <= NG - 1; k++) {
		NodeE_c(0, i, j, k, &cex, &cey, &cez);
		/* 節点 Xn[i] = d*i における一次関数の値 (境界 i=0, i=NG も含む) */
		if (fabs(cex.r - (a0 + (a1 * d * i))) > 1e-5) bad_e++;
		NodeH_c(0, i, j, k, &chx, &chy, &chz);
		if (fabs(chx.r - (b0 + (b1 * d * j) + (b2 * d * k))) > 1e-5) bad_h++;
	}
	}
	}
	CHECK(bad_e == 0);   /* 内点の 2 点平均と境界の一次外挿 (i=0, i=Nx) */
	CHECK(bad_h == 0);   /* H の 4 点平均 */

	/* 境界の外挿が実際に効いていること : i=0 の節点値は
	   最初のセル中心の値とは違う (半セル分だけ外側) */
	NodeE_c(0, 0, 2, 2, &cex, &cey, &cez);
	CHECK_NEAR(cex.r, a0, 1e-5);                       /* Xn[0] = 0 */
	CHECK(fabs(cex.r - (a0 + (a1 * d * 0.5))) > 0.1);  /* Xc[0] の値ではない */
	NodeE_c(0, NG, 2, 2, &cex, &cey, &cez);
	CHECK_NEAR(cex.r, a0 + (a1 * d * NG), 1e-5);       /* Xn[Nx] = 1 */

	free_grid();
}

/* 4. setup_farfield() : 面の幾何と 4 点平均 */
static void test_setup_surface(void)
{
	alloc_grid(1);   /* 不等間隔 (面素が節点間隔の積であることを見る) */

	/* DFT 配列を定数で埋める。節点補間は定数を厳密に保つので、
	   面の 4 点平均もその定数になる (下の期待値の根拠) */
	const float ce = 2.5f, ch = -1.25f;
	for (int64_t n = 0; n < NN; n++) {
		cEx_r[n] = cEy_r[n] = cEz_r[n] = ce;
		cHx_r[n] = cHy_r[n] = cHz_r[n] = ch;
	}

	alloc_farfield();

	/* 面の総枚数 : 6 面それぞれ Nx*Ny 等の 2 倍 */
	CHECK(NSurface == 2 * ((NG * NG) + (NG * NG) + (NG * NG)));

	setup_farfield();

	const double a = gxn[NG] - gxn[0];
	const double b = gyn[NG] - gyn[0];
	const double c = gzn[NG] - gzn[0];

	/* (a) 面素の総和が箱の表面積 2(ab+bc+ca) と一致すること */
	{
		double sum = 0;
		for (int64_t n = 0; n < NSurface; n++) sum += Surface[n].ds;
		CHECK_NEAR(sum, 2 * ((a * b) + (b * c) + (c * a)), 1e-12);
	}

	/* (b) 法線は外向きの単位ベクトルで、軸に平行な 6 種類しかないこと。
	       各向きがちょうど Nx*Ny 枚ずつあること */
	{
		int cnt[6] = {0, 0, 0, 0, 0, 0};
		int bad = 0;
		for (int64_t n = 0; n < NSurface; n++) {
			const double nx = Surface[n].nx, ny = Surface[n].ny, nz = Surface[n].nz;
			if (fabs((nx * nx) + (ny * ny) + (nz * nz) - 1) > 1e-15) bad++;
			if      (nx == -1) cnt[0]++;
			else if (nx == +1) cnt[1]++;
			else if (ny == -1) cnt[2]++;
			else if (ny == +1) cnt[3]++;
			else if (nz == -1) cnt[4]++;
			else if (nz == +1) cnt[5]++;
			else bad++;
		}
		CHECK(bad == 0);
		for (int s = 0; s < 6; s++) CHECK(cnt[s] == NG * NG);
	}

	/* (c) 面の位置 : 法線方向は境界の節点、面内はセル中心に乗ること。
	       法線方向に節点でなくセル中心を使う取り違えを捕まえる。
	       面内の座標がちょうどセル中心の集合と一致することも数えて見る */
	{
		int bad = 0, nc = 0;
		for (int64_t n = 0; n < NSurface; n++) {
			const surface_t *p = &Surface[n];
			if (p->nx != 0) {
				if (fabs(p->x - ((p->nx < 0) ? gxn[0] : gxn[NG])) > 1e-15) bad++;
				for (int j = 0; j < NG; j++) if (fabs(p->y - Yc[j]) < 1e-15) nc++;
			}
			else if (p->ny != 0) {
				if (fabs(p->y - ((p->ny < 0) ? gyn[0] : gyn[NG])) > 1e-15) bad++;
				for (int k = 0; k < NG; k++) if (fabs(p->z - Zc[k]) < 1e-15) nc++;
			}
			else {
				if (fabs(p->z - ((p->nz < 0) ? gzn[0] : gzn[NG])) > 1e-15) bad++;
				for (int i = 0; i < NG; i++) if (fabs(p->x - Xc[i]) < 1e-15) nc++;
			}
		}
		CHECK(bad == 0);
		CHECK(nc == NSurface);   /* 面内座標はすべてセル中心 */
	}

	/* (d) 4 点平均 : 一様な場なら面の値もその定数。
	       面に垂直な成分 (X 面の Ex/Hx など) は 0 に落とされること */
	{
		int bad = 0;
		for (int64_t n = 0; n < NSurface; n++) {
			const surface_t *p = &Surface[n];
			const double ex = SurfaceEx[0][n].r, ey = SurfaceEy[0][n].r, ez = SurfaceEz[0][n].r;
			const double hx = SurfaceHx[0][n].r, hy = SurfaceHy[0][n].r, hz = SurfaceHz[0][n].r;
			if (p->nx != 0) {
				if (ex != 0.0 || hx != 0.0) bad++;
				if (fabs(ey - ce) > 1e-6 || fabs(ez - ce) > 1e-6) bad++;
				if (fabs(hy - ch) > 1e-6 || fabs(hz - ch) > 1e-6) bad++;
			}
			else if (p->ny != 0) {
				if (ey != 0.0 || hy != 0.0) bad++;
				if (fabs(ez - ce) > 1e-6 || fabs(ex - ce) > 1e-6) bad++;
				if (fabs(hz - ch) > 1e-6 || fabs(hx - ch) > 1e-6) bad++;
			}
			else {
				if (ez != 0.0 || hz != 0.0) bad++;
				if (fabs(ex - ce) > 1e-6 || fabs(ey - ce) > 1e-6) bad++;
				if (fabs(hx - ch) > 1e-6 || fabs(hy - ch) > 1e-6) bad++;
			}
			if (SurfaceEy[0][n].i != 0.0 || SurfaceHz[0][n].i != 0.0) bad++;
		}
		CHECK(bad == 0);
	}

	/* (e) PBC 指定でその向きの面だけ面素 0 になること。
	       消えるのは X 面 (面積 b*c) の 2 枚なので、残るのは Y 面と Z 面 */
	{
		PBCx = 1;
		setup_farfield();
		double sx = 0, sother = 0;
		for (int64_t n = 0; n < NSurface; n++) {
			if (Surface[n].nx != 0) sx += Surface[n].ds;
			else                    sother += Surface[n].ds;
		}
		CHECK_NEAR(sx, 0.0, 1e-15);
		CHECK_NEAR(sother, 2 * ((c * a) + (a * b)), 1e-12);
		PBCx = 0;
	}

	for (int ifreq = 0; ifreq < NFreq2; ifreq++) {
		free(SurfaceEx[ifreq]); free(SurfaceEy[ifreq]); free(SurfaceEz[ifreq]);
		free(SurfaceHx[ifreq]); free(SurfaceHy[ifreq]); free(SurfaceHz[ifreq]);
	}
	free(SurfaceEx); free(SurfaceEy); free(SurfaceEz);
	free(SurfaceHx); free(SurfaceHy); free(SurfaceHz);
	free(Surface);
	free_grid();
}

/* ---------------------------------------------------------------- */
/* 5. farfactor() : 正規化係数                                       */
/* ---------------------------------------------------------------- */
static void test_far_factor(void)
{
	NFreq2 = 2;
	Freq2 = (double *)malloc(2 * sizeof(double));
	Freq2[0] = 1e9;
	Freq2[1] = 2e9;

	/* (a) 給電も平面波も無ければ 0 */
	NFeed = 0;
	IPlanewave = 0;
	CHECK_NEAR(farfactor(0), 0.0, 0.0);

	/* (b) 平面波 : 波数に比例する (ffctr = k / sqrt(4π)) */
	IPlanewave = 1;
	{
		const double f0 = farfactor(0);
		const double f1 = farfactor(1);
		CHECK(f0 > 0);
		CHECK_NEAR(f1, 2 * f0, 1e-12);                          /* 周波数 2 倍 */
		CHECK_NEAR(f0, (2 * PI * 1e9 / C) / sqrt(4 * PI), 1e-12);
	}

	/* (c) 給電 : 電力の平方根に反比例し、給電が平面波より優先される */
	NFeed = 1;
	Pin[0] = (double *)malloc(2 * sizeof(double));
	Pin[1] = (double *)malloc(2 * sizeof(double));
	MatchingLoss = 0;
	Pin[0][0] = 1.0;  Pin[0][1] = 4.0;
	Pin[1][0] = 0.25; Pin[1][1] = 1.0;
	{
		const double g0 = farfactor(0);
		CHECK(g0 > 0);
		/* 平面波の値とは違う = 給電側が優先されている */
		CHECK(fabs(g0 - (2 * PI * 1e9 / C) / sqrt(4 * PI)) > 1e-9);
		/* 電力 4 倍・周波数 2 倍 -> k は 2 倍、1/sqrt(P) は 1/2 倍 -> 等倍 */
		CHECK_NEAR(farfactor(1), g0, 1e-12);
		/* MatchingLoss は Pin の行を切り替える。
		   Pin[1] は Pin[0] の 1/4 なので係数は 2 倍になる */
		MatchingLoss = 1;
		CHECK_NEAR(farfactor(0), 2 * g0, 1e-12);
		MatchingLoss = 0;

		/* (d) 複数給電では電力が足し合わされること (2 給電で係数は 1/√2 倍) */
		free(Pin[0]); free(Pin[1]);
		NFeed = 2;
		Pin[0] = (double *)malloc(4 * sizeof(double));
		Pin[1] = (double *)malloc(4 * sizeof(double));
		for (int i = 0; i < 4; i++) { Pin[0][i] = 1.0; Pin[1][i] = 1.0; }
		CHECK_NEAR(farfactor(0), g0 / sqrt(2.0), 1e-12);
	}

	free(Pin[0]); free(Pin[1]);
	Pin[0] = Pin[1] = NULL;
	free(Freq2);
	NFeed = 0;
	IPlanewave = 0;
	MatchingLoss = 0;
}

void test_farfield(void)
{
	test_dipole_pattern();
	test_far_component();
	test_node_interp();
	test_setup_surface();
	test_far_factor();
}
