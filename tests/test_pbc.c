/*
test_pbc.c

周期境界 (pbc = 1 1 1) のテスト。

対象: pbcx() / pbcy() / pbcz()  (sol/pbc[xyz].c)

期待値の出所 (周期性の定義そのもの):
  周期 L の方向では、境界の外側 1 セルの H は反対側の内部の値と等しい。
  H_a は自分の向き a には節点、他の 2 軸にはセル中心に乗っているので、
  他の 2 軸のセル番号について
      -1 → N-1,   N → 0      (wrap)
  と読み替えた値に一致すればよい。

  sol/solve.c は pbcx → pbcy → pbcz の順に呼ぶ。後の方向のループが前の方向の
  ゴースト行も含む範囲を回るので、辺・角のゴースト (例えば Hz(-1, -1, k)) も
  2 回の写しを経て Hz(Nx-1, Ny-1, k) に一致するはず。ここでは全方向周期にして
  辺・角を含むゴーストすべてを wrap と比べ、それ以外の値 (内部・法線成分の
  ゴースト) が書き換わっていないことも確かめる。

  単一方向の周期 (pbc = 1 0 0 など) では、その方向のゴースト面だけが
  変わること (辺の範囲には他方向のゴーストセルを含む)。
*/

#include "ofd.h"
#include "ofd_prototype.h"
#include "ofd_test.h"

#include <stdint.h>
#include <stdlib.h>
#include <string.h>

static int n_[3];

static void pbc_setup(void)
{
	n_[0] = Nx = 4;
	n_[1] = Ny = 5;
	n_[2] = Nz = 3;
	iMin = jMin = kMin = 0;
	iMax = Nx; jMax = Ny; kMax = Nz;

	Nk = 1;
	Nj = (kMax - kMin) + 2 + 1;
	Ni = ((jMax - jMin) + 2 + 1) * Nj;
	N0 = -((iMin - 1) * Ni + (jMin - 1) * Nj + (kMin - 1) * Nk);
	NN = NA(iMax + 1, jMax + 1, kMax + 1) + 1;

	Hx = (real_t *)malloc(NN * sizeof(real_t));
	Hy = (real_t *)malloc(NN * sizeof(real_t));
	Hz = (real_t *)malloc(NN * sizeof(real_t));
}

static void pbc_teardown(void)
{
	free(Hx); free(Hy); free(Hz);
	Hx = Hy = Hz = NULL;
}

/* 各要素に固有の値 (どこから写されたか分かるように) */
static void fill(void)
{
	real_t *h[3] = {Hx, Hy, Hz};
	for (int c = 0; c < 3; c++) {
		for (int64_t q = 0; q < NN; q++) {
			h[c][q] = (real_t)(1000 * (c + 1) + q);
		}
	}
}

static int wrap(int i, int n)
{
	return (i < 0) ? (n - 1) : ((i >= n) ? 0 : i);
}

/* periodic[3] の方向について、ゴーストの期待値を調べる。
   戻り値 : 不一致の個数 */
static int check(const int periodic[3], real_t *orig[3])
{
	real_t *h[3] = {Hx, Hy, Hz};
	int bad = 0;

	for (int a = 0; a < 3; a++) {
		/* 配列の全要素 (余白 1 層を含む) を走査 */
		for (int i = -1; i <= Nx + 1; i++) {
		for (int j = -1; j <= Ny + 1; j++) {
		for (int k = -1; k <= Nz + 1; k++) {
			const int idx[3] = {i, j, k};
			const int64_t q = NA(i, j, k);
			if ((q < 0) || (q >= NN)) continue;

			/* この要素が「周期の写し先」か : a 以外の周期方向でセル番号が
			   -1 か N、かつ全軸が有効範囲 (a は節点 0..N、他はセル -1..N) */
			int valid = 1, ghost = 0;
			int src[3];
			for (int b = 0; b < 3; b++) {
				const int n = n_[b];
				if (b == a) {
					valid &= (idx[b] >= 0) && (idx[b] <= n);
					src[b] = idx[b];
				}
				else {
					valid &= (idx[b] >= -1) && (idx[b] <= n);
					if (periodic[b] && ((idx[b] == -1) || (idx[b] == n))) {
						ghost = 1;
						src[b] = wrap(idx[b], n);
					}
					else {
						src[b] = idx[b];
					}
				}
			}

			/* 写し元は、周期でない方向のゴースト (Mur 等が書く) なら元の値 */
			const real_t expect = (valid && ghost)
				? orig[a][NA(src[0], src[1], src[2])]
				: orig[a][q];
			if (h[a][q] != expect) bad++;
		}
		}
		}
	}

	return bad;
}

static void run_case(int px, int py, int pz)
{
	const int periodic[3] = {px, py, pz};
	real_t *orig[3];

	fill();
	orig[0] = (real_t *)malloc(NN * sizeof(real_t));
	orig[1] = (real_t *)malloc(NN * sizeof(real_t));
	orig[2] = (real_t *)malloc(NN * sizeof(real_t));
	memcpy(orig[0], Hx, NN * sizeof(real_t));
	memcpy(orig[1], Hy, NN * sizeof(real_t));
	memcpy(orig[2], Hz, NN * sizeof(real_t));

	/* sol/solve.c と同じ順 */
	if (px) pbcx();
	if (py) pbcy();
	if (pz) pbcz();

	const int bad = check(periodic, orig);
	CHECK(bad == 0);
	if (bad) printf("    pbc = %d %d %d : %d mismatches\n", px, py, pz, bad);

	free(orig[0]); free(orig[1]); free(orig[2]);
}

void test_pbc(void)
{
	pbc_setup();

	/* 単一方向 / 2 方向 (tpa_slab.ofd などの pbc = 1 1 0) / 全方向 */
	run_case(1, 0, 0);
	run_case(0, 1, 0);
	run_case(0, 0, 1);
	run_case(1, 1, 0);
	run_case(0, 1, 1);
	run_case(1, 0, 1);
	run_case(1, 1, 1);

	/* 全方向の角 : Hz(-1,-1,k) は 2 回の写しで Hz(Nx-1, Ny-1, k) になる */
	fill();
	const real_t corner = HZ(Nx - 1, Ny - 1, 1);
	pbcx(); pbcy(); pbcz();
	CHECK(HZ(-1, -1, 1) == corner);
	CHECK(HX(2, -1, -1) == HX(2, Ny - 1, Nz - 1));
	CHECK(HY(Nx, 3, Nz) == HY(0, 3, 0));

	pbc_teardown();
}
