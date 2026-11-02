/*
test_alloc.c

配列の確保と解放 memalloc1/2/3 / memfree1/2/3 (sol/memallocfree?.c) のテスト。

確保の大きさと解放の漏れは、普通に走らせても結果に現れない (確保が少し
足りなくても隣の領域を壊すだけで "normal end" になり得る、解放漏れは終了時に
OS が回収する)。ここでは構成を変えて (Mur / PML、vector の有無、近傍界の
周波数の有無、給電・観測点) 確保し、各配列を**使われる大きさいっぱいまで
書き込み**、解放する。

  - 書き込む大きさは、その配列を使う側の範囲から決める
    (例: VPoint は観測点 + 差分用の 2 点、PML の損失分布は -L .. N+L-1 の添字)
  - CI の AddressSanitizer のステップ (ユニットテストを ASan で走らせる) では、
    確保が足りなければ heap-buffer-overflow として落ちる
  - 解放後にすべてのグローバルのポインタを NULL にする。memfree が解放し
    忘れた配列は、どこからも指されなくなるので LeakSanitizer が漏れとして報告する
    (ポインタが残っていると「到達可能」とみなされ報告されない)
  - ASan の無い通常のビルドでも、確保・書き込み・解放の一連が通ること
  - 書き込んだ後に initfield() を呼び、時間発展の状態 (E/H、Mur の 1 つ前の値、
    PML の分割成分、収束履歴、給電・観測点の時系列、近傍界 DFT の和) が
    すべて 0 に戻ること (使う側の範囲いっぱいまで)
*/

#include "ofd.h"
#include "ofd_prototype.h"
#include "ofd_test.h"

#include <stdlib.h>
#include <string.h>

#define TOUCH(p, n) do { \
		ok &= ((p) != NULL); \
		if ((p) != NULL) memset((p), 0x5a, (size_t)(n) * sizeof(*(p))); \
	} while (0)

/* p[0 .. n-1] がすべて 0 か (0x5a で埋めた後に initfield で戻ったか) */
#define ZERO(p, n) do { \
		const unsigned char *b_ = (const unsigned char *)(p); \
		const size_t nb_ = (size_t)(n) * sizeof(*(p)); \
		int z_ = ((p) != NULL); \
		for (size_t m_ = 0; z_ && (m_ < nb_); m_++) z_ = (b_[m_] == 0); \
		zero &= z_; \
		if (!z_) printf("    initfield : %s not zeroed\n", #p); \
	} while (0)

static void null_all(void)
{
	RXn = RYn = RZn = RXc = RYc = RZc = NULL;
	iEx = iEy = iEz = iHx = iHy = iHz = NULL;
	C1 = C2 = D1 = D2 = NULL;
	K1Ex = K2Ex = K1Ey = K2Ey = K1Ez = K2Ez = NULL;
	K1Hx = K2Hx = K1Hy = K2Hy = K1Hz = K2Hz = NULL;
	fMurHx = fMurHy = fMurHz = NULL;
	fPmlEx = fPmlEy = fPmlEz = fPmlHx = fPmlHy = fPmlHz = NULL;
	gPmlXn = gPmlYn = gPmlZn = gPmlXc = gPmlYc = gPmlZc = NULL;
	rPmlE = rPmlH = rPml = NULL;
	Eiter = Hiter = VFeed = IFeed = VPoint = NULL;
	Spara = NULL;
	cEdft = cHdft = cFdft = NULL;
	Ex = Ey = Ez = Hx = Hy = Hz = NULL;
	Exy = Exz = Eyz = Eyx = Ezx = Ezy = NULL;
	Hxy = Hxz = Hyz = Hyx = Hzx = Hzy = NULL;
	cEx_r = cEx_i = cEy_r = cEy_i = cEz_r = cEz_i = NULL;
	cHx_r = cHx_i = cHy_r = cHy_i = cHz_r = cHz_i = NULL;
	Xn = Yn = Zn = Xc = Yc = Zc = NULL;
	Feed = NULL; Point = NULL; Freq1 = Freq2 = NULL;
}

static void run_config(int abc, int vector, int nfreq2)
{
	/* 入力相当 (input_data が作るもの) */
	Nx = 5; Ny = 4; Nz = 6;
	iABC = abc;
	cPML.l = 3; cPML.m = 2; cPML.r0 = 1e-5;
	VECTOR = vector;
	Npx = Npy = Npz = 1;
	NMaterial = 4;
	Solver.maxiter = 30;
	NFeed = 2;
	NPoint = 1;
	NFreq1 = 3;
	NFreq2 = nfreq2;
	Xn = (double *)calloc(Nx + 1, sizeof(double)); Xc = (double *)calloc(Nx, sizeof(double));
	Yn = (double *)calloc(Ny + 1, sizeof(double)); Yc = (double *)calloc(Ny, sizeof(double));
	Zn = (double *)calloc(Nz + 1, sizeof(double)); Zc = (double *)calloc(Nz, sizeof(double));
	Feed = (feed_t *)calloc(NFeed, sizeof(feed_t));
	Point = (point_t *)calloc(NPoint + 2, sizeof(point_t));
	Freq1 = (double *)calloc(NFreq1, sizeof(double));
	/* input_data と同じく、frequency2 が無ければ確保しない (memfree1 も NFreq2 > 0 のときだけ解放) */
	Freq2 = (NFreq2 > 0) ? (double *)calloc(NFreq2, sizeof(double)) : NULL;

	/* sol_Main と同じ順 */
	CHECK(setupSize(1, 1, 1, 0) == 0);
	setupABCsize();
	memalloc1();
	memalloc2();
	memalloc3();

	int ok = 1;
	const int64_t nt = Solver.maxiter + 1;
	const int l = (abc == 1) ? cPML.l : 0;

	/* 使う側の範囲いっぱいまで書く */
	TOUCH(RXn, Nx + 1); TOUCH(RYn, Ny + 1); TOUCH(RZn, Nz + 1);
	TOUCH(RXc, Nx);     TOUCH(RYc, Ny);     TOUCH(RZc, Nz);
	TOUCH(iEx, NN); TOUCH(iEy, NN); TOUCH(iEz, NN);
	TOUCH(iHx, NN); TOUCH(iHy, NN); TOUCH(iHz, NN);
	TOUCH(C1, NMaterial); TOUCH(C2, NMaterial); TOUCH(D1, NMaterial); TOUCH(D2, NMaterial);
	if (vector) {
		TOUCH(K1Ex, NN); TOUCH(K2Ex, NN); TOUCH(K1Ey, NN); TOUCH(K2Ey, NN); TOUCH(K1Ez, NN); TOUCH(K2Ez, NN);
		TOUCH(K1Hx, NN); TOUCH(K2Hx, NN); TOUCH(K1Hy, NN); TOUCH(K2Hy, NN); TOUCH(K1Hz, NN); TOUCH(K2Hz, NN);
	}
	if (abc == 0) {
		TOUCH(fMurHx, numMurHx); TOUCH(fMurHy, numMurHy); TOUCH(fMurHz, numMurHz);
		CHECK(numMurHx > 0 && numMurHy > 0 && numMurHz > 0);
	}
	else {
		TOUCH(fPmlEx, numPmlEx); TOUCH(fPmlEy, numPmlEy); TOUCH(fPmlEz, numPmlEz);
		TOUCH(fPmlHx, numPmlHx); TOUCH(fPmlHy, numPmlHy); TOUCH(fPmlHz, numPmlHz);
		TOUCH(Exy, numPmlEx); TOUCH(Exz, numPmlEx); TOUCH(Eyz, numPmlEy); TOUCH(Eyx, numPmlEy);
		TOUCH(Ezx, numPmlEz); TOUCH(Ezy, numPmlEz);
		TOUCH(Hxy, numPmlHx); TOUCH(Hxz, numPmlHx); TOUCH(Hyz, numPmlHy); TOUCH(Hyx, numPmlHy);
		TOUCH(Hzx, numPmlHz); TOUCH(Hzy, numPmlHz);
		/* 損失分布は添字 i + L、i = -L .. N+L-1 (sol/setupPml.c, sol/pml*.c) */
		TOUCH(gPmlXn, Nx + 2 * l); TOUCH(gPmlXc, Nx + 2 * l);
		TOUCH(gPmlYn, Ny + 2 * l); TOUCH(gPmlYc, Ny + 2 * l);
		TOUCH(gPmlZn, Nz + 2 * l); TOUCH(gPmlZc, Nz + 2 * l);
		TOUCH(rPmlE, NMaterial); TOUCH(rPmlH, NMaterial); TOUCH(rPml, NMaterial);
		CHECK(numPmlEx > 0 && numPmlHz > 0);
	}
	TOUCH(Ex, NN); TOUCH(Ey, NN); TOUCH(Ez, NN); TOUCH(Hx, NN); TOUCH(Hy, NN); TOUCH(Hz, NN);
	TOUCH(Eiter, nt); TOUCH(Hiter, nt);
	TOUCH(VFeed, NFeed * nt); TOUCH(IFeed, NFeed * nt);
	TOUCH(VPoint, (NPoint + 2) * nt);          /* 観測点 + 差分用の 2 点 (sol/vpoint.c) */
	TOUCH(Spara, NPoint * NFreq1);
	if (nfreq2 > 0) {
		TOUCH(cEdft, NFreq2 * nt); TOUCH(cHdft, NFreq2 * nt); TOUCH(cFdft, NFreq2);
		const int64_t nf = (int64_t)NFreq2 * NN;
		TOUCH(cEx_r, nf); TOUCH(cEx_i, nf); TOUCH(cEy_r, nf); TOUCH(cEy_i, nf); TOUCH(cEz_r, nf); TOUCH(cEz_i, nf);
		TOUCH(cHx_r, nf); TOUCH(cHx_i, nf); TOUCH(cHy_r, nf); TOUCH(cHy_i, nf); TOUCH(cHz_r, nf); TOUCH(cHz_i, nf);
	}
	CHECK(ok);
	if (!ok) printf("    abc=%d vector=%d nfreq2=%d : NULL array\n", abc, vector, nfreq2);

	/* initfield : 時間発展の状態をすべて 0 に戻す。残っていると結果が
	   確保直後のメモリの中身に依存する (前の計算の値を引き継ぐ) */
	Niter = 12345;
	initfield();
	int zero = 1;
	ZERO(Ex, NN); ZERO(Ey, NN); ZERO(Ez, NN); ZERO(Hx, NN); ZERO(Hy, NN); ZERO(Hz, NN);
	if (abc == 0) {
		int fz = 1;
		for (int64_t n = 0; n < numMurHx; n++) fz &= (fMurHx[n].f == 0);
		for (int64_t n = 0; n < numMurHy; n++) fz &= (fMurHy[n].f == 0);
		for (int64_t n = 0; n < numMurHz; n++) fz &= (fMurHz[n].f == 0);
		CHECK(fz);
	}
	else {
		ZERO(Exy, numPmlEx); ZERO(Exz, numPmlEx); ZERO(Eyz, numPmlEy); ZERO(Eyx, numPmlEy);
		ZERO(Ezx, numPmlEz); ZERO(Ezy, numPmlEz);
		ZERO(Hxy, numPmlHx); ZERO(Hxz, numPmlHx); ZERO(Hyz, numPmlHy); ZERO(Hyx, numPmlHy);
		ZERO(Hzx, numPmlHz); ZERO(Hzy, numPmlHz);
	}
	ZERO(Eiter, nt); ZERO(Hiter, nt);
	ZERO(VFeed, NFeed * nt); ZERO(IFeed, NFeed * nt);
	ZERO(VPoint, (NPoint + 2) * nt);
	if (nfreq2 > 0) {
		const int64_t nf = (int64_t)NFreq2 * NN;
		ZERO(cEx_r, nf); ZERO(cEx_i, nf); ZERO(cEy_r, nf); ZERO(cEy_i, nf); ZERO(cEz_r, nf); ZERO(cEz_i, nf);
		ZERO(cHx_r, nf); ZERO(cHx_i, nf); ZERO(cHy_r, nf); ZERO(cHy_i, nf); ZERO(cHz_r, nf); ZERO(cHz_i, nf);
	}
	CHECK(zero);
	CHECK(Niter == 0);

	/* 解放 : memfree1 は入力由来の Xn..Zc / Feed / Point / Freq1/2 も解放する */
	memfree1();
	memfree2();
	memfree3();
	null_all();
}

void test_alloc(void)
{
	const int abc0 = iABC, vec0 = VECTOR;
	run_config(0, 0, 0);
	run_config(0, 1, 2);
	run_config(1, 0, 2);
	run_config(1, 1, 0);
	iABC = abc0;
	VECTOR = vec0;
	Nx = Ny = Nz = 0;
	NFeed = NPoint = NFreq1 = NFreq2 = 0;
	NMaterial = 0;
}
