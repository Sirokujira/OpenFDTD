/*
test_solve.c

時間ループ solve() (sol/solve.c) のテスト。

ループの中で呼ぶ各処理 (更新・境界・給電・DFT・観測) はそれぞれのテストが
あり、物理の正しさは検証スクリプト (TPA・Mie・PML・熱・利得・定在波) が
解析解と比べている。ここではループそのものの約束を確かめる:

  - 平均値 (Eiter / Hiter) は itime % nout == 0 のステップと最終ステップ
    (itime == maxiter) で記録する (Niter 個)
  - 収束判定は記録のときだけ行い、E と H の両方がそれまでの最大の converg 倍を
    下回った最初の記録で止まる。止まったら "--- converged ---"、最後まで
    回ったら "--- max steps ---" をログに書く
  - 実行したステップ数 Ntime : 収束で止まったら itime + 1、最後まで回ったら
    maxiter + 1 (給電・観測点の時系列の有効長で、DFT の和の長さ)
  - 記録ごとに HDF5 の瞬時値スナップショットを 1 枚出す (hdf5 の interval 既定)
  - 給電の電圧・電流の時系列が実行したステップぶん埋まる

問題は dipole.ofd を小さくしたもの (Z 給電、Mur)。
*/

#include "ofd.h"
#include "ofd_prototype.h"
#include "ofd_hdf5.h"
#include "ofd_test.h"
#include "hdf5.h"

#include <stdlib.h>
#include <string.h>

#define TMPIN  "ofd_test_solve.ofd"
#define TMPLOG "ofd_test_solve.log"

static char logbuf[65536];

/* 入力 → setup → solve → 後始末を sol_Main と同じ順で行い、ログを logbuf に読む。
   解放の前に呼ぶ check() に結果を渡す */
static void run(const char *solver_line, int maxiter, void (*check)(void))
{
	char text[2048];
	snprintf(text, sizeof(text),
		"OpenFDTD 4 2\n"
		"title = solve test\n"
		"xmesh = -0.03 12 0.03\n"
		"ymesh = -0.03 12 0.03\n"
		"zmesh = -0.03 12 0.03\n"
		"geometry = 1 1 0 0 0 0 -0.0125 0.0125\n"
		"feed = Z 0 0 0 1 0 50\n"
		"frequency1 = 3e9 3e9 0\n"
		"frequency2 = 3e9 3e9 0\n"
		"%s\n"
		"end\n", solver_line);
	FILE *fp = fopen(TMPIN, "w");
	if (fp == NULL) { CHECK(0); return; }
	fputs(text, fp);
	fclose(fp);
	fp = fopen(TMPIN, "r");
	const int ierr = (fp != NULL) ? input_data(fp) : 1;
	if (fp != NULL) fclose(fp);
	remove(TMPIN);
	CHECK(ierr == 0);
	if (ierr) return;

	/* input_data は maxiter を nout の倍数に切り上げる (sol/input_data.c)。
	   solve() 自身は倍数でない maxiter も扱う (最終ステップでも記録する) ので、
	   その分岐を見るときは入力の後・確保の前に直接与える */
	if (maxiter > 0) Solver.maxiter = maxiter;

	Npx = Npy = Npz = 1;
	VECTOR = 0;
	CHECK(setupSize(1, 1, 1, 0) == 0);
	setupABCsize();
	memalloc1();
	memalloc2();
	memalloc3();
	setup();

	FILE *fl = fopen(TMPLOG, "w");
	double tdft = 0;
	remove(OFD_H5_FILE);
	solve(1, &tdft, fl);       /* memfree2 は solve の中で呼ばれる */
	if (fl != NULL) fclose(fl);
	fl = fopen(TMPLOG, "r");
	const size_t n = (fl != NULL) ? fread(logbuf, 1, sizeof(logbuf) - 1, fl) : 0;
	logbuf[n] = '\0';
	if (fl != NULL) fclose(fl);
	remove(TMPLOG);
	CHECK(tdft >= 0);

	check();

	remove(OFD_H5_FILE);
	memfree1();
	memfree3();
	/* memfree は解放してもポインタを残す (本番は 1 回で終わるため)。次の
	   input_data() が realloc で受けるポインタは NULL に戻し、memfree が
	   解放しない入力由来の配列はここで解放する */
	Xn = Yn = Zn = Xc = Yc = Zc = NULL;
	RXn = RYn = RZn = RXc = RYc = RZc = NULL;
	Geometry = NULL; NGeometry = 0;
	Feed = NULL; Point = NULL; Freq1 = Freq2 = NULL;
	Eiter = Hiter = VFeed = IFeed = VPoint = NULL;
	Spara = NULL; cEdft = cHdft = cFdft = NULL;
	cEx_r = cEx_i = cEy_r = cEy_i = cEz_r = cEz_i = NULL;
	cHx_r = cHx_i = cHy_r = cHy_i = cHz_r = cHz_i = NULL;
	free(Material); Material = NULL; NMaterial = 0;
	free(Inductor); Inductor = NULL; NInductor = 0;
	free(Tpa); Tpa = NULL; NTpa = 0;
	free(Gline); Gline = NULL;
	free(MGline); MGline = NULL; NGline = 0;
	NFeed = NPoint = NFreq1 = NFreq2 = 0;
}

static int nsnapshot(void)
{
	int v = -1;
	hid_t f = H5Fopen(OFD_H5_FILE, H5F_ACC_RDONLY, H5P_DEFAULT);
	if (f < 0) return -1;
	hid_t d = H5Dopen(f, "/metadata/Nsnapshot", H5P_DEFAULT);
	if (d >= 0) {
		H5Dread(d, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &v);
		H5Dclose(d);
	}
	H5Fclose(f);
	return v;
}

/* (a) 収束する場合 */
static void check_converged(void)
{
	const int nout = Solver.nout;
	CHECK(strstr(logbuf, "--- converged ---") != NULL);
	CHECK(strstr(logbuf, "--- max steps ---") == NULL);
	CHECK((Ntime > 1) && (Ntime < Solver.maxiter + 1));
	/* 止まるのは記録のステップ (nout の倍数) */
	CHECK((Ntime - 1) % nout == 0);
	CHECK(Niter == (Ntime - 1) / nout + 1);
	/* 最後の記録で初めて E と H の両方が閾値を下回る */
	double emax = 0, hmax = 0;
	int first = -1;
	for (int n = 0; n < Niter; n++) {
		emax = MAX(emax, Eiter[n]);
		hmax = MAX(hmax, Hiter[n]);
		if ((first < 0) && (Eiter[n] < emax * Solver.converg) && (Hiter[n] < hmax * Solver.converg)) first = n;
	}
	CHECK(first == Niter - 1);
	CHECK(emax > 0);
	/* ログの収束履歴の行数 = 記録の数 */
	int nline = 0;
	for (const char *p = logbuf; p != NULL; ) {
		int it; double a, b;
		if (sscanf(p, "%d %lf %lf", &it, &a, &b) == 3) nline++;
		p = strchr(p, '\n');
		if (p != NULL) p++;
	}
	CHECK(nline == Niter);
	/* スナップショットは記録ごとに 1 枚 */
	CHECK(nsnapshot() == Niter);
	/* 給電の時系列は実行したステップぶん埋まり、その先は 0 のまま */
	double vmax = 0, vrest = 0;
	for (int n = 0; n < Ntime; n++) vmax = MAX(vmax, fabs(VFeed[n]));
	for (int n = Ntime; n <= Solver.maxiter; n++) vrest = MAX(vrest, fabs(VFeed[n]));
	CHECK((vmax > 0) && (vrest == 0));
}

/* (b) 収束しない (converg = 0) : 最後まで回り、nout の倍数でない最終ステップでも記録する */
static void check_maxsteps(void)
{
	CHECK(strstr(logbuf, "--- max steps ---") != NULL);
	CHECK(strstr(logbuf, "--- converged ---") == NULL);
	CHECK(Ntime == Solver.maxiter + 1);
	/* maxiter = 25, nout = 10 : 0, 10, 20, 25 の 4 回 */
	CHECK(Niter == 4);
	CHECK(strstr(logbuf, "\n     25 ") != NULL);
	CHECK(nsnapshot() == 4);
	CHECK(fabs(VFeed[Solver.maxiter]) > 0);
}

void test_solve(void)
{
	run("solver = 3000 10 1e-3", 0, check_converged);
	run("solver = 30 10 0", 25, check_maxsteps);
}
