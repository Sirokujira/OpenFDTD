/*
test_monitor.c

ログ出力 monitor1..4 (sol/monitor.c) のテスト。

ofd.log の行は契約である: 検証スクリプトが grep で読み (例:
pml_reflection_check.sh は "ABC = PML"、各スクリプトは "normal end")、
OpenFDTD-X も表示に使う。ここでは決まった行がその形で出ることと、
メモリ量・出力ファイルの大きさの見積りが正しい桁であることを確かめる。

見積りは以前 NFreq2 * Nx * Ny * Nz を int で掛けていて、積が 2^31 を超える
大きな問題 (例: 1 周波数 × 1300^3 セル = 2.2e9) で桁あふれしていた。ここでは実際に確保する量
(sol/memallocfree?.c と同じ配列: 場 6 + 材料 id 6 + 近傍界 DFT 12 × NFreq2 面) を
テスト側で 64 ビットで数え、見積りが正で ±5% に入ることを確かめる。
*/

#include "ofd.h"
#include "ofd_prototype.h"
#include "ofd_test.h"

#include <stdlib.h>
#include <string.h>

#define TMPLOG "ofd_test_monitor.log"

static char logbuf[16384];

/* fn() を一時ファイルに書かせて logbuf に読む */
static void capture(void (*fn)(FILE *))
{
	FILE *fp = fopen(TMPLOG, "w");
	if (fp == NULL) { logbuf[0] = '\0'; return; }
	fn(fp);
	fclose(fp);
	fp = fopen(TMPLOG, "r");
	size_t n = (fp != NULL) ? fread(logbuf, 1, sizeof(logbuf) - 1, fp) : 0;
	logbuf[n] = '\0';
	if (fp != NULL) fclose(fp);
	remove(TMPLOG);
}

static int has(const char *line) { return strstr(logbuf, line) != NULL; }

/* "key = <int> [MB]" の値 */
static long long value_mb(const char *key)
{
	const char *p = strstr(logbuf, key);
	if (p == NULL) return -1;
	p = strchr(p, '=');
	return (p != NULL) ? atoll(p + 1) : -1;
}

static void w2(FILE *fp) { monitor2(fp, 0, 1); }
static void w3(FILE *fp) { monitor3(fp, "ofd.log", "ofd.out"); }
static void w4(FILE *fp) { const double cpu[5] = {0, 1, 3, 6, 10}; monitor4(fp, cpu); }

static void set_size(int nx, int ny, int nz, int l)
{
	Nx = nx; Ny = ny; Nz = nz;
	/* setupSize() と同じ (余白 l 層) */
	NN = (int64_t)(Nx + 2 * l + 1) * (Ny + 2 * l + 1) * (Nz + 2 * l + 1);
}

void test_monitor(void)
{
	char title0[256];
	memcpy(title0, Title, sizeof(title0));

	strcpy(Title, "monitor test");
	NMaterial = 3; NGeometry = 2;
	NFeed = 1; NPoint = 2; NFreq1 = 4; NFreq2 = 0;
	IPlanewave = 0;
	iABC = 0; PBCx = PBCy = PBCz = 0;
	VECTOR = 0;
	Dt = 1.5e-12; Tw = 6e-11;
	Solver.maxiter = 1000; Solver.converg = 1e-3;
	set_size(5, 4, 6, 1);

	/* 給電・Mur */
	capture(w2);
	CHECK(has("Title = monitor test\n"));
	CHECK(has("Source = feed\n"));
	CHECK(has("Cells = 5 x 4 x 6 = 120\n"));
	CHECK(has("No. of Materials   = 3\n"));
	CHECK(has("No. of Feeds       = 1\n"));
	CHECK(has("No. of Freq.s (1)  = 4\n"));
	CHECK(has("ABC = Mur-1st\n"));
	CHECK(!has("PBC"));
	CHECK(has("Dt[sec] = 1.5000e-12, Tw[sec] = 6.0000e-11, Tw/Dt = 40.000\n"));
	CHECK(has("Iterations = 1000, Convergence = 1.000e-03\n"));
	CHECK(has("=== iteration start ===\n"));

	/* 平面波・PML・PBC (検証スクリプトが grep する形) */
	NFeed = 0; IPlanewave = 1;
	iABC = 1; cPML.l = 8; cPML.m = 2; cPML.r0 = 1e-6;
	PBCx = 1; PBCz = 1;              /* X と Z (Y と取り違えると "XYZ" などになる) */
	capture(w2);
	CHECK(has("Source = plane wave\n"));
	CHECK(!has("No. of Feeds"));
	CHECK(has("ABC = PML (L=8, M=2.00, R0=1.00e-06)\n"));
	CHECK(has("PBC : XZ\n"));

	/* 大きな問題 : 見積りが正で、実際に確保する量の ±5%。
	   1 周波数 × 1300^3 セル : NFreq2 * Nx * Ny * Nz = 2.2e9 は int を超える。
	   周波数を 1 にして、場と材料 id の分 (全体の半分) も見積りに効くようにする */
	NFeed = 1; IPlanewave = 0;
	iABC = 0; PBCx = PBCy = PBCz = 0;
	NFreq2 = 1;
	set_size(1300, 1300, 1300, 1);
	capture(w2);
	CHECK(has("Cells = 1300 x 1300 x 1300 = 2197000000\n"));
	const double mb = 1024.0 * 1024.0;
	const double alloc = 6.0 * NN * sizeof(real_t) + 6.0 * NN * sizeof(id_t)
	                   + 12.0 * NFreq2 * NN * sizeof(float);
	const long long cpu_mb = value_mb("CPU Memory size");
	CHECK(cpu_mb > 0);
	CHECK(fabs(cpu_mb - alloc / mb) < 0.05 * alloc / mb);
	/* 出力ファイル : 近傍界 DFT 12 面 × NFreq2 (ofd.out の大半) */
	const double out = 12.0 * NFreq2 * NN * sizeof(float);
	const long long out_mb = value_mb("Output filesize");
	CHECK(out_mb > 0);
	CHECK(fabs(out_mb - out / mb) < 0.05 * out / mb);
	if ((cpu_mb <= 0) || (out_mb <= 0)) {
		printf("    CPU Memory %lld MB (expected %.0f), Output %lld MB (expected %.0f)\n",
			cpu_mb, alloc / mb, out_mb, out / mb);
	}

	/* 出力ファイル名と終了 (各スクリプトが "normal end" を grep する) */
	capture(w3);
	CHECK(has("=== output files ===\nofd.log, ofd.out\n"));
	capture(w4);
	CHECK(has("  part-1 :       1.000\n"));
	CHECK(has("  total  :      10.000\n"));
	CHECK(has("=== normal end ===\n"));

	memcpy(Title, title0, sizeof(title0));
	NMaterial = 0; NGeometry = 0;
	NFeed = NPoint = NFreq1 = NFreq2 = 0;
	iABC = 0;
	Nx = Ny = Nz = 0; NN = 0;
}
