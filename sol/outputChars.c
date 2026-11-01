/*
outputChars.c

calculate and output to oth.log
*/

#include "ofd.h"
#include "complex.h"
#include "ofd_prototype.h"


void outputChars(FILE *fp)
{
	// setup far field
	alloc_farfield();
	setup_farfield();

	// input imepedanece
	if (NFeed && NFreq1) {
		calcZin();
		outputZin(fp);
	}

	// input power at frequency2 (far field gain in post, written to ofd.out).
	// frequency1 とは独立 (calcPin は NFeed / NFreq2 を自分で見る)。以前は上の
	// 条件の中にあり、frequency1 の無い給電の入力では Pin が NULL のまま
	// writeout() が NFeed * NFreq2 個を書こうとして落ちていた
	calcPin();

	// S-parameters
	if (NPoint && NFreq1) {
		calcSpara();
		outputSpara(fp);
	}

	// coupling
	if (NFeed && NPoint && NFreq1) {
		outputCoupling(fp);
	}

	// cross section
	if (IPlanewave && NFreq2) {
		outputCross(fp);
	}
}
