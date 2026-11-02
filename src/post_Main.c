/*
OpenFDTD Version 4.2.1

post process
*/

#define MAIN
#include "ofd.h"
#undef MAIN

#include "ofd_prototype.h"

static void args(int, char *[], int *, int *, char [], char []);

static int UseHdf5 = 0;

int main(int argc, char *argv[])
{
	const char errfmt[] = "*** file %s open error.\n";
	FILE *fp_in = NULL, *fp_out = NULL;

	// arguments
	HTML = 0;
	int nthread = 1;
	int prompt = 0;
	char fn_in[BUFSIZ] = "";
	char fn_out[BUFSIZ] = "ofd.out";
	args(argc, argv, &nthread, &prompt, fn_in, fn_out);

	// set number of threads
#ifdef _OPENMP
	omp_set_num_threads(nthread);
#endif

	// input data
	if ((fp_in = fopen(fn_in, "r")) == NULL) {
		printf(errfmt, fn_in);
		if (prompt) getchar();
		exit(1);
	}
	if (post_data(fp_in)) {
		fclose(fp_in);
		if (prompt) getchar();
		exit(1);
	}

	// read ofd.out (既定) / time_series_data.h5 (-h5)
	if (UseHdf5) {
		if (readhdf5()) {
			if (prompt) getchar();
			exit(1);
		}
	}
	else {
		if ((fp_out = fopen(fn_out, "rb")) == NULL) {
			printf(errfmt, fn_out);
			if (prompt) getchar();
			exit(1);
		}
		readout(fp_out);
		fclose(fp_out);
	}

	// NOTE: 正本 (既定) は ofd.out (readout)。-h5 を付けると time_series_data.h5
	// から読む (readhdf5)。ofd.out と同じ値を /metadata・/geometry・/convergence
	// と /post (余白込みの近傍界 DFT・給電・入力電力) から復元するので、後処理の
	// 出力は同じになる (data/sample/hdf5_post_check.sh が確かめる)。coupling と
	// 散乱断面積は ofd_post がその場で計算し直す (または使わない) ので保存値は要らない。

	// post process
	post();

	// prompt
	if (prompt) getchar();

	return 0;
}


static void args(int argc, char *argv[],
	int *nthread, int *prompt, char fn_in[], char fn_out[])
{
	const char usage[] = "Usage : ofd_post [-n <thread>] [-html] [-out <outfile>] [-h5] <datafile>";

	if (argc < 2) {
		printf("%s\n", usage);
		exit(0);
	}

	while (--argc) {
		++argv;
		if (!strcmp(*argv, "-n")) {
			if (--argc) {
				*nthread = atoi(*++argv);
				if (*nthread < 1) *nthread = 1;
			}
			else {
				break;
			}
		}
		else if (!strcmp(*argv, "-html")) {
			HTML = 1;
		}
		else if (!strcmp(*argv, "-prompt")) {
			*prompt = 1;
		}
		else if (!strcmp(*argv, "-h5")) {
			// ofd.out ではなく time_series_data.h5 から読む
			UseHdf5 = 1;
		}
		else if (!strcmp(*argv, "-out")) {
			argc--;
			strcpy(fn_out, *++argv);
		}
		else if (!strcmp(*argv, "--help")) {
			printf("%s\n", usage);
			exit(0);
		}
		else if (!strcmp(*argv, "--version")) {
			printf("%s Ver.%d.%d.%d\n", PROGRAM, VERSION_MAJOR, VERSION_MINOR, VERSION_BUILD);
			exit(0);
		}
		else {
			strcpy(fn_in, *argv);
		}
	}
}
