/*
comm.c (MPI)

MPI routines
*/

#ifdef _MPI
#include <mpi.h>
#endif

#include "ofd.h"
#include "ofd_prototype.h"

// initialize
void mpi_init(int argc, char **argv)
{
#ifdef _MPI
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &commSize);
	MPI_Comm_rank(MPI_COMM_WORLD, &commRank);
#else
	commSize = 1;
	commRank = 0;
	argc = argc;	// dummy
	argv = argv;	// dummy
#endif
}


// close
void mpi_close(void)
{
#ifdef _MPI
	MPI_Finalize();
#endif
}


// check error code
// mode = 0/1 : Bcast/Allreduce
void comm_check(int ierr, int mode, int prompt)
{
#ifdef _MPI
	if (commSize > 1) {
		if (mode == 0) {
			MPI_Bcast(&ierr, 1, MPI_INT, 0, MPI_COMM_WORLD);
		}
		else {
			int g_ierr;
			MPI_Allreduce(&ierr, &g_ierr, 1, MPI_INT, MPI_LOR, MPI_COMM_WORLD);
			ierr = g_ierr;
		}
	}
	if (ierr) {
		MPI_Finalize();
	}
#endif
	mode = mode;  // dummy
	if (ierr) {
		if (prompt && (commRank == 0)) {
			fflush(stdout);
			getchar();
		}
		exit(1);
	}
}


// gather string
void comm_string(const char *lstr, char *gstr)
{
#ifdef _MPI
	char buff[BUFSIZ];
	if (commRank == 0) {
		MPI_Status status;
		strcpy(gstr, lstr);
		for (int i = 1; i < commSize; i++) {
			MPI_Recv(buff, BUFSIZ, MPI_CHAR, i, 0, MPI_COMM_WORLD, &status);
			strcat(gstr, "\n");
			strcat(gstr, buff);
		}
	}
	else {
		strcpy(buff, lstr);
		MPI_Send(buff, BUFSIZ, MPI_CHAR, 0, 0, MPI_COMM_WORLD);
	}
#else
	strcpy(gstr, lstr);
#endif
}


// get cpu time [sec]
double comm_cputime(void)
{
#ifdef _MPI
	MPI_Barrier(MPI_COMM_WORLD);
	return MPI_Wtime();
#else
#ifdef _WIN32
	return (double)clock() / CLOCKS_PER_SEC;
#else
	struct timespec ts;
	clock_gettime(CLOCK_REALTIME, &ts);
	return (ts.tv_sec + (ts.tv_nsec * 1e-9));
#endif  // _WIN32
#endif  // MPI
}


// broadcast input data
void comm_broadcast(void)
{
#ifdef _MPI
	int    *i_buf = NULL;
	double *d_buf = NULL;
	char   *c_buf = NULL;
	int    i_num = 0;
	int    d_num = 0;
	int    c_num = 0;

	// variables to buffers (root only)

	if (commRank == 0) {
		// number of data
		i_num = 18 + (1 * (int)NMaterial) + (2 * (int)NGeometry);
		d_num =  6 + (8 * (int)NMaterial) + (8 * (int)NGeometry) + NFreq1 + NFreq2 + (Nx + 1) + (Ny + 1) + (Nz + 1) + Nx + Ny + Nz;
		c_num =  0;
		if (NFeed > 0) {
			i_num += 3 * NFeed;
			d_num += 5 * NFeed;
			c_num += 1 * NFeed;
		}
		if (IPlanewave) {
			d_num += 13;
		}
		if (NPoint > 0) {
			i_num += 3 * (NPoint + 2);
			d_num += 3 * (NPoint + 2);
			c_num += 1 * (NPoint + 2);
		}
		if (NInductor > 0) {
			i_num += 3 * NInductor;
			d_num += 6 * NInductor;
			c_num += 1 * NInductor;
		}
		// HDF5 出力 (hdf5 キー) : Hdf5Output, Hdf5Interval
		i_num += 2;

		// TPA (二光子吸収) : NTpa + material id、および CW 波源のパラメータ
		i_num += 1 + (1 * NTpa);
		d_num += 2 + (1 * NTpa);

		// alloc
		i_buf = (int *)   malloc(i_num * sizeof(int));
		d_buf = (double *)malloc(d_num * sizeof(double));
		c_buf = (char *)  malloc(c_num * sizeof(char));

		int i_id = 0;
		int d_id = 0;
		int c_id = 0;

		i_buf[i_id++] = Nx;
		i_buf[i_id++] = Ny;
		i_buf[i_id++] = Nz;
		i_buf[i_id++] = (int)NMaterial;
		i_buf[i_id++] = (int)NGeometry;
		i_buf[i_id++] = NFeed;
		i_buf[i_id++] = IPlanewave;
		i_buf[i_id++] = iABC;
		i_buf[i_id++] = cPML.l;
		i_buf[i_id++] = PBCx;
		i_buf[i_id++] = PBCy;
		i_buf[i_id++] = PBCz;
		i_buf[i_id++] = NFreq1;
		i_buf[i_id++] = NFreq2;
		i_buf[i_id++] = Solver.maxiter;
		i_buf[i_id++] = Solver.nout;
		i_buf[i_id++] = NPoint;
		i_buf[i_id++] = NInductor;

		d_buf[d_id++] = rFeed;
		d_buf[d_id++] = cPML.m;
		d_buf[d_id++] = cPML.r0;
		d_buf[d_id++] = Solver.converg;
		d_buf[d_id++] = Dt;
		d_buf[d_id++] = Tw;

		for (int i = 0; i <= Nx; i++) {
			d_buf[d_id++] = Xn[i];
		}
		for (int j = 0; j <= Ny; j++) {
			d_buf[d_id++] = Yn[j];
		}
		for (int k = 0; k <= Nz; k++) {
			d_buf[d_id++] = Zn[k];
		}

		for (int i = 0; i < Nx; i++) {
			d_buf[d_id++] = Xc[i];
		}
		for (int j = 0; j < Ny; j++) {
			d_buf[d_id++] = Yc[j];
		}
		for (int k = 0; k < Nz; k++) {
			d_buf[d_id++] = Zc[k];
		}

		for (int n = 0; n < NMaterial; n++) {
			i_buf[i_id++] = Material[n].type;
			d_buf[d_id++] = Material[n].epsr;
			d_buf[d_id++] = Material[n].esgm;
			d_buf[d_id++] = Material[n].amur;
			d_buf[d_id++] = Material[n].msgm;
			d_buf[d_id++] = Material[n].einf;
			d_buf[d_id++] = Material[n].ae;
			d_buf[d_id++] = Material[n].be;
			d_buf[d_id++] = Material[n].ce;
		}

		for (int n = 0; n < NGeometry; n++) {
			i_buf[i_id++] = (int)Geometry[n].m;
			i_buf[i_id++] = Geometry[n].shape;
			for (int i = 0; i < 8; i++) {
				d_buf[d_id++] = Geometry[n].g[i];
			}
		}

		for (int n = 0; n < NFeed; n++) {
			c_buf[c_id++] = Feed[n].dir;
			i_buf[i_id++] = Feed[n].i;
			i_buf[i_id++] = Feed[n].j;
			i_buf[i_id++] = Feed[n].k;
			d_buf[d_id++] = Feed[n].volt;
			d_buf[d_id++] = Feed[n].delay;
			d_buf[d_id++] = Feed[n].dx;
			d_buf[d_id++] = Feed[n].dy;
			d_buf[d_id++] = Feed[n].dz;
		}

		if (IPlanewave) {
			for (int m = 0; m < 3; m++) {
				d_buf[d_id++] = Planewave.ei[m];
				d_buf[d_id++] = Planewave.hi[m];
				d_buf[d_id++] = Planewave.ri[m];
				d_buf[d_id++] = Planewave.r0[m];
			}
			d_buf[d_id++] = Planewave.ai;
		}

		if (NPoint > 0) {
			for (int n = 0; n < NPoint + 2; n++) {
				c_buf[c_id++] = Point[n].dir;
				i_buf[i_id++] = Point[n].i;
				i_buf[i_id++] = Point[n].j;
				i_buf[i_id++] = Point[n].k;
				d_buf[d_id++] = Point[n].dx;
				d_buf[d_id++] = Point[n].dy;
				d_buf[d_id++] = Point[n].dz;
			}
		}

		for (int n = 0; n < NInductor; n++) {
			c_buf[c_id++] = Inductor[n].dir;
			i_buf[i_id++] = Inductor[n].i;
			i_buf[i_id++] = Inductor[n].j;
			i_buf[i_id++] = Inductor[n].k;
			d_buf[d_id++] = Inductor[n].dx;
			d_buf[d_id++] = Inductor[n].dy;
			d_buf[d_id++] = Inductor[n].dz;
			d_buf[d_id++] = Inductor[n].fctr;
			d_buf[d_id++] = Inductor[n].e;
			d_buf[d_id++] = Inductor[n].esum;
		}

		// HDF5 出力 (hdf5 キー)
		i_buf[i_id++] = Hdf5Output;
		i_buf[i_id++] = Hdf5Interval;

		// TPA (二光子吸収)
		i_buf[i_id++] = NTpa;
		d_buf[d_id++] = WaveAmp;
		d_buf[d_id++] = WaveOmega;
		for (int n = 0; n < NTpa; n++) {
			i_buf[i_id++] = (int)Tpa[n].m;
			d_buf[d_id++] = Tpa[n].beta;
		}

		for (int n = 0; n < NFreq1; n++) {
			d_buf[d_id++] = Freq1[n];
		}

		for (int n = 0; n < NFreq2; n++) {
			d_buf[d_id++] = Freq2[n];
		}

		// check
		assert(i_id == i_num);
		assert(d_id == d_num);
		assert(c_id == c_num);
	}

	// broadcast (root to non-root)

	MPI_Bcast(&i_num, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&d_num, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&c_num, 1, MPI_INT, 0, MPI_COMM_WORLD);
	//printf("%d %d %d %d %d\n", commSize, commRank, i_num, d_num, c_num); fflush(stdout);

	// alloc
	if (commRank > 0) {
		i_buf = (int *)   malloc(i_num * sizeof(int));
		d_buf = (double *)malloc(d_num * sizeof(double));
		c_buf = (char *)  malloc(c_num * sizeof(char));
	}

	MPI_Bcast(i_buf, i_num, MPI_INT,    0, MPI_COMM_WORLD);
	MPI_Bcast(d_buf, d_num, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(c_buf, c_num, MPI_CHAR,   0, MPI_COMM_WORLD);

	// buffers to variables (non-root)

	if (commRank > 0) {
		int i_id = 0;
		//int f_id = 0;
		int d_id = 0;
		int c_id = 0;

		Nx              = i_buf[i_id++];
		Ny              = i_buf[i_id++];
		Nz              = i_buf[i_id++];
		NMaterial       = i_buf[i_id++];
		NGeometry       = i_buf[i_id++];
		NFeed           = i_buf[i_id++];
		IPlanewave      = i_buf[i_id++];
		iABC            = i_buf[i_id++];
		cPML.l          = i_buf[i_id++];
		PBCx            = i_buf[i_id++];
		PBCy            = i_buf[i_id++];
		PBCz            = i_buf[i_id++];
		NFreq1          = i_buf[i_id++];
		NFreq2          = i_buf[i_id++];
		Solver.maxiter  = i_buf[i_id++];
		Solver.nout     = i_buf[i_id++];
		NPoint          = i_buf[i_id++];
		NInductor       = i_buf[i_id++];

		rFeed           = d_buf[d_id++];
		cPML.m          = d_buf[d_id++];
		cPML.r0         = d_buf[d_id++];
		Solver.converg  = d_buf[d_id++];
		Dt              = d_buf[d_id++];
		Tw              = d_buf[d_id++];

		Xn = (double *)malloc((Nx + 1) * sizeof(double));
		Yn = (double *)malloc((Ny + 1) * sizeof(double));
		Zn = (double *)malloc((Nz + 1) * sizeof(double));
		for (int i = 0; i <= Nx; i++) {
			Xn[i] = d_buf[d_id++];
		}
		for (int j = 0; j <= Ny; j++) {
			Yn[j] = d_buf[d_id++];
		}
		for (int k = 0; k <= Nz; k++) {
			Zn[k] = d_buf[d_id++];
		}

		Xc = (double *)malloc(Nx * sizeof(double));
		Yc = (double *)malloc(Ny * sizeof(double));
		Zc = (double *)malloc(Nz * sizeof(double));
		for (int i = 0; i < Nx; i++) {
			Xc[i] = d_buf[d_id++];
		}
		for (int j = 0; j < Ny; j++) {
			Yc[j] = d_buf[d_id++];
		}
		for (int k = 0; k < Nz; k++) {
			Zc[k] = d_buf[d_id++];
		}

		if (NMaterial > 0) {
			Material = (material_t *)malloc(NMaterial * sizeof(material_t));
			for (int n = 0; n < NMaterial; n++) {
				Material[n].type = i_buf[i_id++];
				Material[n].epsr = d_buf[d_id++];
				Material[n].esgm = d_buf[d_id++];
				Material[n].amur = d_buf[d_id++];
				Material[n].msgm = d_buf[d_id++];
				Material[n].einf = d_buf[d_id++];
				Material[n].ae   = d_buf[d_id++];
				Material[n].be   = d_buf[d_id++];
				Material[n].ce   = d_buf[d_id++];
			}
		}

		if (NGeometry > 0){
			Geometry = (geometry_t *)malloc(NGeometry * sizeof(geometry_t));
			for (int n = 0; n < NGeometry; n++) {
				Geometry[n].m     = (id_t)i_buf[i_id++];
				Geometry[n].shape = i_buf[i_id++];
				for (int i = 0; i < 8; i++) {
					Geometry[n].g[i] = d_buf[d_id++];
				}
			}
		}

		if (NFeed > 0) {
			Feed = (feed_t *)malloc(NFeed * sizeof(feed_t));
			for (int n = 0; n < NFeed; n++) {
				Feed[n].dir   = c_buf[c_id++];
				Feed[n].i     = i_buf[i_id++];
				Feed[n].j     = i_buf[i_id++];
				Feed[n].k     = i_buf[i_id++];
				Feed[n].volt  = d_buf[d_id++];
				Feed[n].delay = d_buf[d_id++];
				Feed[n].dx    = d_buf[d_id++];
				Feed[n].dy    = d_buf[d_id++];
				Feed[n].dz    = d_buf[d_id++];
			}
		}

		if (IPlanewave) {
			for (int m = 0; m < 3; m++) {
				Planewave.ei[m] = d_buf[d_id++];
				Planewave.hi[m] = d_buf[d_id++];
				Planewave.ri[m] = d_buf[d_id++];
				Planewave.r0[m] = d_buf[d_id++];
			}
			Planewave.ai = d_buf[d_id++];
		}

		if (NPoint > 0) {
			Point = (point_t *)malloc((NPoint + 2) * sizeof(point_t));
			for (int n = 0; n < NPoint + 2; n++) {
				Point[n].dir = c_buf[c_id++];
				Point[n].i   = i_buf[i_id++];
				Point[n].j   = i_buf[i_id++];
				Point[n].k   = i_buf[i_id++];
				Point[n].dx  = d_buf[d_id++];
				Point[n].dy  = d_buf[d_id++];
				Point[n].dz  = d_buf[d_id++];
			}
		}

		if (NInductor > 0) {
			Inductor = (inductor_t *)malloc(NInductor * sizeof(inductor_t));
			for (int n = 0; n < NInductor; n++) {
				Inductor[n].dir  = c_buf[c_id++];
				Inductor[n].i    = i_buf[i_id++];
				Inductor[n].j    = i_buf[i_id++];
				Inductor[n].k    = i_buf[i_id++];
				Inductor[n].dx   = d_buf[d_id++];
				Inductor[n].dy   = d_buf[d_id++];
				Inductor[n].dz   = d_buf[d_id++];
				Inductor[n].fctr = d_buf[d_id++];
				Inductor[n].e    = d_buf[d_id++];
				Inductor[n].esum = d_buf[d_id++];
			}
		}

		// TPA (二光子吸収)
		// これを配り忘れると rank != 0 で NTpa = 0 となり TPA 処理を丸ごと
		// スキップするため、ランク間で集団操作の呼び出しが食い違って
		// デッドロックする (実際に踏んだ)。
		Hdf5Output   = i_buf[i_id++];
		Hdf5Interval = i_buf[i_id++];

		NTpa      = i_buf[i_id++];
		WaveAmp   = d_buf[d_id++];
		WaveOmega = d_buf[d_id++];
		if (NTpa > 0) {
			Tpa = (tpa_t *)malloc(NTpa * sizeof(tpa_t));
			for (int n = 0; n < NTpa; n++) {
				Tpa[n].m    = (id_t)i_buf[i_id++];
				Tpa[n].beta = d_buf[d_id++];
			}
		}

		if (NFreq1 > 0) {
			Freq1 = (double *)malloc(NFreq1 * sizeof(double));
			for (int n = 0; n < NFreq1; n++) {
				Freq1[n] = d_buf[d_id++];
			}
		}

		if (NFreq2 > 0) {
			Freq2 = (double *)malloc(NFreq2 * sizeof(double));
			for (int n = 0; n < NFreq2; n++) {
				Freq2[n] = d_buf[d_id++];
			}
		}

		// check
		assert(i_id == i_num);
		assert(d_id == d_num);
		assert(c_id == c_num);
	}

	// free
	free(i_buf);
	free(d_buf);
	free(c_buf);

	// debug
	//printf("%d %d %d %d\n", commSize, commRank, iSIMD, nThread);
	//printf("%d %d %d %d %d %d %d\n", Nx, Ny, Nz, NI, NJ, N0, NN);
	//printf("%d %d %d\n", NMaterial, NGeometry, NFeed);
	//printf("%d %d %e\n", Solver.maxiter, Solver.nout, Solver.converg);
	//for (int i = 0; i <= Nx; i++) printf("%d Xn[%d]=%.5f\n", commRank, i, Xn[i] * 1e3);
	//for (int j = 0; j <= Ny; j++) printf("%d Yn[%d]=%.5f\n", commRank, j, Yn[j] * 1e3);
	//for (int k = 0; k <= Nz; k++) printf("%d Zn[%d]=%.5f\n", commRank, k, Zn[k] * 1e3);
	//for (int n = 0; n < NMaterial; n++) printf("%d %d %e %e %e %e\n", commRank, n, Material[n].epsr, Material[n].esgm, Material[n].amur, Material[n].msgm);
	//for (int n = 0; n < NGeometry; n++) printf("%d %d %d %e %e %e %e %e %e\n", commRank, n, Geometry[n].m, Geometry[n].g[0], Geometry[n].g[1], Geometry[n].g[2], Geometry[n].g[3], Geometry[n].g[4], Geometry[n].g[5]);
	//for (int n = 0; n < NFeed; n++) printf("%d %d %c %d %d %d %e\n", commRank, n, Feed[n].dir, Feed[n].i, Feed[n].j, Feed[n].k, Feed[n].volt);
	//for (int n = 0; n < NFreq1; n++) printf("%d %d %e\n", commRank, n, Freq1[n]);
	//for (int n = 0; n < NFreq2; n++) printf("%d %d %e\n", commRank, n, Freq2[n]);
	//for (int n = 0; n < NPoint + 2; n++) printf("%d %d %c %d %d %d\n", commRank, n, Point[n].dir, Point[n].i, Point[n].j, Point[n].k);
	//for (int n = 0; n < NInductor; n++) printf("%d %d %c %d %d %d %e\n", commRank, n, Inductor[n].dir, Inductor[n].i, Inductor[n].j, Inductor[n].k, Inductor[n].fctr);
	fflush(stdout);

#endif
}


// allreduce average
void comm_average(double fsum[])
{
#ifdef _MPI
	double ftmp[2];

	MPI_Allreduce(fsum, ftmp, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

	fsum[0] = ftmp[0];
	fsum[1] = ftmp[1];
#else
	fsum = fsum;	// dummy
#endif
}


// send feed waveform to root process
void comm_feed(void)
{
#ifdef _MPI
	MPI_Status status;

	const int count = Solver.maxiter + 1;
	for (int n = 0; n < NFeed; n++) {
		const int i = Feed[n].i;
		const int j = Feed[n].j;
		const int k = Feed[n].k;
		//printf("%d %d %d\n", i, j, k); fflush(stdout);
		//printf("%d %d\n", commRank, comm_inproc(i, j, k)); fflush(stdout);
		// non-root only
		if      ((commRank == 0) && !comm_inproc(i, j, k)) {
			MPI_Recv(&VFeed[n * count], count, MPI_DOUBLE, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);
			MPI_Recv(&IFeed[n * count], count, MPI_DOUBLE, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);
		}
		else if ((commRank  > 0) &&  comm_inproc(i, j, k)) {
			MPI_Send(&VFeed[n * count], count, MPI_DOUBLE, 0,              0, MPI_COMM_WORLD);
			MPI_Send(&IFeed[n * count], count, MPI_DOUBLE, 0,              0, MPI_COMM_WORLD);
		}
		MPI_Barrier(MPI_COMM_WORLD);
	}
#endif
}


// send point waveform to root process
void comm_point(void)
{
#ifdef _MPI
	MPI_Status status;

	if (NPoint <= 0) return;

	const int count = Solver.maxiter + 1;
	for (int n = 0; n < NPoint + 2; n++) {
		const int i = Point[n].i;
		const int j = Point[n].j;
		const int k = Point[n].k;
		// non-root only
		if      ((commRank == 0) && !comm_inproc(i, j, k)) {
			MPI_Recv(&VPoint[n * count], count, MPI_DOUBLE, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);
		}
		else if ((commRank  > 0) &&  comm_inproc(i, j, k)) {
			MPI_Send(&VPoint[n * count], count, MPI_DOUBLE, 0,              0, MPI_COMM_WORLD);
		}
		MPI_Barrier(MPI_COMM_WORLD);
	}
#endif
}


// send near3d data to root process
void comm_near3d(void)
{
#ifdef _MPI
	const int tag = 0;
	MPI_Status status;
	int isend[11], irecv[11];
	size_t size;
	int64_t g_n, l_n;
	//printf("%d %d %zd\n", commSize, commRank, NN); fflush(stdout);

	if ((NN <= 0) || (NFreq2 <= 0)) return;

	// root : self copy to global array
	if (commRank == 0) {
		// save local index
		const int imin = iMin;
		const int imax = iMax;
		const int jmin = jMin;
		const int jmax = jMax;
		const int kmin = kMin;
		const int kmax = kMax;
		const int64_t ni = Ni;
		const int64_t nj = Nj;
		const int64_t nk = Nk;
		const int64_t n0 = N0;
		const int64_t nn = NN;

		// new global index (Ni, Nj, Nk, N0, NN)
		setupSize(1, 1, 1, 0);

		// alloc global array
		size = NN * NFreq2 * sizeof(float);
		g_cEx_r = (float *)malloc(size);
		g_cEx_i = (float *)malloc(size);
		g_cEy_r = (float *)malloc(size);
		g_cEy_i = (float *)malloc(size);
		g_cEz_r = (float *)malloc(size);
		g_cEz_i = (float *)malloc(size);
		g_cHx_r = (float *)malloc(size);
		g_cHx_i = (float *)malloc(size);
		g_cHy_r = (float *)malloc(size);
		g_cHy_i = (float *)malloc(size);
		g_cHz_r = (float *)malloc(size);
		g_cHz_i = (float *)malloc(size);
		// ゼロ初期化する。Yee 格子では最外縁の一部の成分 (Ex の i=Nx,
		// Ey の j=Ny, Ez の k=Nz) が更新対象外で、以下のコピー/受信でも
		// 埋められない。非 MPI 版では確保時に 0 になっているので、
		// 初期化を省くとその成分だけ MPI 版が未初期化値になる
		// (HDF5 に全域を出すようになって顕在化した)。
		memset(g_cEx_r, 0, size);  memset(g_cEx_i, 0, size);
		memset(g_cEy_r, 0, size);  memset(g_cEy_i, 0, size);
		memset(g_cEz_r, 0, size);  memset(g_cEz_i, 0, size);
		memset(g_cHx_r, 0, size);  memset(g_cHx_i, 0, size);
		memset(g_cHy_r, 0, size);  memset(g_cHy_i, 0, size);
		memset(g_cHz_r, 0, size);  memset(g_cHz_i, 0, size);
		// self copy
		for (int ifreq = 0; ifreq < NFreq2; ifreq++) {
			for (int i = imin; i <  imax; i++) {
			for (int j = jmin; j <= jmax; j++) {
			for (int k = kmin; k <= kmax; k++) {
				g_n = (Ni * i) + (Nj * j) + (Nk * k) + N0 + (ifreq * NN);
				l_n = (ni * i) + (nj * j) + (nk * k) + n0 + (ifreq * nn);
				g_cEx_r[g_n] = cEx_r[l_n];
				g_cEx_i[g_n] = cEx_i[l_n];
			}
			}
			}
			for (int i = imin; i <= imax; i++) {
			for (int j = jmin; j <  jmax; j++) {
			for (int k = kmin; k <= kmax; k++) {
				g_n = (Ni * i) + (Nj * j) + (Nk * k) + N0 + (ifreq * NN);
				l_n = (ni * i) + (nj * j) + (nk * k) + n0 + (ifreq * nn);
				g_cEy_r[g_n] = cEy_r[l_n];
				g_cEy_i[g_n] = cEy_i[l_n];
			}
			}
			}
			for (int i = imin; i <= imax; i++) {
			for (int j = jmin; j <= jmax; j++) {
			for (int k = kmin; k <  kmax; k++) {
				g_n = (Ni * i) + (Nj * j) + (Nk * k) + N0 + (ifreq * NN);
				l_n = (ni * i) + (nj * j) + (nk * k) + n0 + (ifreq * nn);
				g_cEz_r[g_n] = cEz_r[l_n];
				g_cEz_i[g_n] = cEz_i[l_n];
			}
			}
			}
			for (int i = imin - 0; i <= imax; i++) {
			for (int j = jmin - 1; j <= jmax; j++) {
			for (int k = kmin - 1; k <= kmax; k++) {
				g_n = (Ni * i) + (Nj * j) + (Nk * k) + N0 + (ifreq * NN);
				l_n = (ni * i) + (nj * j) + (nk * k) + n0 + (ifreq * nn);
				g_cHx_r[g_n] = cHx_r[l_n];
				g_cHx_i[g_n] = cHx_i[l_n];
			}
			}
			}
			for (int i = imin - 1; i <= imax; i++) {
			for (int j = jmin - 0; j <= jmax; j++) {
			for (int k = kmin - 1; k <= kmax; k++) {
				g_n = (Ni * i) + (Nj * j) + (Nk * k) + N0 + (ifreq * NN);
				l_n = (ni * i) + (nj * j) + (nk * k) + n0 + (ifreq * nn);
				g_cHy_r[g_n] = cHy_r[l_n];
				g_cHy_i[g_n] = cHy_i[l_n];
			}
			}
			}
			for (int i = imin - 1; i <= imax; i++) {
			for (int j = jmin - 1; j <= jmax; j++) {
			for (int k = kmin - 0; k <= kmax; k++) {
				g_n = (Ni * i) + (Nj * j) + (Nk * k) + N0 + (ifreq * NN);
				l_n = (ni * i) + (nj * j) + (nk * k) + n0 + (ifreq * nn);
				g_cHz_r[g_n] = cHz_r[l_n];
				g_cHz_i[g_n] = cHz_i[l_n];
			}
			}
			}
		}
	}

	// loop on freqencies
	for (int ifreq = 0; ifreq < NFreq2; ifreq++) {

		// root
		if (commRank == 0) {
			for (int irank = 1; irank < commSize; irank++) {
				// receive local index
				MPI_Recv(irecv, 11, MPI_INT, irank, tag, MPI_COMM_WORLD, &status);
				const int imin = irecv[0];
				const int imax = irecv[1];
				const int jmin = irecv[2];
				const int jmax = irecv[3];
				const int kmin = irecv[4];
				const int kmax = irecv[5];
				const int ni   = irecv[6];
				const int nj   = irecv[7];
				const int nk   = irecv[8];
				const int n0   = irecv[9];
				const int nn   = irecv[10];
				//printf("%d %d %d %d %d\n", commRank, imin, imax, n0, nn); fflush(stdout);

				// alloc receive buffer
				size = nn * sizeof(float);
				float *recv_r = (float *)malloc(size);
				float *recv_i = (float *)malloc(size);

				// === copy from local array to global array ===

				// Ex
				MPI_Recv(recv_r, nn, MPI_FLOAT, irank, tag, MPI_COMM_WORLD, &status);
				MPI_Recv(recv_i, nn, MPI_FLOAT, irank, tag, MPI_COMM_WORLD, &status);
				for (int i = imin; i <  imax; i++) {
				for (int j = jmin; j <= jmax; j++) {
				for (int k = kmin; k <= kmax; k++) {
					g_n = (Ni * i) + (Nj * j) + (Nk * k) + N0 + (ifreq * NN);
					l_n = (ni * i) + (nj * j) + (nk * k) + n0;
					g_cEx_r[g_n] = recv_r[l_n];
					g_cEx_i[g_n] = recv_i[l_n];
				}
				}
				}

				// Ey
				MPI_Recv(recv_r, nn, MPI_FLOAT, irank, tag, MPI_COMM_WORLD, &status);
				MPI_Recv(recv_i, nn, MPI_FLOAT, irank, tag, MPI_COMM_WORLD, &status);
				for (int i = imin; i <= imax; i++) {
				for (int j = jmin; j <  jmax; j++) {
				for (int k = kmin; k <= kmax; k++) {
					g_n = (Ni * i) + (Nj * j) + (Nk * k) + N0 + (ifreq * NN);
					l_n = (ni * i) + (nj * j) + (nk * k) + n0;
					g_cEy_r[g_n] = recv_r[l_n];
					g_cEy_i[g_n] = recv_i[l_n];
				}
				}
				}

				// Ez
				MPI_Recv(recv_r, nn, MPI_FLOAT, irank, tag, MPI_COMM_WORLD, &status);
				MPI_Recv(recv_i, nn, MPI_FLOAT, irank, tag, MPI_COMM_WORLD, &status);
				for (int i = imin; i <= imax; i++) {
				for (int j = jmin; j <= jmax; j++) {
				for (int k = kmin; k <  kmax; k++) {
					g_n = (Ni * i) + (Nj * j) + (Nk * k) + N0 + (ifreq * NN);
					l_n = (ni * i) + (nj * j) + (nk * k) + n0;
					g_cEz_r[g_n] = recv_r[l_n];
					g_cEz_i[g_n] = recv_i[l_n];
				}
				}
				}

				// Hx
				MPI_Recv(recv_r, nn, MPI_FLOAT, irank, tag, MPI_COMM_WORLD, &status);
				MPI_Recv(recv_i, nn, MPI_FLOAT, irank, tag, MPI_COMM_WORLD, &status);
				for (int i = imin - 0; i <= imax; i++) {
				for (int j = jmin - 1; j <= jmax; j++) {
				for (int k = kmin - 1; k <= kmax; k++) {
					g_n = (Ni * i) + (Nj * j) + (Nk * k) + N0 + (ifreq * NN);
					l_n = (ni * i) + (nj * j) + (nk * k) + n0;
					g_cHx_r[g_n] = recv_r[l_n];
					g_cHx_i[g_n] = recv_i[l_n];
				}
				}
				}

				// Hy
				MPI_Recv(recv_r, nn, MPI_FLOAT, irank, tag, MPI_COMM_WORLD, &status);
				MPI_Recv(recv_i, nn, MPI_FLOAT, irank, tag, MPI_COMM_WORLD, &status);
				for (int i = imin - 1; i <= imax; i++) {
				for (int j = jmin - 0; j <= jmax; j++) {
				for (int k = kmin - 1; k <= kmax; k++) {
					g_n = (Ni * i) + (Nj * j) + (Nk * k) + N0 + (ifreq * NN);
					l_n = (ni * i) + (nj * j) + (nk * k) + n0;
					g_cHy_r[g_n] = recv_r[l_n];
					g_cHy_i[g_n] = recv_i[l_n];
				}
				}
				}

				// Hz
				MPI_Recv(recv_r, nn, MPI_FLOAT, irank, tag, MPI_COMM_WORLD, &status);
				MPI_Recv(recv_i, nn, MPI_FLOAT, irank, tag, MPI_COMM_WORLD, &status);
				for (int i = imin - 1; i <= imax; i++) {
				for (int j = jmin - 1; j <= jmax; j++) {
				for (int k = kmin - 0; k <= kmax; k++) {
					g_n = (Ni * i) + (Nj * j) + (Nk * k) + N0 + (ifreq * NN);
					l_n = (ni * i) + (nj * j) + (nk * k) + n0;
					g_cHz_r[g_n] = recv_r[l_n];
					g_cHz_i[g_n] = recv_i[l_n];
				}
				}
				}

				// free
				free(recv_r);
				free(recv_i);
			}
		}

		// non-root
		else {
			// send to root
			// index
			isend[0]  = iMin;
			isend[1]  = iMax;
			isend[2]  = jMin;
			isend[3]  = jMax;
			isend[4]  = kMin;
			isend[5]  = kMax;
			isend[6]  = (int)Ni;
			isend[7]  = (int)Nj;
			isend[8]  = (int)Nk;
			isend[9]  = (int)N0;
			isend[10] = (int)NN;
			MPI_Send(isend, 11, MPI_INT, 0, tag, MPI_COMM_WORLD);

			// data
			const int64_t ns = ifreq * NN;
			const int nn = (int)NN;
			assert(nn > 0);
			MPI_Send(&cEx_r[ns], nn, MPI_FLOAT, 0, tag, MPI_COMM_WORLD);
			MPI_Send(&cEx_i[ns], nn, MPI_FLOAT, 0, tag, MPI_COMM_WORLD);
			MPI_Send(&cEy_r[ns], nn, MPI_FLOAT, 0, tag, MPI_COMM_WORLD);
			MPI_Send(&cEy_i[ns], nn, MPI_FLOAT, 0, tag, MPI_COMM_WORLD);
			MPI_Send(&cEz_r[ns], nn, MPI_FLOAT, 0, tag, MPI_COMM_WORLD);
			MPI_Send(&cEz_i[ns], nn, MPI_FLOAT, 0, tag, MPI_COMM_WORLD);
			MPI_Send(&cHx_r[ns], nn, MPI_FLOAT, 0, tag, MPI_COMM_WORLD);
			MPI_Send(&cHx_i[ns], nn, MPI_FLOAT, 0, tag, MPI_COMM_WORLD);
			MPI_Send(&cHy_r[ns], nn, MPI_FLOAT, 0, tag, MPI_COMM_WORLD);
			MPI_Send(&cHy_i[ns], nn, MPI_FLOAT, 0, tag, MPI_COMM_WORLD);
			MPI_Send(&cHz_r[ns], nn, MPI_FLOAT, 0, tag, MPI_COMM_WORLD);
			MPI_Send(&cHz_i[ns], nn, MPI_FLOAT, 0, tag, MPI_COMM_WORLD);

		}
	}

	// copy pointer : global -> root
	if (commRank == 0) {
		free(cEx_r);
		free(cEx_i);
		free(cEy_r);
		free(cEy_i);
		free(cEz_r);
		free(cEz_i);
		free(cHx_r);
		free(cHx_i);
		free(cHy_r);
		free(cHy_i);
		free(cHz_r);
		free(cHz_i);

		cEx_r = g_cEx_r;
		cEx_i = g_cEx_i;
		cEy_r = g_cEy_r;
		cEy_i = g_cEy_i;
		cEz_r = g_cEz_r;
		cEz_i = g_cEz_i;
		cHx_r = g_cHx_r;
		cHx_i = g_cHx_i;
		cHy_r = g_cHy_r;
		cHy_i = g_cHy_i;
		cHz_r = g_cHz_r;
		cHz_i = g_cHz_i;
	}
	//printf("%d %d\n", commSize, commRank); fflush(stdout);
#endif
}


// my rank include index (i, j, k) ?
int comm_inproc(int i, int j, int k)
{
#ifdef _MPI
	const int b = (iMin <= i) && (i < iMax)
               && (jMin <= j) && (j < jMax)
               && (kMin <= k) && (k < kMax);
	return b || ((commRank == commSize - 1) && ((i == Nx) || (j == Ny) || (k == Nz)));
#else
	return i && j && k;  // dummy
#endif
}


// 瞬時値スナップショットを rank 0 の全域配列 (g_Ex..g_Hz) に集める
//
// comm_near3d() は解析の最後に 1 度だけ呼ばれ、rank 0 の索引を
// setupSize(1,1,1,0) で全域に張り替えて終わる (以後ローカル索引は使わない)。
// スナップショットは時間ループの途中で何度も呼ぶのでその手は使えず、
// ローカル索引を保ったまま全域配列へ詰める。
//
// **全ランクが必ず同じ回数・同じ順序で呼ぶこと。** rank 条件の中に入れると
// 送受信の相手が食い違ってデッドロックする (TPA で実際に踏んだ轍)。
//
// 各成分の範囲は comm_near3d() と同一にしてある。Yee 格子で更新対象外に
// なる最外縁 (Ex の i=Nx, Ey の j=Ny, Ez の k=Nz) は書かれないので、
// 確保時のゼロ初期化がそのまま残る (非 MPI 版と同じ値になる)。
void comm_snapshot(void)
{
#ifdef _MPI
	const int tag = 0;
	MPI_Status status;
	int isend[11], irecv[11];
	int64_t g_n, l_n;

	if (NN <= 0) return;

	// 全域の索引 (Ni, Nj, Nk, N0, NN) を 1 度だけ取得して覚えておく。
	// setupSize() を使うので添字の計算式を二重に持たない。
	static int64_t G_Ni = 0, G_Nj = 0, G_Nk = 0, G_N0 = 0, G_NN = 0;
	if (G_NN <= 0) {
		setupSize(1, 1, 1, 0);
		G_Ni = Ni;  G_Nj = Nj;  G_Nk = Nk;  G_N0 = N0;  G_NN = NN;
		setupSize(Npx, Npy, Npz, commRank);   // ローカル索引に戻す

		if (commRank == 0) {
			const size_t size = G_NN * sizeof(real_t);
			g_Ex = (real_t *)malloc(size);
			g_Ey = (real_t *)malloc(size);
			g_Ez = (real_t *)malloc(size);
			g_Hx = (real_t *)malloc(size);
			g_Hy = (real_t *)malloc(size);
			g_Hz = (real_t *)malloc(size);
			if ((g_Ex == NULL) || (g_Ey == NULL) || (g_Ez == NULL) ||
			    (g_Hx == NULL) || (g_Hy == NULL) || (g_Hz == NULL)) {
				fprintf(stderr, "*** snapshot buffer malloc error (NN=%zu)\n", (size_t)G_NN);
				fflush(stderr);
				G_NN = 0;
				return;
			}
			memset(g_Ex, 0, size);  memset(g_Ey, 0, size);  memset(g_Ez, 0, size);
			memset(g_Hx, 0, size);  memset(g_Hy, 0, size);  memset(g_Hz, 0, size);
		}
	}
	if (G_NN <= 0) return;

	// root : 自分の分をコピーする
	if (commRank == 0) {
		for (int i = iMin; i <  iMax; i++) {
		for (int j = jMin; j <= jMax; j++) {
		for (int k = kMin; k <= kMax; k++) {
			g_n = (G_Ni * i) + (G_Nj * j) + (G_Nk * k) + G_N0;
			g_Ex[g_n] = Ex[NA(i, j, k)];
		}
		}
		}
		for (int i = iMin; i <= iMax; i++) {
		for (int j = jMin; j <  jMax; j++) {
		for (int k = kMin; k <= kMax; k++) {
			g_n = (G_Ni * i) + (G_Nj * j) + (G_Nk * k) + G_N0;
			g_Ey[g_n] = Ey[NA(i, j, k)];
		}
		}
		}
		for (int i = iMin; i <= iMax; i++) {
		for (int j = jMin; j <= jMax; j++) {
		for (int k = kMin; k <  kMax; k++) {
			g_n = (G_Ni * i) + (G_Nj * j) + (G_Nk * k) + G_N0;
			g_Ez[g_n] = Ez[NA(i, j, k)];
		}
		}
		}
		for (int i = iMin - 0; i <= iMax; i++) {
		for (int j = jMin - 1; j <= jMax; j++) {
		for (int k = kMin - 1; k <= kMax; k++) {
			g_n = (G_Ni * i) + (G_Nj * j) + (G_Nk * k) + G_N0;
			g_Hx[g_n] = Hx[NA(i, j, k)];
		}
		}
		}
		for (int i = iMin - 1; i <= iMax; i++) {
		for (int j = jMin - 0; j <= jMax; j++) {
		for (int k = kMin - 1; k <= kMax; k++) {
			g_n = (G_Ni * i) + (G_Nj * j) + (G_Nk * k) + G_N0;
			g_Hy[g_n] = Hy[NA(i, j, k)];
		}
		}
		}
		for (int i = iMin - 1; i <= iMax; i++) {
		for (int j = jMin - 1; j <= jMax; j++) {
		for (int k = kMin - 0; k <= kMax; k++) {
			g_n = (G_Ni * i) + (G_Nj * j) + (G_Nk * k) + G_N0;
			g_Hz[g_n] = Hz[NA(i, j, k)];
		}
		}
		}
	}

	// root : 他ランクから受け取る
	if (commRank == 0) {
		for (int irank = 1; irank < commSize; irank++) {
			MPI_Recv(irecv, 11, MPI_INT, irank, tag, MPI_COMM_WORLD, &status);
			const int imin = irecv[0];
			const int imax = irecv[1];
			const int jmin = irecv[2];
			const int jmax = irecv[3];
			const int kmin = irecv[4];
			const int kmax = irecv[5];
			const int ni   = irecv[6];
			const int nj   = irecv[7];
			const int nk   = irecv[8];
			const int n0   = irecv[9];
			const int nn   = irecv[10];

			real_t *recv = (real_t *)malloc(nn * sizeof(real_t));
			if (recv == NULL) continue;

			// Ex
			MPI_Recv(recv, nn, MPI_REAL_T, irank, tag, MPI_COMM_WORLD, &status);
			for (int i = imin; i <  imax; i++) {
			for (int j = jmin; j <= jmax; j++) {
			for (int k = kmin; k <= kmax; k++) {
				g_n = (G_Ni * i) + (G_Nj * j) + (G_Nk * k) + G_N0;
				l_n = (ni * i) + (nj * j) + (nk * k) + n0;
				g_Ex[g_n] = recv[l_n];
			}
			}
			}

			// Ey
			MPI_Recv(recv, nn, MPI_REAL_T, irank, tag, MPI_COMM_WORLD, &status);
			for (int i = imin; i <= imax; i++) {
			for (int j = jmin; j <  jmax; j++) {
			for (int k = kmin; k <= kmax; k++) {
				g_n = (G_Ni * i) + (G_Nj * j) + (G_Nk * k) + G_N0;
				l_n = (ni * i) + (nj * j) + (nk * k) + n0;
				g_Ey[g_n] = recv[l_n];
			}
			}
			}

			// Ez
			MPI_Recv(recv, nn, MPI_REAL_T, irank, tag, MPI_COMM_WORLD, &status);
			for (int i = imin; i <= imax; i++) {
			for (int j = jmin; j <= jmax; j++) {
			for (int k = kmin; k <  kmax; k++) {
				g_n = (G_Ni * i) + (G_Nj * j) + (G_Nk * k) + G_N0;
				l_n = (ni * i) + (nj * j) + (nk * k) + n0;
				g_Ez[g_n] = recv[l_n];
			}
			}
			}

			// Hx
			MPI_Recv(recv, nn, MPI_REAL_T, irank, tag, MPI_COMM_WORLD, &status);
			for (int i = imin - 0; i <= imax; i++) {
			for (int j = jmin - 1; j <= jmax; j++) {
			for (int k = kmin - 1; k <= kmax; k++) {
				g_n = (G_Ni * i) + (G_Nj * j) + (G_Nk * k) + G_N0;
				l_n = (ni * i) + (nj * j) + (nk * k) + n0;
				g_Hx[g_n] = recv[l_n];
			}
			}
			}

			// Hy
			MPI_Recv(recv, nn, MPI_REAL_T, irank, tag, MPI_COMM_WORLD, &status);
			for (int i = imin - 1; i <= imax; i++) {
			for (int j = jmin - 0; j <= jmax; j++) {
			for (int k = kmin - 1; k <= kmax; k++) {
				g_n = (G_Ni * i) + (G_Nj * j) + (G_Nk * k) + G_N0;
				l_n = (ni * i) + (nj * j) + (nk * k) + n0;
				g_Hy[g_n] = recv[l_n];
			}
			}
			}

			// Hz
			MPI_Recv(recv, nn, MPI_REAL_T, irank, tag, MPI_COMM_WORLD, &status);
			for (int i = imin - 1; i <= imax; i++) {
			for (int j = jmin - 1; j <= jmax; j++) {
			for (int k = kmin - 0; k <= kmax; k++) {
				g_n = (G_Ni * i) + (G_Nj * j) + (G_Nk * k) + G_N0;
				l_n = (ni * i) + (nj * j) + (nk * k) + n0;
				g_Hz[g_n] = recv[l_n];
			}
			}
			}

			free(recv);
		}
	}
	// root 以外 : 自分の領域を送る
	else {
		isend[0] = iMin;
		isend[1] = iMax;
		isend[2] = jMin;
		isend[3] = jMax;
		isend[4] = kMin;
		isend[5] = kMax;
		isend[6] = (int)Ni;
		isend[7] = (int)Nj;
		isend[8] = (int)Nk;
		isend[9] = (int)N0;
		isend[10] = (int)NN;
		MPI_Send(isend, 11, MPI_INT, 0, tag, MPI_COMM_WORLD);

		MPI_Send(Ex, (int)NN, MPI_REAL_T, 0, tag, MPI_COMM_WORLD);
		MPI_Send(Ey, (int)NN, MPI_REAL_T, 0, tag, MPI_COMM_WORLD);
		MPI_Send(Ez, (int)NN, MPI_REAL_T, 0, tag, MPI_COMM_WORLD);
		MPI_Send(Hx, (int)NN, MPI_REAL_T, 0, tag, MPI_COMM_WORLD);
		MPI_Send(Hy, (int)NN, MPI_REAL_T, 0, tag, MPI_COMM_WORLD);
		MPI_Send(Hz, (int)NN, MPI_REAL_T, 0, tag, MPI_COMM_WORLD);
	}
#endif
}
