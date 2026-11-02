/*
setupDispersion.c

setup dispersion data
*/

#include "ofd.h"

static void setupDispersion_id(int);


/*
RC (recursive convolution) factors of dispersion material (type = 2)

eps(w) = einf + ae / (jw) + be / (ce + jw)
chi(t) = ae + be * exp(-ce * t)  (t >= 0)

xi0  = integral of chi(t) over [0, dt]
dxi0 = xi(0) - xi(1) (be term only, the ae term is constant)
ke   = exp(-ce * dt) = dxi(m + 1) / dxi(m)

ce = 0 is the limit ce -> 0 (be / ce * (1 - ke) -> be * dt, dxi0 -> 0) :
"material = 2 einf ae 0 0" (conductivity only) must not give 0 / 0.
*/
void dispersionFactor(double ae, double be, double ce, double dt, double *xi0, double *dxi0, double *ke)
{
	const double k = exp(-ce * dt);
	const double b = (ce != 0) ? (be / ce) * (1 - k) : be * dt;  // integral of be * exp(-ce * t)
	*xi0 = (ae * dt) + b;
	*dxi0 = (ce != 0) ? b * (1 - k) : 0;
	*ke = k;
}


void setupDispersion(void)
{
	for (int m = 0; m < NMaterial; m++) {
		if (Material[m].type == 2) {
			const double einf = Material[m].einf;
			double xi0, dxi0, ke;
			dispersionFactor(Material[m].ae, Material[m].be, Material[m].ce, Dt, &xi0, &dxi0, &ke);
			Material[m].edisp[0] = (real_t)(1 / (einf + xi0));
			Material[m].edisp[1] = (real_t)dxi0;
			Material[m].edisp[2] = (real_t)ke;
		}
	}

	// setup dispersion id

	numDispersionEx = numDispersionEy = numDispersionEz = 0;
	setupDispersion_id(0);
	//printf("%zd %zd %zd\n", numDispersionEx, numDispersionEy, numDispersionEz);

	if (numDispersionEx > 0) {
		mDispersionEx = (dispersion_t *)malloc(numDispersionEx * sizeof(dispersion_t));
		DispersionEx =        (real_t *)malloc(numDispersionEx * sizeof(real_t));
		memset(DispersionEx, 0,                numDispersionEx * sizeof(real_t));
	}
	if (numDispersionEy > 0) {
		mDispersionEy = (dispersion_t *)malloc(numDispersionEy * sizeof(dispersion_t));
		DispersionEy =        (real_t *)malloc(numDispersionEy * sizeof(real_t));
		memset(DispersionEy, 0,                numDispersionEy * sizeof(real_t));
	}
	if (numDispersionEz > 0) {
		mDispersionEz = (dispersion_t *)malloc(numDispersionEz * sizeof(dispersion_t));
		DispersionEz =        (real_t *)malloc(numDispersionEz * sizeof(real_t));
		memset(DispersionEz, 0,                numDispersionEz * sizeof(real_t));
	}

	setupDispersion_id(1);
}


// mode = 0/1  (OpenMP : NG for CUDA)
static void setupDispersion_id(int mode)
{
	// Ex
	int64_t nx = 0;
	for (int i = iMin; i <  iMax; i++) {
	for (int j = jMin; j <= jMax; j++) {
	for (int k = kMin; k <= kMax; k++) {
		const id_t m = IEX(i, j, k);
		if (Material[m].type == 2) {
			if (mode == 1) {
				mDispersionEx[nx].i = i;
				mDispersionEx[nx].j = j;
				mDispersionEx[nx].k = k;
				mDispersionEx[nx].f1 = Material[m].edisp[0];
				mDispersionEx[nx].f2 = Material[m].edisp[1];
				mDispersionEx[nx].f3 = Material[m].edisp[2];
			}
			nx++;
		}
	}
	}
	}
	if (mode == 0) {
		numDispersionEx = nx;
	}

	// Ey
	int64_t ny = 0;
	for (int i = iMin; i <= iMax; i++) {
	for (int j = jMin; j <  jMax; j++) {
	for (int k = kMin; k <= kMax; k++) {
		const id_t m = IEY(i, j, k);
		if (Material[m].type == 2) {
			if (mode == 1) {
				mDispersionEy[ny].i = i;
				mDispersionEy[ny].j = j;
				mDispersionEy[ny].k = k;
				mDispersionEy[ny].f1 = Material[m].edisp[0];
				mDispersionEy[ny].f2 = Material[m].edisp[1];
				mDispersionEy[ny].f3 = Material[m].edisp[2];
			}
			ny++;
		}
	}
	}
	}
	if (mode == 0) {
		numDispersionEy = ny;
	}

	// Ez
	int64_t nz = 0;
	for (int i = iMin; i <= iMax; i++) {
	for (int j = jMin; j <= jMax; j++) {
	for (int k = kMin; k <  kMax; k++) {
		const id_t m = IEZ(i, j, k);
		if (Material[m].type == 2) {
			if (mode == 1) {
				mDispersionEz[nz].i = i;
				mDispersionEz[nz].j = j;
				mDispersionEz[nz].k = k;
				mDispersionEz[nz].f1 = Material[m].edisp[0];
				mDispersionEz[nz].f2 = Material[m].edisp[1];
				mDispersionEz[nz].f3 = Material[m].edisp[2];
			}
			nz++;
		}
	}
	}
	}
	if (mode == 0) {
		numDispersionEz = nz;
	}
}
