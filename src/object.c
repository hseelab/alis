/**********************************************************************
 * ALiS is Copyright (C) 2009-2019 by Ho-Seok Ee <hsee@kongju.ac.kr>. *
 * Redistribution and use with or without modification, are permitted *
 * under the terms of the Artistic License version 2.                 *
 **********************************************************************/

#include <string.h>
#include "alis_c.h"

#define FOR2D(i,iMin,iMax,j,jMin,jMax) \
_Pragma("omp parallel for schedule (guided)") for (int i=iMin; i<=iMax; i++) \
_Pragma("GCC ivdep") for (int j=jMin; j<=jMax; j++)

#define FOR3D(i,iMin,iMax,j,jMin,jMax,k,kMin,kMax) \
_Pragma("omp parallel for schedule (guided)") for (int i=iMin; i<=iMax; i++) for (int j=jMin; j<=jMax; j++) \
_Pragma("GCC ivdep") for (int k=kMin; k<=kMax; k++)


typedef struct {
	int iMin, iMax, jMin, jMax, kMin, kMax;
} bounds;

static float complex Gamma(world W, matter M, int p)
{
	return M.P[p].gamma + I * M.P[p].omega;
}

static float complex Sigma(world W, matter M, int p)
{
	return 2 * I * M.P[p].omega * (M.P[p].freal - I * M.P[p].fimag);
}

static float CJa(world W, matter M, int p)
{
	return (1 - PI * M.P[p].gamma * W->dt / W->eV) / (1 + PI * M.P[p].gamma * W->dt / W->eV);
}

static float CJb(world W, matter M, int p)
{
	return sq(2 * PI * M.P[p].omega * W->dt / W->eV) / (1 + PI * M.P[p].gamma * W->dt / W->eV);
}

static float complex CKa(world W, matter M, int p)
{
	return (1 - PI * Gamma(W, M, p) * W->dt / W->eV) / (1 + PI * Gamma(W, M, p) * W->dt / W->eV);
}

static float complex CKb(world W, matter M, int p)
{
	return - Gamma(W, M, p) * Sigma(W, M, p) * sq(2 * PI * W->dt / W->eV) / (1 + PI * Gamma(W, M, p) * W->dt / W->eV);
}

static float Ca(world W, matter M, int q)
{
	int m = M.e.y || M.e.z ? 1 : 0;
	float Conductivity=0;
	for (int p=0; p<MAXPOLES; p+=1+2*m) {
		if (!M.P[p+q*m].omega && M.P[p+q*m].gamma > 0) Conductivity += M.P[p+q*m].gamma;
		if (M.P[p+q*m].fimag) Conductivity += crealf(Sigma(W, M, p+q*m));
	}
	return Conductivity / W->eV;
}

static float Cb(world W, matter M, int q)
{
	return ((q==2 && M.e.z) ? M.e.z : (q==1 && M.e.y) ? M.e.y : M.e.x) + PI * W->dt * Ca(W, M, q);
}


static float inThisObject(world W, object *O, int o, float i, float j, float k)
{
	if (!(O[o].S(O[o], W->itox(W, i), W->jtoy(W, j), W->ktoz(W, k)))) return 0;
	while (o--) if (O[o].S(O[o], W->itox(W, i), W->jtoy(W, j), W->ktoz(W, k))) return 0;
	return 1;
}


static float fillFactor(world W, object *O, int o, float i, float j, float k)
{
	float n = inThisObject(W, O, o, i, j, k);
	if (W->subN < 2) {
		if (n || W->subN >= 0) {
			return n;
		}
		else if (floorf(i) < ceilf(i)) {
			if (inThisObject(W, O, o, i-0.5, j-0.5, k) && inThisObject(W, O, o, i+0.5, j-0.5, k)) return 1;
			if (inThisObject(W, O, o, i-0.5, j+0.5, k) && inThisObject(W, O, o, i+0.5, j+0.5, k)) return 1;
			if (inThisObject(W, O, o, i-0.5, j, k-0.5) && inThisObject(W, O, o, i+0.5, j, k-0.5)) return 1;
			if (inThisObject(W, O, o, i-0.5, j, k+0.5) && inThisObject(W, O, o, i+0.5, j, k+0.5)) return 1;
		}
		else if (floorf(j) < ceilf(j)) {
			if (inThisObject(W, O, o, i-0.5, j-0.5, k) && inThisObject(W, O, o, i-0.5, j+0.5, k)) return 1;
			if (inThisObject(W, O, o, i+0.5, j-0.5, k) && inThisObject(W, O, o, i+0.5, j+0.5, k)) return 1;
			if (inThisObject(W, O, o, i, j-0.5, k-0.5) && inThisObject(W, O, o, i, j+0.5, k-0.5)) return 1;
			if (inThisObject(W, O, o, i, j-0.5, k+0.5) && inThisObject(W, O, o, i, j+0.5, k+0.5)) return 1;
		}
		else if (floorf(k) < ceilf(k)) {
			if (inThisObject(W, O, o, i-0.5, j, k-0.5) && inThisObject(W, O, o, i-0.5, j, k+0.5)) return 1;
			if (inThisObject(W, O, o, i+0.5, j, k-0.5) && inThisObject(W, O, o, i+0.5, j, k+0.5)) return 1;
			if (inThisObject(W, O, o, i, j-0.5, k-0.5) && inThisObject(W, O, o, i, j-0.5, k+0.5)) return 1;
			if (inThisObject(W, O, o, i, j+0.5, k-0.5) && inThisObject(W, O, o, i, j+0.5, k+0.5)) return 1;
		}
		return 0;
	}

	float outmost = (W->subN-1.0) / (W->subN*2.0);
	float nnn = inThisObject(W, O, o, i-outmost, j-outmost, k-outmost);
	float nnp = inThisObject(W, O, o, i-outmost, j-outmost, k+outmost);
	float npn = inThisObject(W, O, o, i-outmost, j+outmost, k-outmost);
	float npp = inThisObject(W, O, o, i-outmost, j+outmost, k+outmost);
	float pnn = inThisObject(W, O, o, i+outmost, j-outmost, k-outmost);
	float pnp = inThisObject(W, O, o, i+outmost, j-outmost, k+outmost);
	float ppn = inThisObject(W, O, o, i+outmost, j+outmost, k-outmost);
	float ppp = inThisObject(W, O, o, i+outmost, j+outmost, k+outmost);
	if (n==nnn && nnn==nnp && nnp==npn && npn==npp && npp==pnn && pnn==pnp && pnp==ppn && ppn==ppp) return n;

	float Portion = 0;
	for (int di=1-W->subN; di<W->subN; di+=2) for (int dj=1-W->subN; dj<W->subN; dj+=2) for (int dk=1-W->subN; dk<W->subN; dk+=2)
		Portion += inThisObject(W, O, o, i+0.5*di/W->subN, j+0.5*dj/W->subN, k+0.5*dk/W->subN);
	return Portion / (W->subN * W->subN * W->subN);
}


static void makeCoeffs(world W, vfield *F, coeffs *C, float Cx, float Cy, float Cz)
{
	for (int m=0; m<=W->complexField; m++) {
		MALLOC(F[m].Jx, 0, C->N);
		MALLOC(F[m].Jy, 0, C->N);
		MALLOC(F[m].Jz, 0, C->N);
		MALLOC(F[m].Kx, 0, C->N);
		MALLOC(F[m].Ky, 0, C->N);
		MALLOC(F[m].Kz, 0, C->N);
	}

	MALLOC(C->NJ, 0, C->N);
	MALLOC(C->NK, 0, C->N);
	MALLOC(C->iMin, 0, C->N);
	MALLOC(C->iMax, 0, C->N);
	MALLOC(C->jMin, 0, C->N);
	MALLOC(C->jMax, 0, C->N);
	MALLOC(C->kMin, 0, C->N);
	MALLOC(C->kMax, 0, C->N);
	MALLOC(C->Jx, 0, C->N);
	MALLOC(C->Jy, 0, C->N);
	MALLOC(C->Jz, 0, C->N);
	MALLOC(C->ax, 0, C->N);
	MALLOC(C->ay, 0, C->N);
	MALLOC(C->az, 0, C->N);
	MALLOC(C->JJx, 0, C->N);
	MALLOC(C->JJy, 0, C->N);
	MALLOC(C->JJz, 0, C->N);
	MALLOC(C->JKx, 0, C->N);
	MALLOC(C->JKy, 0, C->N);
	MALLOC(C->JKz, 0, C->N);
	MALLOC(C->JEx, 0, C->N);
	MALLOC(C->JEy, 0, C->N);
	MALLOC(C->JEz, 0, C->N);
	MALLOC(C->KJx, 0, C->N);
	MALLOC(C->KJy, 0, C->N);
	MALLOC(C->KJz, 0, C->N);
	MALLOC(C->KKx, 0, C->N);
	MALLOC(C->KKy, 0, C->N);
	MALLOC(C->KKz, 0, C->N);
	MALLOC(C->KEx, 0, C->N);
	MALLOC(C->KEy, 0, C->N);
	MALLOC(C->KEz, 0, C->N);

	C->x = makeField(W->iMIN, W->iMAX+1, W->jMIN, W->jMAX, W->kMIN, W->kMAX);
	C->y = makeField(W->iMIN, W->iMAX, W->jMIN, W->jMAX+1, W->kMIN, W->kMAX);
	C->z = makeField(W->iMIN, W->iMAX, W->jMIN, W->jMAX, W->kMIN, W->kMAX+1);

	FOR3D(i,W->iMIN,W->iMAX+1, j,W->jMIN,W->jMAX, k,W->kMIN,W->kMAX) C->x[i][j][k] = 1/Cx;
	FOR3D(i,W->iMIN,W->iMAX, j,W->jMIN,W->jMAX+1, k,W->kMIN,W->kMAX) C->y[i][j][k] = 1/Cy;
	FOR3D(i,W->iMIN,W->iMAX, j,W->jMIN,W->jMAX, k,W->kMIN,W->kMAX+1) C->z[i][j][k] = 1/Cz;
}


static void makeAuxiliaryFields(world W, vfield *F, coeffs *C, int n, bounds B)
{
	C->iMin[n] = B.iMin, C->iMax[n] = B.iMax + (B.iMax == W->iMAX ? 1 : 0);
	C->jMin[n] = B.jMin, C->jMax[n] = B.jMax + (B.jMax == W->jMAX ? 1 : 0);
	C->kMin[n] = B.kMin, C->kMax[n] = B.kMax + (B.kMax == W->kMAX ? 1 : 0);

	MALLOC(C->JJx[n], 0, C->NK[n]);
	MALLOC(C->JJy[n], 0, C->NK[n]);
	MALLOC(C->JJz[n], 0, C->NK[n]);
	MALLOC(C->JKx[n], 0, C->NK[n]);
	MALLOC(C->JKy[n], 0, C->NK[n]);
	MALLOC(C->JKz[n], 0, C->NK[n]);
	MALLOC(C->JEx[n], 0, C->NK[n]);
	MALLOC(C->JEy[n], 0, C->NK[n]);
	MALLOC(C->JEz[n], 0, C->NK[n]);
	MALLOC(C->KJx[n], C->NJ[n], C->NK[n] - C->NJ[n]);
	MALLOC(C->KJy[n], C->NJ[n], C->NK[n] - C->NJ[n]);
	MALLOC(C->KJz[n], C->NJ[n], C->NK[n] - C->NJ[n]);
	MALLOC(C->KKx[n], C->NJ[n], C->NK[n] - C->NJ[n]);
	MALLOC(C->KKy[n], C->NJ[n], C->NK[n] - C->NJ[n]);
	MALLOC(C->KKz[n], C->NJ[n], C->NK[n] - C->NJ[n]);
	MALLOC(C->KEx[n], C->NJ[n], C->NK[n] - C->NJ[n]);
	MALLOC(C->KEy[n], C->NJ[n], C->NK[n] - C->NJ[n]);
	MALLOC(C->KEz[n], C->NJ[n], C->NK[n] - C->NJ[n]);

	C->Jx[n] = makeField(C->iMin[n], C->iMax[n], C->jMin[n], C->jMax[n], C->kMin[n], C->kMax[n]);
	C->Jy[n] = makeField(C->iMin[n], C->iMax[n], C->jMin[n], C->jMax[n], C->kMin[n], C->kMax[n]);
	C->Jz[n] = makeField(C->iMin[n], C->iMax[n], C->jMin[n], C->jMax[n], C->kMin[n], C->kMax[n]);

	for (int m=0; m<=W->complexField; m++) {
		MALLOC(F[m].Jx[n], 0, C->NK[n]);
		MALLOC(F[m].Jy[n], 0, C->NK[n]);
		MALLOC(F[m].Jz[n], 0, C->NK[n]);
		MALLOC(F[m].Kx[n], C->NJ[n], C->NK[n] - C->NJ[n]);
		MALLOC(F[m].Ky[n], C->NJ[n], C->NK[n] - C->NJ[n]);
		MALLOC(F[m].Kz[n], C->NJ[n], C->NK[n] - C->NJ[n]);

		for (int p=0; p<C->NK[n]; p++) {
			F[m].Jx[n][p] = makeField(C->iMin[n], C->iMax[n], C->jMin[n], C->jMax[n], C->kMin[n], C->kMax[n]);
			F[m].Jy[n][p] = makeField(C->iMin[n], C->iMax[n], C->jMin[n], C->jMax[n], C->kMin[n], C->kMax[n]);
			F[m].Jz[n][p] = makeField(C->iMin[n], C->iMax[n], C->jMin[n], C->jMax[n], C->kMin[n], C->kMax[n]);
		}
		for (int p=C->NJ[n]; p<C->NK[n]; p++) {
			F[m].Kx[n][p] = makeField(C->iMin[n], C->iMax[n], C->jMin[n], C->jMax[n], C->kMin[n], C->kMax[n]);
			F[m].Ky[n][p] = makeField(C->iMin[n], C->iMax[n], C->jMin[n], C->jMax[n], C->kMin[n], C->kMax[n]);
			F[m].Kz[n][p] = makeField(C->iMin[n], C->iMax[n], C->jMin[n], C->jMax[n], C->kMin[n], C->kMax[n]);
		}
	}
}


static bounds *getObjectsBounds(world W, object *O, int N)
{
	bounds *B;
	MALLOC(B, 0, N);

	for (int o=0; o<N; o++) {
		B[o].iMin = W->iMAX, B[o].iMax = W->iMIN;
		B[o].jMin = W->jMAX, B[o].jMax = W->jMIN;
		B[o].kMin = W->kMAX, B[o].kMax = W->kMIN;
	}
	#pragma omp parallel
	{
		bounds *b;
		MALLOC(b, 0, N);

		for (int o=0; o<N; o++) {
			b[o].iMin = W->iMAX, b[o].iMax = W->iMIN;
			b[o].jMin = W->jMAX, b[o].jMax = W->jMIN;
			b[o].kMin = W->kMAX, b[o].kMax = W->kMIN;
		}
		#pragma omp for
		for (int i=W->iMIN; i<=W->iMAX; i++) for (int j=W->jMIN; j<=W->jMAX; j++) for (int k=W->kMIN; k<=W->kMAX; k++) {
			int o = 0;
			while (O[o].S && !(O[o].S(O[o], W->itox(W, i), W->jtoy(W, j), W->ktoz(W, k)))) o++;
			if (o < N) {
				b[o].iMin = MIN(i, b[o].iMin), b[o].iMax = MAX(i, b[o].iMax);
				b[o].jMin = MIN(j, b[o].jMin), b[o].jMax = MAX(j, b[o].jMax);
				b[o].kMin = MIN(k, b[o].kMin), b[o].kMax = MAX(k, b[o].kMax);
			}
		}
		#pragma omp critical
		for (int o=0; o<N; o++) {
			B[o].iMin = MIN(b[o].iMin, B[o].iMin), B[o].iMax = MAX(b[o].iMax, B[o].iMax);
			B[o].jMin = MIN(b[o].jMin, B[o].jMin), B[o].jMax = MAX(b[o].jMax, B[o].jMax);
			B[o].kMin = MIN(b[o].kMin, B[o].kMin), B[o].kMax = MAX(b[o].kMax, B[o].kMax);
		}
		free(b);
	}
	for (int o=0; o<N; o++) {
		if (B[o].iMin <= B[o].iMax && B[o].jMin <= B[o].jMax && B[o].kMin <= B[o].kMax) {
			B[o].iMin = MAX(B[o].iMin-1, W->iMIN), B[o].iMax = MIN(B[o].iMax+2, W->iMAX);
			B[o].jMin = MAX(B[o].jMin-1, W->jMIN), B[o].jMax = MIN(B[o].jMax+2, W->jMAX);
			B[o].kMin = MAX(B[o].kMin-1, W->kMIN), B[o].kMax = MIN(B[o].kMax+2, W->kMAX);
		}
		else {
			printf("\rVolume of object %d/%d is zero!\n", o+1, N);
			exit(0);
		}
	}
	return B;
}


static int auxiliaryFieldRequired(matter M)
{
	for (int i=0; i<3; i++) for (int j=0; j<6; j++) if (M.d[i][j]) return 2;
	if (M.P[0].omega || M.P[0].gamma) return 1;
	return 0;
}


static void setAuxiliaryCoeffs(world W, coeffs *C, matter *M, object *O, int N, int n, int o, int m)
{
	FOR3D(i,C->iMin[n],C->iMax[n],j,C->jMin[n],C->jMax[n],k,C->kMin[n],C->kMax[n]) {
		C->Jx[n][i][j][k] = fillFactor(W, O, o, i-0.5, j, k);
		C->Jy[n][i][j][k] = fillFactor(W, O, o, i, j-0.5, k);
		C->Jz[n][i][j][k] = fillFactor(W, O, o, i, j, k-0.5);
	}

	C->ax[n] = PI * W->dt * Ca(W, M[o], 0) / (M[o].e.x ? M[o].e.x : M[o].e.x);
	C->ay[n] = PI * W->dt * Ca(W, M[o], 1) / (M[o].e.y ? M[o].e.y : M[o].e.x);
	C->az[n] = PI * W->dt * Ca(W, M[o], 2) / (M[o].e.z ? M[o].e.z : M[o].e.x);

	for (int p=0, j=0, k=C->NJ[n]; p<MAXPOLES; p+=1+2*m) {
		if (M[o].P[p].omega) {
			if (!M[o].P[p].freal && !M[o].P[p].fimag) {
				C->JJx[n][j] = CJa(W, M[o], p+0*m);
				C->JJy[n][j] = CJa(W, M[o], p+1*m);
				C->JJz[n][j] = CJa(W, M[o], p+2*m);
				C->JEx[n][j] = CJb(W, M[o], p+0*m) / Cb(W, M[o], 0);
				C->JEy[n][j] = CJb(W, M[o], p+1*m) / Cb(W, M[o], 1);
				C->JEz[n][j] = CJb(W, M[o], p+2*m) / Cb(W, M[o], 2);
				j++;
			}
			else {
				C->JJx[n][k] = creal(CKa(W, M[o], p+0*m));
				C->JJy[n][k] = creal(CKa(W, M[o], p+1*m));
				C->JJz[n][k] = creal(CKa(W, M[o], p+2*m));
				C->JKx[n][k] =-cimag(CKa(W, M[o], p+0*m));
				C->JKy[n][k] =-cimag(CKa(W, M[o], p+1*m));
				C->JKz[n][k] =-cimag(CKa(W, M[o], p+2*m));
				C->JEx[n][k] = creal(CKb(W, M[o], p+0*m)) / Cb(W, M[o], 0);
				C->JEy[n][k] = creal(CKb(W, M[o], p+1*m)) / Cb(W, M[o], 1);
				C->JEz[n][k] = creal(CKb(W, M[o], p+2*m)) / Cb(W, M[o], 2);
				C->KJx[n][k] =-C->JKx[n][k] / C->JJx[n][k];
				C->KJy[n][k] =-C->JKy[n][k] / C->JJy[n][k];
				C->KJz[n][k] =-C->JKz[n][k] / C->JJz[n][k];
				C->KKx[n][k] = C->JKx[n][k] * C->JKx[n][k] / C->JJx[n][k] + C->JJx[n][k];
				C->KKy[n][k] = C->JKy[n][k] * C->JKy[n][k] / C->JJy[n][k] + C->JJy[n][k];
				C->KKz[n][k] = C->JKz[n][k] * C->JKz[n][k] / C->JJz[n][k] + C->JJz[n][k];
				C->KEx[n][k] = C->JEx[n][k] * C->JKx[n][k] / C->JJx[n][k] + cimag(CKb(W, M[o], p+0*m)) / Cb(W, M[o], 0);
				C->KEy[n][k] = C->JEy[n][k] * C->JKy[n][k] / C->JJy[n][k] + cimag(CKb(W, M[o], p+1*m)) / Cb(W, M[o], 1);
				C->KEz[n][k] = C->JEz[n][k] * C->JKz[n][k] / C->JJz[n][k] + cimag(CKb(W, M[o], p+2*m)) / Cb(W, M[o], 2);
				k++;
			}
		}
	}
}


void putObjectArray(world W, matter *M, object *O)
{
	removeObjects(W);
	int N;
	for (N=0; O[N].S; N++) if (auxiliaryFieldRequired(M[N])) W->CE->N++;
	printf("\rPreparing to put %d objects.\n", N);
	W->subN = W->Res.spec.subN = W->Res.spec.subN ? W->Res.spec.subN : W->CE->N ? 1 : 12;
	bounds *B = getObjectsBounds(W, O, N);
	makeCoeffs(W, W->E, W->CE, Cb(W, M[N], 0), Cb(W, M[N], 1), Cb(W, M[N], 2));
	MALLOC(W->CE->O, 0, N+1);
	MALLOC(W->CE->M, 0, N+1);
	W->CE->M[N] = M[N];

	for (int o=N-1, n=W->CE->N; o>=0; o--) {
		W->CE->O[o] = O[o];
		W->CE->M[o] = M[o];
		float cx = Cb(W, M[o], 0) - Cb(W, M[N], 0);
		float cy = Cb(W, M[o], 1) - Cb(W, M[N], 1);
		float cz = Cb(W, M[o], 2) - Cb(W, M[N], 2);

		printf("\rPutting objects: %d left.", o+1);
		fflush(stdout);

		if (auxiliaryFieldRequired(M[o])) {
			n--;
			int m = M[o].e.y || M[o].e.z ? 1 : 0;

			for (int p=0; p<MAXPOLES; p++) {
				if (M[o].P[p].omega && !M[o].P[p].freal && !M[o].P[p].fimag) W->CE->NJ[n]++;
				if (M[o].P[p].omega) W->CE->NK[n]++;
			}
			if (m) W->CE->NJ[n] /= 3, W->CE->NK[n] /= 3;
			makeAuxiliaryFields(W, W->E, W->CE, n, B[o]);
			setAuxiliaryCoeffs(W, W->CE, M, O, N, n, o, m);

			FOR3D(i,B[o].iMin,B[o].iMax,j,B[o].jMin,B[o].jMax,k,B[o].kMin,B[o].kMax) {
				W->CE->x[i][j][k] = 1 / (1 / W->CE->x[i][j][k] + cx * W->CE->Jx[n][i][j][k]);
				W->CE->y[i][j][k] = 1 / (1 / W->CE->y[i][j][k] + cy * W->CE->Jy[n][i][j][k]);
				W->CE->z[i][j][k] = 1 / (1 / W->CE->z[i][j][k] + cz * W->CE->Jz[n][i][j][k]);
			}
			if (B[o].iMax == W->iMAX) FOR2D(j,B[o].jMin,B[o].jMax,k,B[o].kMin,B[o].kMax)
				W->CE->x[B[o].iMax+1][j][k] = 1 / (1 / W->CE->x[B[o].iMax+1][j][k] + cx * W->CE->Jx[n][B[o].iMax+1][j][k]);
			if (B[o].jMax == W->jMAX) FOR2D(i,B[o].iMin,B[o].iMax,k,B[o].kMin,B[o].kMax)
				W->CE->y[i][B[o].jMax+1][k] = 1 / (1 / W->CE->y[i][B[o].jMax+1][k] + cy * W->CE->Jy[n][i][B[o].jMax+1][k]);
			if (B[o].kMax == W->kMAX) FOR2D(i,B[o].iMin,B[o].iMax,j,B[o].jMin,B[o].jMax)
				W->CE->z[i][j][B[o].kMax+1] = 1 / (1 / W->CE->z[i][j][B[o].kMax+1] + cz * W->CE->Jz[n][i][j][B[o].kMax+1]);
		}

		else {
			FOR3D(i,B[o].iMin,B[o].iMax,j,B[o].jMin,B[o].jMax,k,B[o].kMin,B[o].kMax) {
				W->CE->x[i][j][k] = 1 / (1 / W->CE->x[i][j][k] + cx * fillFactor(W, O, o, i-0.5, j, k));
				W->CE->y[i][j][k] = 1 / (1 / W->CE->y[i][j][k] + cy * fillFactor(W, O, o, i, j-0.5, k));
				W->CE->z[i][j][k] = 1 / (1 / W->CE->z[i][j][k] + cz * fillFactor(W, O, o, i, j, k-0.5));
			}
			if (B[o].iMax == W->iMAX) FOR2D(j,B[o].jMin,B[o].jMax,k,B[o].kMin,B[o].kMax)
				W->CE->x[B[o].iMax+1][j][k] = 1 / (1 / W->CE->x[B[o].iMax+1][j][k] + cx * fillFactor(W, O, o, B[o].iMax+0.5, j, k));
			if (B[o].jMax == W->jMAX) FOR2D(i,B[o].iMin,B[o].iMax,k,B[o].kMin,B[o].kMax)
				W->CE->y[i][B[o].jMax+1][k] = 1 / (1 / W->CE->y[i][B[o].jMax+1][k] + cy * fillFactor(W, O, o, i, B[o].jMax+0.5, k));
			if (B[o].kMax == W->kMAX) FOR2D(i,B[o].iMin,B[o].iMax,j,B[o].jMin,B[o].jMax)
				W->CE->z[i][j][B[o].kMax+1] = 1 / (1 / W->CE->z[i][j][B[o].kMax+1] + cz * fillFactor(W, O, o, i, j, B[o].kMax+0.5));
		}
	}

	free(B);
	printf("\rPutting objects: %d done.\n", N);
	fflush(stdout);
}


void (putObjects)(world W, ...)
{
	int N=0;
	va_list ap;
	matter *M;
	object *O;

	va_start(ap, W);
	while (va_arg(ap, matter), va_arg(ap, object).S) N++;
	va_end(ap);
	MALLOC(M, 0, N+1);
	MALLOC(O, 0, N+1);
	va_start(ap, W);
	for (int o=0; o<=N; o++) {
		M[o] = va_arg(ap, matter);
		O[o] = va_arg(ap, object);
	}
	va_end(ap);
	putObjectArray(W, M, O);
	free(M);
	free(O);
}


static void removeAuxiliaryFields(world W, vfield *F, coeffs *C)
{
	for (int m=0; m<=W->complexField; m++) {
		for (int n=0; n<C->N; n++) {
			for (int p=0; p<C->NK[n]; p++) {
				removeField(F[m].Jx[n][p], C->iMin[n], C->jMin[n], C->kMin[n]);
				removeField(F[m].Jy[n][p], C->iMin[n], C->jMin[n], C->kMin[n]);
				removeField(F[m].Jz[n][p], C->iMin[n], C->jMin[n], C->kMin[n]);
			}
			for (int p=C->NJ[n]; p<C->NK[n]; p++) {
				removeField(F[m].Kx[n][p], C->iMin[n], C->jMin[n], C->kMin[n]);
				removeField(F[m].Ky[n][p], C->iMin[n], C->jMin[n], C->kMin[n]);
				removeField(F[m].Kz[n][p], C->iMin[n], C->jMin[n], C->kMin[n]);
			}
			free(F[m].Jx[n]); free(F[m].Kx[n]+C->NJ[n]);
			free(F[m].Jy[n]); free(F[m].Ky[n]+C->NJ[n]);
			free(F[m].Jz[n]); free(F[m].Kz[n]+C->NJ[n]);
		}
		free(F[m].Jx); free(F[m].Kx);
		free(F[m].Jy); free(F[m].Ky);
		free(F[m].Jz); free(F[m].Kz);
	}
}


static void removeCoeffs(world W, coeffs *C) {
	for (int n=0; n<C->N; n++) {
		removeField(C->Jx[n], C->iMin[n], C->jMin[n], C->kMin[n]);
		removeField(C->Jy[n], C->iMin[n], C->jMin[n], C->kMin[n]);
		removeField(C->Jz[n], C->iMin[n], C->jMin[n], C->kMin[n]);
		free(C->JJx[n]); free(C->JKx[n]); free(C->JEx[n]);
		free(C->JJy[n]); free(C->JKy[n]); free(C->JEy[n]);
		free(C->JJz[n]); free(C->JKz[n]); free(C->JEz[n]);
		free(C->KJx[n]+C->NJ[n]); free(C->KKx[n]+C->NJ[n]); free(C->KEx[n]+C->NJ[n]);
		free(C->KJy[n]+C->NJ[n]); free(C->KKy[n]+C->NJ[n]); free(C->KEy[n]+C->NJ[n]);
		free(C->KJz[n]+C->NJ[n]); free(C->KKz[n]+C->NJ[n]); free(C->KEz[n]+C->NJ[n]);
	}
	removeField(C->x, W->iMIN, W->jMIN, W->kMIN);
	removeField(C->y, W->iMIN, W->jMIN, W->kMIN);
	removeField(C->z, W->iMIN, W->jMIN, W->kMIN);
	free(C->NJ);
	free(C->NK);
	free(C->iMin); free(C->jMax);
	free(C->iMax); free(C->kMin);
	free(C->jMin); free(C->kMax);
	free(C->Jx); free(C->ax);
	free(C->Jy); free(C->ay);
	free(C->Jz); free(C->az);
	free(C->JJx); free(C->JKx); free(C->JEx); free(C->KJx); free(C->KKx); free(C->KEx);
	free(C->JJy); free(C->JKy); free(C->JEy); free(C->KJy); free(C->KKy); free(C->KEy);
	free(C->JJz); free(C->JKz); free(C->JEz); free(C->KJz); free(C->KKz); free(C->KEz);
	free(C->O);
	free(C->M);
	C->O = 0;
	C->M = 0;
	C->N = 0;
}


void removeObjects(world W)
{
	if (W->CE->O && W->CE->M) {
		removeAuxiliaryFields(W, W->E, W->CE);
		removeCoeffs(W, W->CE);
	}
}


matter LT(matter M, float T)
{
	for(int i=0; i<MAXPOLES; i++) {
		M.P[i].gamma *= T;
		M.P[i].fimag *= T;
	}
	return M;
}


matter nkFile(char *name, float l, float f)
{
	int N=0;
	char str[1024];
	FILE *file;
	float n=0, k=0, *lArr=0, *nArr=0, *kArr=0;

	if (file=fopen(name, "r")) {
		for (int i=0; fgets(str, 256, file); N++) {
			if (i <= N) {
				i += 256;
				lArr = realloc(lArr, i*sizeof(float));
				nArr = realloc(nArr, i*sizeof(float));
				kArr = realloc(kArr, i*sizeof(float));
			}
			char *s=str;
			lArr[N] = atof(s);
			s = strchr(s,'\t') ? strchr(s,'\t')+1 : strchr(s,' ') ? strchr(s,' ')+1 : strchr(s,',')+1;
			nArr[N] = atof(s);                                                                                 
			s = strchr(s,'\t') ? strchr(s,'\t')+1 : strchr(s,' ') ? strchr(s,' ')+1 : strchr(s,',')+1;
			kArr[N] = atof(s);
		}

		int i=1;
		while (i<N-1 && l>lArr[i]) i++;
		n = (nArr[i-1]*(lArr[i]-l) + nArr[i]*(l-lArr[i-1])) / (lArr[i]-lArr[i-1]);
		k = (kArr[i-1]*(lArr[i]-l) + kArr[i]*(l-lArr[i-1])) / (lArr[i]-lArr[i-1]);

		free(lArr);
		free(nArr);
		free(kArr);
	}

	else {
		printf("The (n,k) file not found.\n");
		exit(0);
	}

	return (matter) {{(n)*(n)-(k)*(k)},{0,2*n*k*((float)f)}};
}
