/**********************************************************************
 * ALiS is Copyright (C) 2009-2019 by Ho-Seok Ee <hsee@kongju.ac.kr>. *
 * Redistribution and use with or without modification, are permitted *
 * under the terms of the Artistic License version 2.                 *
 **********************************************************************/

#include <sys/stat.h>
#include <string.h>
#include "alis_c.h"
#define FOR1D(i,iMin,iMax) \
_Pragma("omp parallel for") for (int i=P->iMin; i<=P->iMax; i++)


void CylindericalMap(int x, int y, int X, int Y, float angle[2])
{
	angle[0] = 180.0 - 180.0 * y / Y;
	angle[1] = 360.0 * x / X;
}


void NorthernAzimuthalMap(int x, int y, int X, int Y, float angle[2])
{
	angle[0] = 180 * sqrtf(sq(x/(X-1.0)-0.5) + sq(y/(Y-1.0)-0.5));
	angle[1] = 180 * atan2f(y/(Y-1.0)-0.5, x/(X-1.0)-0.5) / PI;
	if (angle[0] > 90) angle[0] = INF;
}


void SouthernAzimuthalMap(int x, int y, int X, int Y, float angle[2])
{
	angle[0] = 180 - 180 * sqrtf(sq(x/(X-1.0)-0.5) + sq(y/(Y-1.0)-0.5));
	angle[1] = 180 * atan2f(y/(Y-1.0)-0.5, x/(X-1.0)-0.5) / PI;
	if (angle[0] < 90) angle[0] = INF;
}


void AzimuthalMap(int x, int y, int X, int Y, float angle[2])
{
	if (x < (X+1)/2) NorthernAzimuthalMap(x, y, (X+1)/2, Y, angle);
	if (x > X/2) SouthernAzimuthalMap(x-X/2, y, (X+1)/2, Y, angle);
}


static float rotateA(float theta, float a, float b)
{
	return acosf(cosf(a*PI/180)*cosf(theta*PI/180)-(sinf(a*PI/180)*cosf(b*PI/180))*sinf(theta*PI/180))*180/PI;
}

static float rotateB(float theta, float a, float b)
{
	return atan2f(sinf(a*PI/180)*sinf(b*PI/180), (sinf(a*PI/180)*cosf(b*PI/180))*cosf(theta*PI/180)+cosf(a*PI/180)*sinf(theta*PI/180))*180/PI;
}


static int auxiliaryFieldRequired(matter M)
{
	for (int i=0; i<3; i++) for (int j=0; j<6; j++) if (M.d[i][j]) return 2;
	if (M.P[0].omega || M.P[0].gamma) return 1;
	return 0;
}


static float complex pole(world W, matter M, float omega, int p)
{
	if (M.P[p].freal || M.P[p].fimag)
		return M.P[p].omega
			*((M.P[p].freal + I * M.P[p].fimag) / (M.P[p].omega + omega + I * M.P[p].gamma)
			+ (M.P[p].freal - I * M.P[p].fimag) / (M.P[p].omega - omega - I * M.P[p].gamma));
	if (M.P[p].omega) return - sq(M.P[p].omega) / (sq(omega) + I * omega * M.P[p].gamma);
	return -I * M.P[p].gamma / omega;
}


void setFarfieldCoeffs(phaser P, world W) {
	P->kMIN = (W->zMinSur==SYM) ? -P->kMax : P->kMin;
	MALLOC(P->dE, P->kMIN, 2-P->kMIN+P->kMax);
	MALLOC(P->dH, P->kMIN, 2-P->kMIN+P->kMax);
	MALLOC(P->ex, P->kMIN, 2-P->kMIN+P->kMax);
	MALLOC(P->ez, P->kMIN, 2-P->kMIN+P->kMax);

	#pragma omp parallel for
	for (int k=P->kMin; k<=P->kMax+1; k++) {
		P->dE[k] = 1 / W->CE->dtdz[k];
		P->dH[k] = 1 / W->CH->dtdz[k];
		P->ex[k] = W->CE->x ? 1/W->CE->x[P->iMax][P->jMax][k] : 1;
		P->ez[k] = W->CE->z ? 1/W->CE->z[P->iMax][P->jMax][k] : 1;

		for (int n=0, o=0; n<W->CE->N; n++, o++) {
			while (!auxiliaryFieldRequired(W->CE->M[o])) o++;

			if (W->CE->iMin[n]<=P->iMax && P->iMax<=W->CE->iMax[n] && W->CE->jMin[n]<=P->jMax && P->jMax<=W->CE->jMax[n]) {
				if (W->CE->kMin[n]<=k && k<=W->CE->kMax[n]) {
					int m = W->CE->M[o].e.y || W->CE->M[o].e.z ? 1 : 0;

					P->ex[k] += -W->CE->ax[n] * (W->CE->M[o].e.x ? W->CE->M[o].e.x : W->CE->M[o].e.x) * W->CE->Jx[n][P->iMax][P->jMax][k];
					P->ez[k] += -W->CE->az[n] * (W->CE->M[o].e.z ? W->CE->M[o].e.z : W->CE->M[o].e.x) * W->CE->Jz[n][P->iMax][P->jMax][k];

					for (int p=0; p<MAXPOLES; p+=1+2*m) {
						P->ex[k] += pole(W, W->CE->M[o], W->eV * P->f, p+0*m) * W->CE->Jx[n][P->iMax][P->jMax][k];
						P->ez[k] += pole(W, W->CE->M[o], W->eV * P->f, p+2*m) * W->CE->Jz[n][P->iMax][P->jMax][k];
					}
				}
			}
		}
	}

	if (W->zMinSur==SYM) {
		#pragma omp parallel for
//		#pragma GCC ivdep
		for (int k=P->kMIN; k<P->kMin; k++) {
			P->dE[k] = P->dE[-k];
			P->dH[k] = P->dH[1-k];
			P->ex[k] = P->ex[-k];
			P->ez[k] = P->ez[1-k];
		}
	}
}


tfield Unpol(phaser P, float theta, float phi)
{
	return (tfield) {0};
}


tfield Ppol(phaser P, float theta, float phi)
{
	tfield T = {0};
	float sinth, sin2th, sinphi = sinf(phi*PI/180), cosphi = cosf(phi*PI/180);

	MALLOC(T.Ex, P->kMIN, 2-P->kMIN+P->kMax);
	MALLOC(T.Hy, P->kMIN, 2-P->kMIN+P->kMax);

	if (theta < 90 || 270 < theta) {
		sinth = sinf(theta*PI/180) * sqrtf(crealf(P->ex[P->kMax])), sin2th = sq(sinth);
		T.Hy[P->kMIN] = -1;
		T.Ex[P->kMIN] = csqrtf((1-sin2th/P->ez[P->kMIN]) / P->ex[P->kMIN]);
		T.Ex[P->kMIN] *= cexpf(-I * casinf(sinf(PI*P->f*P->W->dt) * csqrtf(P->ex[P->kMIN]-sin2th) * P->dE[P->kMIN]));

		for (int k=P->kMIN; k<=P->kMax; k++) {
			T.Hy[k+1] = T.Hy[k] + 2 * I * sinf(PI*P->f*P->W->dt) * P->dE[k] * T.Ex[k] * P->ex[k];
			T.Ex[k+1] = T.Ex[k] + 2 * I * sinf(PI*P->f*P->W->dt) * P->dH[k+1] * T.Hy[k+1] * (1-sin2th/P->ez[k+1]);
			if (cabsf(T.Hy[k+1])>1E10) for (int kk=P->kMIN; kk<=k+1; kk++) T.Ex[kk]/=1E10, T.Hy[kk]/=1E10;
		}
		T.A = sqrtf(sqrtf(P->ex[P->kMax]));
		T.A /= cabsf(0.5*(T.Hy[P->kMax]+T.Hy[P->kMax+1]) / sqrtf(P->ex[P->kMax]) - T.Ex[P->kMax] / cosf(theta*PI/180));
	}

	if (90 < theta && theta < 270) {
		sinth = sinf(theta*PI/180) * sqrtf(crealf(P->ex[P->kMIN])), sin2th = sq(sinth);
		T.Hy[P->kMax+1] = 1;
		T.Ex[P->kMax+1] = csqrtf((1-sin2th/P->ez[P->kMax+1]) / P->ex[P->kMax+1]);
		T.Ex[P->kMax+1] *= cexpf(I * casinf(sinf(PI*P->f*P->W->dt) * csqrtf(P->ex[P->kMax+1]-sin2th) * P->dE[P->kMax+1]));

		for (int k=P->kMax; k>=P->kMIN; k--) {
			T.Ex[k] = T.Ex[k+1] - 2 * I * sinf(PI*P->f*P->W->dt) * P->dH[k+1] * T.Hy[k+1] * (1-sin2th/P->ez[k+1]);
			T.Hy[k] = T.Hy[k+1] - 2 * I * sinf(PI*P->f*P->W->dt) * P->dE[k] * T.Ex[k] * P->ex[k];
			if (cabsf(T.Hy[k+1])>1E10) for (int kk=P->kMax; kk>=k; kk--) T.Ex[kk]/=1E10, T.Hy[kk]/=1E10;
		}
		T.A = sqrtf(sqrtf(P->ex[P->kMIN]));
		T.A /= cabsf(0.5*(T.Hy[P->kMIN]+T.Hy[P->kMIN+1]) / sqrtf(P->ex[P->kMIN]) - T.Ex[P->kMIN] / cosf(theta*PI/180));
	}

	if (theta) {
	   	MALLOC(T.Ez, P->kMIN, 2-P->kMIN+P->kMax);
		#pragma omp parallel for
//		#pragma GCC ivdep
		for (int k=P->kMIN; k<=P->kMax+1; k++) T.Ez[k] = T.Hy[k] * sinth/P->ez[k];
	}

	if (phi) {
		MALLOC(T.Hx, P->kMIN, 2-P->kMIN+P->kMax);
		MALLOC(T.Ey, P->kMIN, 2-P->kMIN+P->kMax);
		#pragma omp parallel for
//		#pragma GCC ivdep
		for (int k=P->kMIN; k<=P->kMax+1; k++) {
			T.Ey[k] = T.Ex[k] * sinphi;
			T.Ex[k] = T.Ex[k] * cosphi;
			T.Hx[k] =-T.Hy[k] * sinphi;
			T.Hy[k] = T.Hy[k] * cosphi;
		}
	}

	return T;
}


tfield Spol(phaser P, float theta, float phi)
{
	tfield T = {0};
	float sinth, sin2th, sinphi = sinf(phi*PI/180), cosphi = cosf(phi*PI/180);

	MALLOC(T.Hx, P->kMIN, 2-P->kMIN+P->kMax);
	MALLOC(T.Ey, P->kMIN, 2-P->kMIN+P->kMax);

	if (theta < 90 || 270 < theta) {
		sinth = sinf(theta*PI/180) * sqrtf(crealf(P->ex[P->kMax])), sin2th = sq(sinth);
		T.Hx[P->kMIN] = csqrtf(P->ex[P->kMIN] - sin2th);
		T.Ey[P->kMIN] = cexpf(-I * casinf(sinf(PI*P->f*P->W->dt) * csqrtf(P->ex[P->kMIN]-sin2th) * P->dE[P->kMIN]));

		for (int k=P->kMIN; k<=P->kMax; k++) {
			T.Hx[k+1] = T.Hx[k] - 2 * I * sinf(PI*P->f*P->W->dt) * P->dE[k] * T.Ey[k] * (P->ex[k] - sin2th);
			T.Ey[k+1] = T.Ey[k] - 2 * I * sinf(PI*P->f*P->W->dt) * P->dH[k+1] * T.Hx[k+1];
			if (cabsf(T.Hx[k+1])>1E10) for (int kk=P->kMIN; kk<=k+1; kk++) T.Ey[kk]/=1E10, T.Hx[kk]/=1E10;
		}
		T.A = sqrtf(sqrtf(P->ex[P->kMax]));
		T.A /= cabsf(T.Ey[P->kMax] + 0.5*(T.Hx[P->kMax]+T.Hx[P->kMax+1]) / (sqrtf(P->ex[P->kMax]) * cosf(theta*PI/180)));
	}

	if (90 < theta && theta < 270) {
		sinth = sinf(theta*PI/180) * sqrtf(crealf(P->ex[P->kMIN])), sin2th = sq(sinth);
		T.Ey[P->kMax+1] = cexpf(I * casinf(sinf(PI*P->f*P->W->dt) * csqrtf(P->ex[P->kMax+1]-sin2th) * P->dE[P->kMax+1]));
		T.Hx[P->kMax+1] = -csqrtf(P->ex[P->kMax+1] - sin2th);

		for (int k=P->kMax; k>=P->kMIN; k--) {
			T.Ey[k] = T.Ey[k+1] + 2 * I * sinf(PI*P->f*P->W->dt) * P->dH[k+1] * T.Hx[k+1];
			T.Hx[k] = T.Hx[k+1] + 2 * I * sinf(PI*P->f*P->W->dt) * P->dE[k] * T.Ey[k] * (P->ex[k] - sin2th);
			if (cabsf(T.Hx[k+1])>1E10) for (int kk=P->kMax; kk>=k; kk--) T.Ey[kk]/=1E10, T.Hx[kk]/=1E10;
		}
		T.A = sqrtf(sqrtf(P->ex[P->kMIN]));
		T.A /= cabsf(T.Ey[P->kMIN] + 0.5*(T.Hx[P->kMIN]+T.Hx[P->kMIN+1]) / (sqrtf(P->ex[P->kMIN]) * cosf(theta*PI/180)));
	}

	if (theta) {
		MALLOC(T.Hz, P->kMIN, 2-P->kMIN+P->kMax);
		#pragma omp parallel for
//		#pragma GCC ivdep
		for (int k=P->kMIN; k<=P->kMax+1; k++) T.Hz[k] = -T.Ey[k] * sinth;
	}

	if (phi) {
		MALLOC(T.Ex, P->kMIN, 2-P->kMIN+P->kMax);
		MALLOC(T.Hy, P->kMIN, 2-P->kMIN+P->kMax);
		#pragma omp parallel for
//		#pragma GCC ivdep
		for (int k=P->kMIN; k<=P->kMax+1; k++) {
			T.Ex[k] =-T.Ey[k] * sinphi;
			T.Ey[k] = T.Ey[k] * cosphi;
			T.Hy[k] = T.Hx[k] * sinphi;
			T.Hx[k] = T.Hx[k] * cosphi;
		}
	}

	return T;
}


static float complex farFieldX(phaser P, cfield S, int i, tfield T, float kx, float ky)
{
	float complex *Y, *Z;
	float complex EyHz=0, EzHy=0, HyEz=0, HzEy=0, Sum=0;

	MALLOC(Y, P->jMin, 1-P->jMin+P->jMax);
	MALLOC(Z, P->kMin, 1-P->kMin+P->kMax);

	if (T.Hz || T.Ey) {
		FOR1D(j,jMin+1,jMax) Y[j] = P->W->CH->dy[j] * cexpf(I*ky*P->W->jtoy(P->W,j-0.5));
		if (P->W->yMinSur==SYM) FOR1D(j,jMin+1,jMax) Y[j] += conjf(Y[j]) * P->W->cosky;
	}

	if (T.Hz) {
		FOR1D(k,kMin,kMax) Z[k] = P->W->CE->dz[k] * T.Hz[k];
		if (P->W->zMinSur==SYM) FOR1D(k,kMin,kMax) Z[k] -= P->W->CE->dz[k] * T.Hz[-k] * P->W->coskz;
		#pragma omp parallel for reduction (+:EyHz)
		for (int j=P->jMin+1; j<=P->jMax; j++) {
			float complex sum = 0.5 * (S.Ey[i][j][P->kMin] * Z[P->kMin] + S.Ey[i][j][P->kMax] * Z[P->kMax]);
			for (int k=P->kMin+1; k<P->kMax; k++) sum += S.Ey[i][j][k] * Z[k];
			EyHz += Y[j] * sum;
		}
	}

	if (T.Ey) {
		FOR1D(k,kMin,kMax) Z[k] = P->W->CE->dz[k] * T.Ey[k];
		if (P->W->zMinSur==SYM) FOR1D(k,kMin,kMax) Z[k] -= P->W->CE->dz[k] * T.Ey[-k] * P->W->coskz;
		#pragma omp parallel for reduction (+:HzEy)
		for (int j=P->jMin+1; j<=P->jMax; j++) {
			float complex sum = 0.5 * (S.Hz[i][j][P->kMin] * Z[P->kMin] + S.Hz[i][j][P->kMax] * Z[P->kMax]);
			for (int k=P->kMin+1; k<P->kMax; k++) sum += S.Hz[i][j][k] * Z[k];
			HzEy += Y[j] * sum;
		}
	}

	if (T.Hy || T.Ez) {
		FOR1D(j,jMin,jMax) Y[j] = P->W->CE->dy[j] * cexpf(I*ky*P->W->jtoy(P->W,j));
		if (P->W->yMinSur==SYM) FOR1D(j,jMin,jMax) Y[j] -= conjf(Y[j]) * P->W->cosky;
	}

	if (T.Hy) {
		FOR1D(k,kMin+1,kMax) Z[k] = P->W->CH->dz[k] * T.Hy[k];
		if (P->W->zMinSur==SYM) FOR1D(k,kMin+1,kMax) Z[k] += P->W->CH->dz[k] * T.Hy[1-k] * P->W->coskz;
		#pragma omp parallel for reduction (+:EzHy)
		for (int k=P->kMin+1; k<=P->kMax; k++) {
			float complex sum = 0.5 * (S.Ez[i][P->jMin][k] * Y[P->jMin] + S.Ez[i][P->jMax][k] * Y[P->jMax]);
			for (int j=P->jMin+1; j<P->jMax; j++) sum += S.Ez[i][j][k] * Y[j];
			EzHy += Z[k] * sum;
		}
	}

	if (T.Ez) {
		FOR1D(k,kMin+1,kMax) Z[k] = P->W->CH->dz[k] * T.Ez[k];
		if (P->W->zMinSur==SYM) FOR1D(k,kMin+1,kMax) Z[k] += P->W->CH->dz[k] * T.Ez[1-k] * P->W->coskz;
		#pragma omp parallel for reduction (+:HyEz)
		for (int k=P->kMin+1; k<=P->kMax; k++) {
			float complex sum = 0.5 * (S.Hy[i][P->jMin][k] * Y[P->jMin] + S.Hy[i][P->jMax][k] * Y[P->jMax]);
			for (int j=P->jMin+1; j<P->jMax; j++) sum += S.Hy[i][j][k] * Y[j];
			HyEz += Z[k] * sum;
		}
	}

	Sum += (EyHz-HzEy-EzHy+HyEz) * cexpf(I*kx*P->W->itox(P->W,i));
	if (P->W->xMinSur==SYM) Sum += (EyHz+HzEy-EzHy-HyEz) * cexpf(-I*kx*P->W->itox(P->W,i)) * P->W->coskx;

	free(&Y[P->jMin]);
	free(&Z[P->kMin]);
	return Sum;
}


static float complex farFieldY(phaser P, cfield S, int j, tfield T, float kx, float ky)
{
	float complex *X, *Z;
	float complex ExHz=0, EzHx=0, HxEz=0, HzEx=0, Sum=0;

	MALLOC(X, P->iMin, 1-P->iMin+P->iMax);
	MALLOC(Z, P->kMin, 1-P->kMin+P->kMax);

	if (T.Hx || T.Ez) {
		FOR1D(i,iMin,iMax) X[i] = P->W->CE->dx[i] * cexpf(I*kx*P->W->itox(P->W,i));
		if (P->W->xMinSur==SYM) FOR1D(i,iMin,iMax) X[i] -= conjf(X[i]) * P->W->coskx;
	}

	if (T.Hx) {
		FOR1D(k,kMin+1,kMax) Z[k] = P->W->CH->dz[k] * T.Hx[k];
		if (P->W->zMinSur==SYM) FOR1D(k,kMin+1,kMax) Z[k] += P->W->CH->dz[k] * T.Hx[1-k] * P->W->coskz;
		#pragma omp parallel for reduction (+:EzHx)
		for (int k=P->kMin+1; k<=P->kMax; k++) {
			float complex sum = 0.5 * (S.Ez[P->iMin][j][k] * X[P->iMin] + S.Ez[P->iMax][j][k] * X[P->iMax]);
			for (int i=P->iMin+1; i<P->iMax; i++) sum += S.Ez[i][j][k] * X[i];
			EzHx += Z[k] * sum;
		}
	}

	if (T.Ez) {
		FOR1D(k,kMin+1,kMax) Z[k] = P->W->CH->dz[k] * T.Ez[k];
		if (P->W->zMinSur==SYM) FOR1D(k,kMin+1,kMax) Z[k] += P->W->CH->dz[k] * T.Ez[1-k] * P->W->coskz;
		#pragma omp parallel for reduction (+:HxEz)
		for (int k=P->kMin+1; k<=P->kMax; k++) {
			float complex sum = 0.5 * (S.Hx[P->iMin][j][k] * X[P->iMin] + S.Hx[P->iMax][j][k] * X[P->iMax]);
			for (int i=P->iMin+1; i<P->iMax; i++) sum += S.Hx[i][j][k] * X[i];
			HxEz += Z[k] * sum;
		}
	}

	if (T.Hz || T.Ex) {
		FOR1D(i,iMin+1,iMax) X[i] = P->W->CH->dx[i] * cexpf(I*kx*P->W->itox(P->W,i-0.5));
		if (P->W->xMinSur==SYM) FOR1D(i,iMin+1,iMax) X[i] += conjf(X[i]) * P->W->coskx;
	}

	if (T.Hz) {
		FOR1D(k,kMin,kMax) Z[k] = P->W->CE->dz[k] * T.Hz[k];
		if (P->W->zMinSur==SYM) FOR1D(k,kMin,kMax) Z[k] -= P->W->CE->dz[k] * T.Hz[-k] * P->W->coskz;
		#pragma omp parallel for reduction (+:ExHz)
		for (int i=P->iMin+1; i<=P->iMax; i++) {
			float complex sum = 0.5 * (S.Ex[i][j][P->kMin] * Z[P->kMin] + S.Ex[i][j][P->kMax] * Z[P->kMax]);
			for (int k=P->kMin+1; k<P->kMax; k++) sum += S.Ex[i][j][k] * Z[k];
			ExHz += X[i] * sum;
		}
	}

	if (T.Ex) {
		FOR1D(k,kMin,kMax) Z[k] = P->W->CE->dz[k] * T.Ex[k];
		if (P->W->zMinSur==SYM) FOR1D(k,kMin,kMax) Z[k] -= P->W->CE->dz[k] * T.Ex[-k] * P->W->coskz;
		#pragma omp parallel for reduction (+:HzEx)
		for (int i=P->iMin+1; i<=P->iMax; i++) {
			float complex sum = 0.5 * (S.Hz[i][j][P->kMin] * Z[P->kMin] + S.Hz[i][j][P->kMax] * Z[P->kMax]);
			for (int k=P->kMin+1; k<P->kMax; k++) sum += S.Hz[i][j][k] * Z[k];
			HzEx += X[i] * sum;
		}
	}

	Sum += (EzHx-HxEz-ExHz+HzEx) * cexpf(I*ky*P->W->jtoy(P->W,j));
	if (P->W->yMinSur==SYM) Sum += (EzHx+HxEz-ExHz-HzEx) * cexpf(-I*ky*P->W->jtoy(P->W,j)) * P->W->cosky;

	free(&X[P->iMin]);
	free(&Z[P->kMin]);
	return Sum;
}


static float complex farFieldZ(phaser P, cfield S, int k, tfield T, float kx, float ky)
{
	float complex *X, *Y;
	float complex ExHy=0, EyHx=0, HxEy=0, HyEx=0, Sum=0;

	MALLOC(X, P->iMin, 1-P->iMin+P->iMax);
	MALLOC(Y, P->jMin, 1-P->jMin+P->jMax);

	if (T.Hy || T.Ex) {
		FOR1D(i,iMin+1,iMax) X[i] = P->W->CH->dx[i] * cexpf(I*kx*P->W->itox(P->W,i-0.5));
		if (P->W->xMinSur==SYM) FOR1D(i,iMin+1,iMax) X[i] += conjf(X[i]) * P->W->coskx;
		FOR1D(j,jMin,jMax) Y[j] = P->W->CE->dy[j] * cexpf(I*ky*P->W->jtoy(P->W,j));
		if (P->W->yMinSur==SYM) FOR1D(j,jMin,jMax) Y[j] -= conjf(Y[j]) * P->W->cosky;
	}

	if (T.Hy) {
		#pragma omp parallel for reduction (+:ExHy)
		for (int i=P->iMin+1; i<=P->iMax; i++) {
			float complex sum = 0.5 * (S.Ex[i][P->jMin][k] * Y[P->jMin] + S.Ex[i][P->jMax][k] * Y[P->jMax]);
			for (int j=P->jMin+1; j<P->jMax; j++) sum += S.Ex[i][j][k] * Y[j];
			ExHy += X[i] * sum;
		}
		Sum += ExHy * 0.5 * (T.Hy[k]+T.Hy[k+1]);
		if (P->W->zMinSur==SYM) Sum += ExHy * 0.5 * (T.Hy[1-k]+T.Hy[-k]) * P->W->coskz;
	}

	if (T.Ex) {
		#pragma omp parallel for reduction (+:HyEx)
		for (int i=P->iMin+1; i<=P->iMax; i++) {
			float complex sum = 0.5 * (S.Hy[i][P->jMin][k] * Y[P->jMin] + S.Hy[i][P->jMax][k] * Y[P->jMax]);
			for (int j=P->jMin+1; j<P->jMax; j++) sum += S.Hy[i][j][k] * Y[j];
			HyEx += X[i] * sum;
		}
		Sum -= HyEx * T.Ex[k];
		if (P->W->zMinSur==SYM) Sum += HyEx * T.Ex[-k] * P->W->coskz;
	}

	if (T.Hx || T.Ey) {
		FOR1D(i,iMin,iMax) X[i] = P->W->CE->dx[i] * cexpf(I*kx*P->W->itox(P->W,i));
		if (P->W->xMinSur==SYM) FOR1D(i,iMin,iMax) X[i] -= conjf(X[i]) * P->W->coskx;
		FOR1D(j,jMin+1,jMax) Y[j] = P->W->CH->dy[j] * cexpf(I*ky*P->W->jtoy(P->W,j-0.5));
		if (P->W->yMinSur==SYM) FOR1D(j,jMin+1,jMax) Y[j] += conjf(Y[j]) * P->W->cosky;
	}

	if (T.Hx) {
		#pragma omp parallel for reduction (+:EyHx)
		for (int j=P->jMin+1; j<=P->jMax; j++) {
			float complex sum = 0.5 * (S.Ey[P->iMin][j][k] * X[P->iMin] + S.Ey[P->iMax][j][k] * X[P->iMax]);
			for (int i=P->iMin+1; i<P->iMax; i++) sum += S.Ey[i][j][k] * X[i];
			EyHx += Y[j] * sum;
		}
		Sum -= EyHx * 0.5 * (T.Hx[k]+T.Hx[k+1]);
		if (P->W->zMinSur==SYM) Sum -= EyHx * 0.5 * (T.Hx[1-k]+T.Hx[-k]) * P->W->coskz;
	}

	if (T.Ey) {
		#pragma omp parallel for reduction (+:HxEy)
		for (int j=P->jMin+1; j<=P->jMax; j++) {
			float complex sum = 0.5 * (S.Hx[P->iMin][j][k] * X[P->iMin] + S.Hx[P->iMax][j][k] * X[P->iMax]);
			for (int i=P->iMin+1; i<P->iMax; i++) sum += S.Hx[i][j][k] * X[i];
			HxEy += Y[j] * sum;
		}
		Sum += HxEy * T.Ey[k];
		if (P->W->zMinSur==SYM) Sum -= HxEy * T.Ey[-k] * P->W->coskz;
	}

	free(&X[P->iMin]);
	free(&Y[P->jMin]);
	return Sum;
}


float complex farField(phaser P, tfield pol(phaser,float,float), float theta, float phi)
{
	float sinth, kx, ky;
	float complex Sum = 0;

	if (!P->dE) setFarfieldCoeffs(P, P->W);
	if (theta == INF || phi == INF) return 0;
	while (theta < 0) theta += 360;
	while (theta >= 360) theta -= 360;
	while (phi < 0) phi += 360;
	while (phi >= 360) phi -= 360;
	if ( 89 < theta && theta <  91) return 0.5 * (farField(P, pol,  89, phi) + farField(P, pol,  91, phi));
	if (269 < theta && theta < 271) return 0.5 * (farField(P, pol, 269, phi) + farField(P, pol, 271, phi));
	if (theta < 90 || 270 < theta) {
		if (cimagf(P->ex[P->kMax])) return 0;
		else sinth = sinf(theta*PI/180) * sqrtf(crealf(P->ex[P->kMax]));
	}
	if (90 < theta && theta < 270) {
		if (cimagf(P->ex[P->kMIN])) return 0;
		else sinth = sinf(theta*PI/180) * sqrtf(crealf(P->ex[P->kMIN]));
	}

	tfield T = pol(P, theta, phi);
	kx = - 2 * PI * P->f * sinth * cosf(phi*PI/180);
	ky = - 2 * PI * P->f * sinth * sinf(phi*PI/180);

	if (P->W->xMinSur==PML) Sum -= farFieldX(P, P->xMin, P->iMin, T, kx, ky);
	if (P->W->xMaxSur==PML) Sum += farFieldX(P, P->xMax, P->iMax, T, kx, ky);
	if (P->W->yMinSur==PML) Sum -= farFieldY(P, P->yMin, P->jMin, T, kx, ky);
	if (P->W->yMaxSur==PML) Sum += farFieldY(P, P->yMax, P->jMax, T, kx, ky);
	if (P->W->zMinSur==PML) Sum -= farFieldZ(P, P->zMin, P->kMin, T, kx, ky);
	if (P->W->zMaxSur==PML) Sum += farFieldZ(P, P->zMax, P->kMax, T, kx, ky);

	if (T.Ex) free(&T.Ex[P->kMIN]);
	if (T.Ey) free(&T.Ey[P->kMIN]);
	if (T.Ez) free(&T.Ez[P->kMIN]);
	if (T.Hx) free(&T.Hx[P->kMIN]);
	if (T.Hy) free(&T.Hy[P->kMIN]);
	if (T.Hz) free(&T.Hz[P->kMIN]);

	if (P->xNUM) {
		float complex s = 0, q = 1, p = P->expikx * cexpf(I*kx*(P->W->xMax-P->W->xMin));
		for (int x=0; x<P->xNUM; x++) s += (q *= p);
		Sum *= s / P->xNUM;
	}
	if (P->yNUM) {
		float complex s = 0, q = 1, p = P->expiky * cexpf(I*ky*(P->W->yMax-P->W->yMin));
		for (int y=0; y<P->yNUM; y++) s += (q *= p);
		Sum *= s / P->yNUM;
	}

	return Sum * T.A * P->A * P->f;
}


float farFieldI(phaser P, tfield pol(phaser,float,float), float theta, float phi)
{
	if (pol == Unpol) return farFieldI(P, Ppol, theta, phi) + farFieldI(P, Spol, theta, phi);
	float complex f = farField(P, pol, theta, phi);
	return conjf(f) * f;
}


void farFieldProfile(world W, phaser P, projmap projMap, int X, int Y, ...)
{
	char *format, fieldname[4], filename[1024], fullname[1024];
	va_list ap;
	va_start(ap, Y);
	method M = va_arg(ap, method);
	unsigned long int *cmap = (M == png) ? va_arg(ap, unsigned long int *) : 0;
	float Max = (M == png) ? va_arg(ap, double) : 1;
	float max=0, angle[2];
	format = va_arg(ap, char*);
	vsprintf(filename, format, ap);
	strcat(strcpy(fullname, W->ID), filename);
	strcpy(fieldname, strstr(filename, "Log") ? "Log" : "");
	if (filename[0]=='/') mkdir(W->ID, 0755);
	if (fullname[strlen(fullname)-1]=='/') mkdir(fullname, 0755);

	slice S;
	MALLOC(S, 0, 1);
	*S = (struct slice) {0, 0, 0, X-1, 0, Y-1, 0, 0, 0, 0, X-1, X, 0, Y-1, Y};
	float ***Up = makeField(S->iMIN, S->iMAX, S->jMIN, S->jMAX, S->kMIN, S->kMAX);
	float ***Pp = makeField(S->iMIN, S->iMAX, S->jMIN, S->jMAX, S->kMIN, S->kMAX);
	float ***Sp = makeField(S->iMIN, S->iMAX, S->jMIN, S->jMAX, S->kMIN, S->kMAX);
	float ***Lp = makeField(S->iMIN, S->iMAX, S->jMIN, S->jMAX, S->kMIN, S->kMAX);
	float ***Rp = makeField(S->iMIN, S->iMAX, S->jMIN, S->jMAX, S->kMIN, S->kMAX);
	float ***Xp = makeField(S->iMIN, S->iMAX, S->jMIN, S->jMAX, S->kMIN, S->kMAX);
	float ***Yp = makeField(S->iMIN, S->iMAX, S->jMIN, S->jMAX, S->kMIN, S->kMAX);

	int n=1, N=0;
	for (int x=0; x<X; x++) for (int y=0; y<Y; y++) {
		projMap(x, y, X, Y, angle);
		if (angle[0] != INF) N++;
	}
	printf("Performing Near-to-Far-Field Transformation.\n");
	timer(0, N);
	for (int x=0; x<X; x++) for (int y=0; y<Y; y++) {
		projMap(x, y, X, Y, angle);
		float complex p = farField(P, Ppol, angle[0], angle[1]);
		float complex s = farField(P, Spol, angle[0], angle[1]);
		if (fieldname[0]) {
			Xp[0][x][y] = log10f(conjf(p*cosf(angle[1]*PI/180)-s*sinf(angle[1]*PI/180))*(p*cosf(angle[1]*PI/180)-s*sinf(angle[1]*PI/180)));
			Yp[0][x][y] = log10f(conjf(s*cosf(angle[1]*PI/180)+p*sinf(angle[1]*PI/180))*(s*cosf(angle[1]*PI/180)+p*sinf(angle[1]*PI/180)));
			Lp[0][x][y] = log10f(sq(crealf(p)+cimagf(s))+sq(cimagf(p)-crealf(s)));
			Rp[0][x][y] = log10f(sq(crealf(p)-cimagf(s))+sq(cimagf(p)+crealf(s)));
			Pp[0][x][y] = log10f(conjf(p)*p);
			Sp[0][x][y] = log10f(conjf(s)*s);
			Up[0][x][y] = log10f(conjf(p)*p+conjf(s)*s);
		}
		else {
			Xp[0][x][y] = conjf(p*cosf(angle[1]*PI/180)-s*sinf(angle[1]*PI/180))*(p*cosf(angle[1]*PI/180)-s*sinf(angle[1]*PI/180));
			Yp[0][x][y] = conjf(s*cosf(angle[1]*PI/180)+p*sinf(angle[1]*PI/180))*(s*cosf(angle[1]*PI/180)+p*sinf(angle[1]*PI/180));
			Lp[0][x][y] = sq(crealf(p)+cimagf(s))+sq(cimagf(p)-crealf(s));
			Rp[0][x][y] = sq(crealf(p)-cimagf(s))+sq(cimagf(p)+crealf(s));
			Pp[0][x][y] = conjf(p)*p;
			Sp[0][x][y] = conjf(s)*s;
			Up[0][x][y] = Pp[0][x][y]+Sp[0][x][y];
		}
		max = MAX(max, Up[0][x][y]);
		if (angle[0] != INF) timer(++n, N);
	}

	if (Max <=0) Max = max;
	S->F = Xp; M(W, S, Max, fieldname, strcat(strcpy(filename, fullname), "X-pol"), cmap); removeField(Xp, 0, 0, 0);
	S->F = Yp; M(W, S, Max, fieldname, strcat(strcpy(filename, fullname), "Y-pol"), cmap); removeField(Yp, 0, 0, 0);
	S->F = Lp; M(W, S, Max, fieldname, strcat(strcpy(filename, fullname), "CL-pol"), cmap); removeField(Lp, 0, 0, 0);
	S->F = Rp; M(W, S, Max, fieldname, strcat(strcpy(filename, fullname), "CR-pol"), cmap); removeField(Rp, 0, 0, 0);
	S->F = Pp; M(W, S, Max, fieldname, strcat(strcpy(filename, fullname), "P-pol"), cmap); removeField(Pp, 0, 0, 0);
	S->F = Sp; M(W, S, Max, fieldname, strcat(strcpy(filename, fullname), "S-pol"), cmap); removeField(Sp, 0, 0, 0);
	S->F = Up; M(W, S, Max, fieldname, strcat(strcpy(filename, fullname), "Unpol"), cmap);
	deleteSlice(S);
}


void farFieldTheta(world W, phaser P, float phi, char *format, ...)
{
	char filename[1024];
	va_list ap;
	va_start(ap, format);
	vsprintf(filename, format, ap);

	for (float theta=0; theta<360; theta+=1) {
		float complex p = farField(P, Ppol, theta, phi);
		float complex s = farField(P, Spol, theta, phi);
		writeRow(W, filename, theta, sq(cabsf(p))+sq(cabsf(s)), sq(cabsf(p)), sq(cabsf(s)));
	}
}


void farFieldPhi(world W, phaser P, float theta, char *format, ...)
{
	char filename[1024];
	va_list ap;
	va_start(ap, format);
	vsprintf(filename, format, ap);

	for (float phi=0; phi<360; phi+=1) {
		float complex p = farField(P, Ppol, theta, phi);
		float complex s = farField(P, Spol, theta, phi);
		writeRow(W, filename, phi, sq(cabsf(p))+sq(cabsf(s)), sq(cabsf(p)), sq(cabsf(s)));
	}
}


float (farFieldFlux)(phaser P, tfield pol(phaser,float,float), float theta, float phi, float aMin, float aMax, float delta, ...)
{
	if (!delta) delta=1;
	if (!aMax) aMax=aMin, aMin=0;

	float Sum = 0;
	int nMin = floor(aMin/delta+0.5);
   	int nMax = ceil(aMax/delta-0.5);

	if (!nMin) {
		Sum = 2 * PI * (cosf(aMin*PI/180) - cosf(MIN(aMax,0.5*delta)*PI/180)) * farFieldI(P,pol,theta,phi);
		aMin = 0.5 * delta;
	   	nMin++;
	}
	if (nMax) {
		for (int n=nMin; n<=nMax; n++) {
			int N = n<floor(60/delta+0.5) ? n*6 : floor(60/delta+0.5)*6;
			float dth = cosf(MAX(aMin,delta*(n-0.5))*PI/180) - cosf(MIN(aMax,delta*(n+0.5))*PI/180);

			for (int m=0; m<N; m++) {
				float a = rotateA(theta, delta*n, 360.0*m/N);
				float b = rotateB(theta, delta*n, 360.0*m/N)+phi;
				float c = rotateB(-a, theta, phi-b)+180-360.0*m/N;
				float complex p = farField(P, Ppol, a, b);
				float complex s = farField(P, Spol, a, b);
				float pp = sq(cabsf((p * cosf(c*PI/180) - s * sinf(c*PI/180))));
				float ss = sq(cabsf((s * cosf(c*PI/180) + p * sinf(c*PI/180))));

				if (pol == Ppol) Sum += dth * (pp * (PI/N + sinf(PI/N)) + ss * (PI/N - sinf(PI/N)));
				if (pol == Spol) Sum += dth * (pp * (PI/N - sinf(PI/N)) + ss * (PI/N + sinf(PI/N)));
				if (pol == Unpol) Sum += 2 * PI/N * dth * (pp + ss);
			}
		}
	}
	return Sum;
}
