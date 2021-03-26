/**********************************************************************
 * ALiS is Copyright (C) 2009-2019 by Ho-Seok Ee <hsee@kongju.ac.kr>. *
 * Redistribution and use with or without modification, are permitted *
 * under the terms of the Artistic License version 2.                 *
 **********************************************************************/

#include <sys/stat.h>
#include <string.h>
#include "alis_c.h"


#define makeComplexFields(P,S,E1,E2,H1,H2,iMin,iMax,jMin,jMax,kMin,kMax) \
{ \
	P.S.E1 = makeComplexField(P.iMin, P.iMax, P.jMin, P.jMax, P.kMin, P.kMax); \
	P.S.E2 = makeComplexField(P.iMin, P.iMax, P.jMin, P.jMax, P.kMin, P.kMax); \
	P.S.H1 = makeComplexField(P.iMin, P.iMax, P.jMin, P.jMax, P.kMin, P.kMax); \
	P.S.H2 = makeComplexField(P.iMin, P.iMax, P.jMin, P.jMax, P.kMin, P.kMax); \
}


phaser (createPhasers)(world W, float LMin, float LMax, float dL, float xMin, float xMax, float yMin, float yMax, float zMin, float zMax, ...)
{
	int nMin = floor(LMin/dL), nMax = ceil(LMax/dL);
	phaser P;
	source S = W->SE && W->N==W->SE->N ? W->SE : W->SH && W->N==W->SH->N ? W->SH : 0;
	for (world w = W; !S && w->S && w->S->W; S = w->SE && W->N==w->SE->N ? w->SE : w->SH && W->N==w->SH->N ? w->SH : 0) w = W->S->W;
	MALLOC(P, 0, 2-nMin+nMax);

	for (int n=0; n<=nMax-nMin; n++) {
		P[n].W = W;
		P[n].N = nMin * nMax * dL / W->dt;
		P[n].f = (nMax-n) / (dL * nMin * nMax);
		P[n].A = 2.0 / (P->W->N ? P->W->N : P->N);
		P[n].iMin = MAX(W->iMin, W->xtoi(W, xMin)), P[n].iMax = MIN(W->iMax, W->xtoi(W, xMax));
		P[n].jMin = MAX(W->jMin, W->ytoj(W, yMin)), P[n].jMax = MIN(W->jMax, W->ytoj(W, yMax));
		P[n].kMin = MAX(W->kMin, W->ztok(W, zMin)), P[n].kMax = MIN(W->kMax, W->ztok(W, zMax));

		if (S && W->N) {
			float complex c = 0;
			float complex iwdt = 2 * PI * I * P[n].f * W->dt;
			for (int n=1; n<W->N; n++) c += S->waveForm(W, S, W->dt*n, S->phase) * cexpf(iwdt*n);
			P[n].A = 1 / cabsf(c);
			if (W->sourceType == 2) P[n].A /= sqrtf(P[n].f * S->wavelength);
			if (W->sourceType == 3) P[n].A /= P[n].f * S->wavelength;
		}

		if ((xMax-xMin < INF && xMin < W->xMin && W->xMax < xMax) && (W->xMaxSur==SYM || W->xMaxSur==PBC || W->xMaxSur==BBC)) {
			P[n].xNUM = floor((xMax - xMin) / (W->xMax - W->xMin) + 0.5);
			P[n].expikx = (W->xMaxSur==BBC) ? (W->coskx+I*W->sinkx) : 1;
		}

		if ((yMax-yMin < INF && yMin < W->yMin && W->yMax < yMax) && (W->yMaxSur==SYM || W->yMaxSur==PBC || W->yMaxSur==BBC)) {
			P[n].yNUM = floor((yMax - yMin) / (W->yMax - W->yMin) + 0.5);
			P[n].expiky = (W->yMaxSur==BBC) ? (W->cosky+I*W->sinky) : 1;
		}

		if (W->xMinSur==PML && P[n].iMin != P[n].iMax) makeComplexFields(P[n], xMin, Ey, Ez, Hy, Hz, iMin, iMin, jMin, jMax, kMin, kMax);
		if (W->yMinSur==PML && P[n].jMin != P[n].jMax) makeComplexFields(P[n], yMin, Ez, Ex, Hz, Hx, iMin, iMax, jMin, jMin, kMin, kMax);
		if (W->zMinSur==PML && P[n].kMin != P[n].kMax) makeComplexFields(P[n], zMin, Ex, Ey, Hx, Hy, iMin, iMax, jMin, jMax, kMin, kMin);
		if (W->xMaxSur==PML) makeComplexFields(P[n], xMax, Ey, Ez, Hy, Hz, iMax, iMax, jMin, jMax, kMin, kMax);
		if (W->yMaxSur==PML) makeComplexFields(P[n], yMax, Ez, Ex, Hz, Hx, iMin, iMax, jMax, jMax, kMin, kMax);
		if (W->zMaxSur==PML) makeComplexFields(P[n], zMax, Ex, Ey, Hx, Hy, iMin, iMax, jMin, jMax, kMax, kMax);
	}

	return P;
}


phaser (createPhaser)(world W, float L, float xMin, float xMax, float yMin, float yMax, float zMin, float zMax, ...)
{
	return (createPhasers)(W, L, L, L, xMin, xMax, yMin, yMax, zMin, zMax);
}


phaser choosePhaser(phaser P, float wavelength)
{
	do P++;
	while (1/wavelength < P->f);
	return (1/wavelength - P->f < (P-1)->f - 1/wavelength) ? P : P-1;
}


#define update(W,P,S,E1,E2,H1,H2,F1,F2,expiwtE,expiwtH,i,iMin,iMax,j,jMin,jMax,ii,jj,kk,iii,jjj,kkk) \
if (P->S.E1 && P->S.E2 && P->S.H1 && P->S.H2) \
_Pragma("omp parallel for") for (int i=iMin; i<=iMax; i++) \
_Pragma("GCC ivdep") for (int j=jMin; j<=jMax; j++) { \
	P->S.E1[ii][jj][kk] += expiwtE * W->E->F1[ii][jj][kk]; \
	P->S.E2[ii][jj][kk] += expiwtE * W->E->F2[ii][jj][kk]; \
	P->S.H1[ii][jj][kk] += expiwtH * (W->H->F1[ii][jj][kk] + W->H->F1[iii][jjj][kkk]) * 0.5; \
	P->S.H2[ii][jj][kk] += expiwtH * (W->H->F2[ii][jj][kk] + W->H->F2[iii][jjj][kkk]) * 0.5; \
}


void updatePhaser(world W, phaser P)
{
	do {
		float complex expiwtE = cexpf(2 * PI * I * P->f * W->E->t);
		float complex expiwtH = cexpf(2 * PI * I * P->f * W->H->t);

		update(W, P, xMin, Ey,Ez, Hy,Hz, y,z, expiwtE, expiwtH, j,P->jMin,P->jMax, k,P->kMin,P->kMax, P->iMin,j,k, P->iMin+1,j,k);
		update(W, P, xMax, Ey,Ez, Hy,Hz, y,z, expiwtE, expiwtH, j,P->jMin,P->jMax, k,P->kMin,P->kMax, P->iMax,j,k, P->iMax+1,j,k);
		update(W, P, yMin, Ez,Ex, Hz,Hx, z,x, expiwtE, expiwtH, i,P->iMin,P->iMax, k,P->kMin,P->kMax, i,P->jMin,k, i,P->jMin+1,k);
		update(W, P, yMax, Ez,Ex, Hz,Hx, z,x, expiwtE, expiwtH, i,P->iMin,P->iMax, k,P->kMin,P->kMax, i,P->jMax,k, i,P->jMax+1,k);
		update(W, P, zMin, Ex,Ey, Hx,Hy, x,y, expiwtE, expiwtH, i,P->iMin,P->iMax, j,P->jMin,P->jMax, i,j,P->kMin, i,j,P->kMin+1);
		update(W, P, zMax, Ex,Ey, Hx,Hy, x,y, expiwtE, expiwtH, i,P->iMin,P->iMax, j,P->jMin,P->jMax, i,j,P->kMax, i,j,P->kMax+1);

		P++;
	}
	while (P->f);
}


#define removeFields(P,S,E1,E2,H1,H2,iMin,jMin,kMin) \
{ \
	if (P->S.E1) removeField(P->S.E1, P->iMin, P->jMin, P->kMin); \
	if (P->S.E2) removeField(P->S.E2, P->iMin, P->jMin, P->kMin); \
	if (P->S.H1) removeField(P->S.H1, P->iMin, P->jMin, P->kMin); \
	if (P->S.H2) removeField(P->S.H2, P->iMin, P->jMin, P->kMin); \
}


void deletePhaser(phaser P)
{
	do {
		removeFields(P, xMin, Ey, Ez, Hy, Hz, iMin, jMin, kMin);
		removeFields(P, xMax, Ey, Ez, Hy, Hz, iMax, jMin, kMin);
		removeFields(P, yMin, Ez, Ex, Hz, Hx, iMin, jMin, kMin);
		removeFields(P, yMax, Ez, Ex, Hz, Hx, iMin, jMax, kMin);
		removeFields(P, zMin, Ex, Ey, Hx, Hy, iMin, jMin, kMin);
		removeFields(P, zMax, Ex, Ey, Hx, Hy, iMin, jMin, kMax);
		if (&P->dE[P->kMIN]) free(&P->dE[P->kMIN]);
		if (&P->dH[P->kMIN]) free(&P->dH[P->kMIN]);
		if (&P->ex[P->kMIN]) free(&P->ex[P->kMIN]);
		if (&P->ez[P->kMIN]) free(&P->ez[P->kMIN]);
		P++;
	}
	while (P->f);
	free(P);
}


static float phaserPoyntingX(phaser P, cfield S, int i)
{
	float Sum=0;

	if (S.Ey && S.Hz) {
		#pragma omp parallel for reduction(+:Sum)
		for (int j=P->jMin+1; j<=P->jMax; j++) {
			float sum=0;
			sum += 0.5 * P->W->CE->dz[P->kMin] * S.Ey[i][j][P->kMin] * conjf(S.Hz[i][j][P->kMin]);
			sum += 0.5 * P->W->CE->dz[P->kMax] * S.Ey[i][j][P->kMax] * conjf(S.Hz[i][j][P->kMax]);
			for (int k=P->kMin+1; k<P->kMax; k++) sum += P->W->CE->dz[k] * S.Ey[i][j][k] * conjf(S.Hz[i][j][k]);
			Sum += P->W->CH->dy[j] * sum;
		}
	}
	if (S.Ez && S.Hy) {
		#pragma omp parallel for reduction(+:Sum)
		for (int k=P->kMin+1; k<=P->kMax; k++) {
			float sum=0;
			sum -= 0.5 * P->W->CE->dy[P->jMin] * S.Ez[i][P->jMin][k] * conjf(S.Hy[i][P->jMin][k]);
			sum -= 0.5 * P->W->CE->dy[P->jMax] * S.Ez[i][P->jMax][k] * conjf(S.Hy[i][P->jMax][k]);
			for (int j=P->jMin+1; j<P->jMax; j++) sum -= P->W->CE->dy[j] * S.Ez[i][j][k] * conjf(S.Hy[i][j][k]);
			Sum += P->W->CH->dz[k] * sum;
		}
	}
	return Sum;
}


static float phaserPoyntingY(phaser P, cfield S, int j)
{
	float Sum=0;

	if (S.Ez && S.Hx) {
		#pragma omp parallel for reduction(+:Sum)
		for (int k=P->kMin+1; k<=P->kMax; k++) {
			float sum=0;
			sum += 0.5 * P->W->CE->dx[P->iMin] * S.Ez[P->iMin][j][k] * conjf(S.Hx[P->iMin][j][k]);
			sum += 0.5 * P->W->CE->dx[P->iMax] * S.Ez[P->iMax][j][k] * conjf(S.Hx[P->iMax][j][k]);
			for (int i=P->iMin+1; i<P->iMax; i++) sum += P->W->CE->dx[i] * S.Ez[i][j][k] * conjf(S.Hx[i][j][k]);
			Sum += P->W->CH->dz[k] * sum;
		}
	}
	if (S.Ex && S.Hz) {
		#pragma omp parallel for reduction(+:Sum)
		for (int i=P->iMin+1; i<=P->iMax; i++) {
			float sum=0;
			sum -= 0.5 * P->W->CE->dz[P->kMin] * S.Ex[i][j][P->kMin] * conjf(S.Hz[i][j][P->kMin]);
			sum -= 0.5 * P->W->CE->dz[P->kMax] * S.Ex[i][j][P->kMax] * conjf(S.Hz[i][j][P->kMax]);
			for (int k=P->kMin+1; k<P->kMax; k++) sum -= P->W->CE->dz[k] * S.Ex[i][j][k] * conjf(S.Hz[i][j][k]);
			Sum += P->W->CH->dx[i] * sum;
		}
	}
	return Sum;
}


static float phaserPoyntingZ(phaser P, cfield S, int k)
{
	float Sum=0;

	if (S.Ex && S.Hy) {
		#pragma omp parallel for reduction(+:Sum)
		for (int i=P->iMin+1; i<=P->iMax; i++) {
			float sum=0;
			sum += 0.5 * P->W->CE->dy[P->jMin] * S.Ex[i][P->jMin][k] * conjf(S.Hy[i][P->jMin][k]);
			sum += 0.5 * P->W->CE->dy[P->jMax] * S.Ex[i][P->jMax][k] * conjf(S.Hy[i][P->jMax][k]);
			for (int j=P->jMin+1; j<P->jMax; j++) sum += P->W->CE->dy[j] * S.Ex[i][j][k] * conjf(S.Hy[i][j][k]);
			Sum += P->W->CH->dx[i] * sum;
		}
	}
	if (S.Ey && S.Hx) {
		#pragma omp parallel for reduction(+:Sum)
		for (int j=P->jMin+1; j<=P->jMax; j++) {
			float sum=0;
			sum -= 0.5 * P->W->CE->dx[P->iMin] * S.Ey[P->iMin][j][k] * conjf(S.Hx[P->iMin][j][k]);
			sum -= 0.5 * P->W->CE->dx[P->iMax] * S.Ey[P->iMax][j][k] * conjf(S.Hx[P->iMax][j][k]);
			for (int i=P->iMin+1; i<P->iMax; i++) sum -= P->W->CE->dx[i] * S.Ey[i][j][k] * conjf(S.Hx[i][j][k]);
			Sum += P->W->CH->dy[j] * sum;
		}
	}
	return Sum;
}


float phaserPoynting(phaser P) { return fabsf(phaserPoyntingOut(P)); }
float phaserPoyntingIn(phaser P) { return -phaserPoyntingOut(P); }
float phaserPoyntingOut(phaser P)
{
	float Sum = 0;

	Sum += phaserPoyntingX(P, P->xMax, P->iMax) - phaserPoyntingX(P, P->xMin, P->iMin);
	Sum += phaserPoyntingY(P, P->yMax, P->jMax) - phaserPoyntingY(P, P->yMin, P->jMin);
	Sum += phaserPoyntingZ(P, P->zMax, P->kMax) - phaserPoyntingZ(P, P->zMin, P->kMin);

	if (P->W->Dom.x.Min==P->W->Dom.x.Max) Sum /= P->W->CH->dx[P->iMax]+P->W->CH->dx[P->iMax-1];
	if (P->W->Dom.y.Min==P->W->Dom.y.Max) Sum /= P->W->CH->dy[P->jMax]+P->W->CH->dy[P->jMax-1];
	if (P->W->Dom.z.Min==P->W->Dom.z.Max) Sum /= P->W->CH->dz[P->kMax]+P->W->CH->dz[P->kMax-1];

	if (P->W->xMinSur==SYM) Sum *= 2;
	if (P->W->yMinSur==SYM) Sum *= 2;
	if (P->W->zMinSur==SYM) Sum *= 2;

	return Sum * sq(P->A);
}


void (writePoyntingSpectrum)(world W, char *format, ...)
{
	int M=0;
	char str[1024], filename[1024];
	va_list ap;
	phaser *P;

	va_start(ap, format);
	vsprintf(filename, format, ap);
	for (phaser p=va_arg(ap, phaser); p; M++) p=va_arg(ap, phaser);

	va_start(ap, format);
	vsprintf(filename, format, ap);
	MALLOC(P, 0, M);
	for (int m=0; m<M; m++) P[m] = va_arg(ap, phaser);
	for (int n=0; P[0][n].f; n++) {
		int s = sprintf(str, "%g", 1/P[0][n].f);
		for (int m=0; m<M; m++) s += sprintf(str+s, "\t%g", phaserPoynting(P[m]+n));
		sprintf(str+s, "\n");
		writeTxt(W, filename, str);
	}
	free(P);
}


void (writeFarFieldSpectrum)(world W, tfield pol(phaser,float,float), float theta, float phi, char *format, ...)
{
	int M=0;
	char str[1024], filename[1024];
	va_list ap;
	phaser *P;

	va_start(ap, format);
	vsprintf(filename, format, ap);
	for (phaser p=va_arg(ap, phaser); p; M++) p=va_arg(ap, phaser);
	va_start(ap, format);
	vsprintf(filename, format, ap);

	MALLOC(P, 0, M);
	for (int m=0; m<M; m++) P[m] = va_arg(ap, phaser);
	for (int n=0; P[0][n].f; n++) {
		int s = sprintf(str, "%g", 1/P[0][n].f);
		for (int m=0; m<M; m++) s += sprintf(str+s, "\t%g", farFieldI(P[m]+n, pol, theta, phi));
		sprintf(str+s, "\n");
		writeTxt(W, filename, str);
	}
	free(P);
}
