/**********************************************************************
 * ALiS is Copyright (C) 2009-2019 by Ho-Seok Ee <hsee@kongju.ac.kr>. *
 * Redistribution and use with or without modification, are permitted *
 * under the terms of the Artistic License version 2.                 *
 **********************************************************************/

#include <sys/time.h>
#include <sys/stat.h>
#include <string.h>
#include "alis_c.h"


typedef struct queue *queue;
struct queue {
	world W;
	char *ID;
	FILE *File;
	unsigned long int N;
	unsigned long long int timeStamp;
	float complex **Spectrum;
	queue Q;
};


static float jx(world W, int n, int i, int j, int k)
{
	float f=0;
	if (W->CE->ax[n]) f += 2 * W->CE->Jx[n][i][j][k] * W->E->x[i][j][k] / (1 + 1 / W->CE->ax[n]);
	for (int p=0; p<W->CE->NK[n]; p++) if (W->E->Jx[n][p]) f += W->E->Jx[n][p][i][j][k];
	return f / W->CE->x[i][j][k];
}

static float jy(world W, int n, int i, int j, int k)
{
	float f=0;
	if (W->CE->ay[n]) f += 2 * W->CE->Jy[n][i][j][k] * W->E->y[i][j][k] / (1 + 1 / W->CE->ay[n]);
	for (int p=0; p<W->CE->NK[n]; p++) if (W->E->Jy[n][p]) f += W->E->Jy[n][p][i][j][k];
	return f / W->CE->y[i][j][k];
}

static float jz(world W, int n, int i, int j, int k)
{
	float f=0;
	if (W->CE->az[n]) f += 2 * W->CE->Jz[n][i][j][k] * W->E->z[i][j][k] / (1 + 1 / W->CE->az[n]);
	for (int p=0; p<W->CE->NK[n]; p++) if (W->E->Jz[n][p]) f += W->E->Jz[n][p][i][j][k];
	return f / W->CE->z[i][j][k];
}


float get(world W, field F, float x, float y, float z)
{
	return F(W, W->xtoi(W, x), W->ytoj(W, y), W->ztok(W, z));
}


float worldMax(world W, field F)
{
	float Max=0;
	#pragma omp parallel
	{
		float f, max=0;
		#pragma omp for
		for (int i=W->iMin; i<W->iMax; i++) for (int j=W->jMin; j<W->jMax; j++) for (int k=W->kMin; k<W->kMax; k++)
			if (max < (f=fabsf(F(W, i, j, k)))) max = f;
		#pragma omp critical
		if (Max < max) Max = max;
	}
	return Max;
}


static float worldSumZ(world W, field F, int i, int j)
{
	float Sum=0;
	for (int k=W->kMin+1; k<W->kMax; k++) Sum += W->CE->dz[k] * F(W, i, j, k);
	Sum += 0.5 * W->CE->dz[W->kMin] * F(W, i, j, W->kMin);
	Sum += 0.5 * W->CE->dz[W->kMax] * F(W, i, j, W->kMax);
	return Sum;
}

static float worldSumY(world W, field F, int i)
{
	float Sum=0;
	for (int j=W->jMin+1; j<W->jMax; j++) Sum += W->CE->dy[j] * worldSumZ(W, F, i, j);
	Sum += 0.5 * W->CE->dy[W->jMin] * worldSumZ(W, F, i, W->jMin);
	Sum += 0.5 * W->CE->dy[W->jMax] * worldSumZ(W, F, i, W->jMax);
	return Sum;
}

float worldSum(world W, field F)
{
	float Sum=0;
	#pragma omp parallel for reduction(+:Sum)
	for (int i=W->iMin+1; i<W->iMax; i++) Sum += W->CE->dx[i] * worldSumY(W, F, i);
	Sum += 0.5 * W->CE->dx[W->iMin] * worldSumY(W, F, W->iMin);
	Sum += 0.5 * W->CE->dx[W->iMax] * worldSumY(W, F, W->iMax);
	if (W->xMinSur==SYM) Sum *= 2;
	if (W->yMinSur==SYM) Sum *= 2;
	if (W->zMinSur==SYM) Sum *= 2;
	return Sum;
}


static queue searchQueue(world W, char *ID, queue *staticQ)
{
	if (!*staticQ) {
		MALLOC(*staticQ, 0, 1);
		**staticQ = (struct queue) {0};
	}
	queue Q = *staticQ;

	while (Q->Q && (W != Q->W || strcmp(ID, Q->ID))) Q = Q->Q;
	if (Q->Q) {
		strcat(strcpy(ID, W->ID), Q->ID);
	}
	else {
		Q->W = W;
		MALLOC(Q->ID, 0, 1+strlen(ID));
		MALLOC(Q->Q, 0, 1);
		strcpy(Q->ID, ID);
		strcat(strcpy(ID, W->ID), Q->ID);
		*Q->Q = (struct queue) {0};
		if (Q->ID[0]=='/') mkdir(W->ID, 0755);
		if (ID[strlen(ID)-1]=='/') mkdir(ID, 0755);
	}
	Q->N++;
	return Q;
}


void (writeTxt)(world W, char *format, ...)
{
	static queue staticQ;
	char str[1024], filename[1024];
	va_list ap;
	va_start(ap, format);
	vsprintf(filename, format, ap);
	format = va_arg(ap, char*);
	vsprintf(str, format, ap);
	queue Q = searchQueue(W, filename, &staticQ);
	struct timeval now;
	gettimeofday(&now, 0);
	unsigned long long int time = now.tv_sec * 1000000LL + now.tv_usec;

	if (!Q->timeStamp) {
		Q->File = fopen(strcat(filename, ".txt"), "w");
	   	Q->timeStamp = time;
	}
	fprintf(Q->File, str);
	if (time - Q->timeStamp >= 300*1000000LL) {
		fflush(Q->File);
		Q->timeStamp = time;
	}
}


void (writeRow)(world W, char *format, ...)
{
	char str[1024], filename[1024];
	va_list ap;
	va_start(ap, format);
	vsprintf(filename, format, ap);
	double f = va_arg(ap, double);
	if (!NAN2INT(f)) {
		int n = sprintf(str, "%g", f);
		f = va_arg(ap, double);
		while (!NAN2INT(f)) {
			n += sprintf(str+n, "\t%g", f);
			f = va_arg(ap, double);
		}
		sprintf(str+n, "\n");
		writeTxt(W, filename, str);
	}
}


void (writeSpectrum)(world W, int N, float LMin, float LMax, char *format, ...)
{
	static queue staticQ;
	char filename[1024];
	double d;
	float dL = LMin * LMax / (N * W->dt);
	int M = 0, nMin = floor(LMin/dL), nMax = ceil(LMax/dL);
	va_list ap;

	va_start(ap, format);
	vsprintf(filename, format, ap);
	do d = va_arg(ap, double), M++;
	while (!NAN2INT(d));
	va_start(ap, format);
	vsprintf(filename, format, ap);
	queue Q = searchQueue(W, filename, &staticQ);

	if (!Q->Spectrum) {
		source S = W->SE && W->N==W->SE->N ? W->SE : W->SH && W->N==W->SH->N ? W->SH : 0;
		for (world w = W; !S && w->S && w->S->W; S = w->SE && W->N==w->SE->N ? w->SE : w->SH && W->N==w->SH->N ? w->SH : 0) w = W->S->W;

		MALLOC(Q->Spectrum, 0, M);
		MALLOC(Q->Spectrum[0], 0, M*(1-nMin+nMax));
		for (int m=1; m<M; m++) Q->Spectrum[m] = Q->Spectrum[0]+m*(1-nMin+nMax);

		if (S && W->N) {
			for (int n=0; n<=nMax-nMin; n++) {
				float f = (nMax-n) / (dL * nMin * nMax);
				for (int m=1; m<W->N; m++) Q->Spectrum[0][n] += S->waveForm(W, S, W->dt*m, S->phase) * cexpf(2*PI*I*f*W->dt*m);
				if (cabsf(Q->Spectrum[0][n]) < 1E-3) Q->Spectrum[0][n] = 0;
				else if (W->sourceType == 2) Q->Spectrum[0][n] *= sqrtf(f * S->wavelength);
				else if (W->sourceType == 3) Q->Spectrum[0][n] *= f * S->wavelength;
			}
		}
	}

	for (int m=1; m<M; m++) {
		float t = W->t / (dL * nMin * nMax);
		d = va_arg(ap, double);
		#pragma omp parallel for
		for (int n=0; n<=nMax-nMin; n++) Q->Spectrum[m][n] += (float) d * cexpf(2*PI*I*t*(nMax-n));
	}

	if (Q->N == N) {
		Q->File = fopen(strcat(filename, ".txt"), "w");
		for (int n=0; n<=nMax-nMin; n++) {
			fprintf(Q->File, "%g", (dL * nMin * nMax) / (nMax-n));
			for (int m=1; m<M; m++) fprintf(Q->File, "\t%g", Q->Spectrum[0][n] ? cabsf(Q->Spectrum[m][n]/Q->Spectrum[0][n]) : 0);
			for (int m=1; m<M; m++) fprintf(Q->File, "\t%g", Q->Spectrum[0][n] ? cargf(Q->Spectrum[m][n]/Q->Spectrum[0][n])/PI : 0);
			fprintf(Q->File, "\n");
		}
		fclose(Q->File);
	}
}


float (poyntingX)(world W, float x, float yMin, float yMax, float zMin, float zMax, ...)
{
	float Sum=0;
	int i = W->xtoi(W, x<0 && W->xMinSur==SYM ? -x : x);
	int jMin = yMin==yMax ? MAX(W->jMin, W->ytoj(W, yMin)-1) : MAX(W->jMin, W->ytoj(W, yMin));
	int jMax = yMin==yMax ? MIN(W->jMax, W->ytoj(W, yMax)+1) : MIN(W->jMax, W->ytoj(W, yMax));
	int kMin = zMin==zMax ? MAX(W->kMin, W->ztok(W, zMin)-1) : MAX(W->kMin, W->ztok(W, zMin));
	int kMax = zMin==zMax ? MIN(W->kMax, W->ztok(W, zMax)+1) : MIN(W->kMax, W->ztok(W, zMax));

	#pragma omp parallel for reduction(+:Sum)
	for (int j=jMin+1; j<=jMax; j++) {
		float sum=0;
		sum += 0.5 * W->CE->dz[kMin] * W->E->y[i][j][kMin] * (W->H->z[i][j][kMin] + W->H->z[i+1][j][kMin]);
		sum += 0.5 * W->CE->dz[kMax] * W->E->y[i][j][kMax] * (W->H->z[i][j][kMax] + W->H->z[i+1][j][kMax]);
		for (int k=kMin+1; k<kMax; k++) sum += W->CE->dz[k] * W->E->y[i][j][k] * (W->H->z[i][j][k] + W->H->z[i+1][j][k]);
		Sum += W->CH->dy[j] * sum;
	}
	#pragma omp parallel for reduction(+:Sum)
	for (int k=kMin+1; k<=kMax; k++) {
		float sum=0;
		sum -= 0.5 * W->CE->dy[jMin] * W->E->z[i][jMin][k] * (W->H->y[i][jMin][k] + W->H->y[i+1][jMin][k]);
		sum -= 0.5 * W->CE->dy[jMax] * W->E->z[i][jMax][k] * (W->H->y[i][jMax][k] + W->H->y[i+1][jMax][k]);
		for (int j=jMin+1; j<jMax; j++) sum -= W->CE->dy[j] * W->E->z[i][j][k] * (W->H->y[i][j][k] + W->H->y[i+1][j][k]);
		Sum += W->CH->dz[k] * sum;
	}

	if (yMin==yMax || W->Dom.y.Min==W->Dom.y.Max) Sum /= W->CH->dy[jMax]+W->CH->dy[jMax-1];
	if (zMin==zMax || W->Dom.z.Min==W->Dom.z.Max) Sum /= W->CH->dz[kMax]+W->CH->dz[kMax-1];
	if (W->yMinSur==SYM) Sum *= 2;
	if (W->zMinSur==SYM) Sum *= 2;
	return (x<0 && W->xMinSur==SYM ? -0.5 : 0.5) * Sum;
}


float (poyntingY)(world W, float y, float xMin, float xMax, float zMin, float zMax, ...)
{
	float Sum=0;
	int j = W->ytoj(W, y<0 && W->yMinSur==SYM ? -y : y);
	int iMin = xMin==xMax ? MAX(W->iMin, W->xtoi(W, xMin)-1) : MAX(W->iMin, W->xtoi(W, xMin));
	int iMax = xMin==xMax ? MIN(W->iMax, W->xtoi(W, xMax)+1) : MIN(W->iMax, W->xtoi(W, xMax));
	int kMin = zMin==zMax ? MAX(W->kMin, W->ztok(W, zMin)-1) : MAX(W->kMin, W->ztok(W, zMin));
	int kMax = zMin==zMax ? MIN(W->kMax, W->ztok(W, zMax)+1) : MIN(W->kMax, W->ztok(W, zMax));

	#pragma omp parallel for reduction(+:Sum)
	for (int k=kMin+1; k<=kMax; k++) {
		float sum=0;
		sum += 0.5 * W->CE->dx[iMin] * W->E->z[iMin][j][k] * (W->H->x[iMin][j][k] + W->H->x[iMin][j+1][k]);
		sum += 0.5 * W->CE->dx[iMax] * W->E->z[iMax][j][k] * (W->H->x[iMax][j][k] + W->H->x[iMax][j+1][k]);
		for (int i=iMin+1; i<iMax; i++) sum += W->CE->dx[i] * W->E->z[i][j][k] * (W->H->x[i][j][k] + W->H->x[i][j+1][k]);
		Sum += W->CH->dz[k] * sum;
	}
	#pragma omp parallel for reduction(+:Sum)
	for (int i=iMin+1; i<=iMax; i++) {
		float sum=0;
		sum -= 0.5 * W->CE->dz[kMin] * W->E->x[i][j][kMin] * (W->H->z[i][j][kMin] + W->H->z[i][j+1][kMin]);
		sum -= 0.5 * W->CE->dz[kMax] * W->E->x[i][j][kMax] * (W->H->z[i][j][kMax] + W->H->z[i][j+1][kMax]);
		for (int k=kMin+1; k<kMax; k++) sum -= W->CE->dz[k] * W->E->x[i][j][k] * (W->H->z[i][j][k] + W->H->z[i][j+1][k]);
		Sum += W->CH->dx[i] * sum;
	}

	if (xMin==xMax || W->Dom.x.Min==W->Dom.x.Max) Sum /= W->CH->dx[iMax]+W->CH->dx[iMax-1];
	if (zMin==zMax || W->Dom.z.Min==W->Dom.z.Max) Sum /= W->CH->dz[kMax]+W->CH->dz[kMax-1];
	if (W->xMinSur==SYM) Sum *= 2;
	if (W->zMinSur==SYM) Sum *= 2;
	return (y<0 && W->yMinSur==SYM ? -0.5 : 0.5) * Sum;
}


float (poyntingZ)(world W, float z, float xMin, float xMax, float yMin, float yMax, ...)
{
	float Sum=0;
	int k = W->ztok(W, z<0 && W->zMinSur==SYM ? -z : z);
	int iMin = xMin==xMax ? MAX(W->iMin, W->xtoi(W, xMin)-1) : MAX(W->iMin, W->xtoi(W, xMin));
	int iMax = xMin==xMax ? MIN(W->iMax, W->xtoi(W, xMax)+1) : MIN(W->iMax, W->xtoi(W, xMax));
	int jMin = yMin==yMax ? MAX(W->jMin, W->ytoj(W, yMin)-1) : MAX(W->jMin, W->ytoj(W, yMin));
	int jMax = yMin==yMax ? MIN(W->jMax, W->ytoj(W, yMax)+1) : MIN(W->jMax, W->ytoj(W, yMax));

	#pragma omp parallel for reduction(+:Sum)
	for (int i=iMin+1; i<=iMax; i++) {
		float sum=0;
		sum += 0.5 * W->CE->dy[jMin] * W->E->x[i][jMin][k] * (W->H->y[i][jMin][k] + W->H->y[i][jMin][k+1]);
		sum += 0.5 * W->CE->dy[jMax] * W->E->x[i][jMax][k] * (W->H->y[i][jMax][k] + W->H->y[i][jMax][k+1]);
		for (int j=jMin+1; j<jMax; j++) sum += W->CE->dy[j] * W->E->x[i][j][k] * (W->H->y[i][j][k] + W->H->y[i][j][k+1]);
		Sum += W->CH->dx[i] * sum;
	}
	#pragma omp parallel for reduction(+:Sum)
	for (int j=jMin+1; j<=jMax; j++) {
		float sum=0;
		sum -= 0.5 * W->CE->dx[iMin] * W->E->y[iMin][j][k] * (W->H->x[iMin][j][k] + W->H->x[iMin][j][k+1]);
		sum -= 0.5 * W->CE->dx[iMax] * W->E->y[iMax][j][k] * (W->H->x[iMax][j][k] + W->H->x[iMax][j][k+1]);
		for (int i=iMin+1; i<iMax; i++) sum -= W->CE->dx[i] * W->E->y[i][j][k] * (W->H->x[i][j][k] + W->H->x[i][j][k+1]);
		Sum += W->CH->dy[j] * sum;
	}

	if (xMin==xMax || W->Dom.x.Min==W->Dom.x.Max) Sum /= W->CH->dx[iMax]+W->CH->dx[iMax-1];
	if (yMin==yMax || W->Dom.y.Min==W->Dom.y.Max) Sum /= W->CH->dy[jMax]+W->CH->dy[jMax-1];
	if (W->xMinSur==SYM) Sum *= 2;
	if (W->yMinSur==SYM) Sum *= 2;
	return (z<0 && W->zMinSur==SYM ? -0.5 : 0.5) * Sum;
}


float poyntingOut(world W, float xMin, float xMax, float yMin, float yMax, float zMin, float zMax)
{
	float Sx=0, Sy=0, Sz=0;

	if (xMin > W->xMin || xMax < W->xMax || W->xMinSur != BBC) {
		Sx = (poyntingX)(W, xMax, yMin, yMax, zMin, zMax);
		if (xMin < 0 && W->xMinSur==SYM) Sx *= 2;
		else Sx -= (poyntingX)(W, xMin, yMin, yMax, zMin, zMax);
	}
	if (yMin > W->yMin || yMax < W->yMax || W->yMinSur != BBC) {
		Sy = (poyntingY)(W, yMax, xMin, xMax, zMin, zMax);
		if (yMin < 0 && W->yMinSur==SYM) Sy *= 2;
		else Sy -= (poyntingY)(W, yMin, xMin, xMax, zMin, zMax);
	}
	if (zMin > W->zMin || zMax < W->zMax || W->zMinSur != BBC) {
		Sz = (poyntingZ)(W, zMax, xMin, xMax, yMin, yMax);
		if (zMin < 0 && W->zMinSur==SYM) Sz *= 2;
		else Sz -= (poyntingZ)(W, zMin, xMin, xMax, yMin, yMax);
	}

	return Sx + Sy + Sz;
}


static int objcmp(object O1, object O2)
{
	if (O1.S != O2.S || O1.p != O2.p) return 0;
	for (int i=0; i<3; i++) for (int j=0; j<3; j++) if (O1.f[i][j] != O2.f[i][j]) return 0;
	return 1;
}


static int auxiliaryFieldRequired(matter M)
{
	for (int i=0; i<3; i++) for (int j=0; j<6; j++) if (M.d[i][j]) return 2;
	if (M.P[0].omega || M.P[0].gamma) return 1;
	return 0;
}


static float partialAbsorption(world W, int n, int i, int j, float x, float y, float z)
{
	float Sum=0, xsum=0, ysum=0, zsum=0;
	int kMinJ = MAX(W->kMIN, W->CE->kMin[n]), kMaxJ = MIN(W->kMAX, W->CE->kMax[n]);

	for (int k=kMinJ; k<=kMaxJ; k++) {
		if (W->CE->x[i][j][k]) xsum += W->CE->dz[k] * W->E->x[i][j][k] * jx(W, n, i, j, k);
		if (W->CE->y[i][j][k]) ysum += W->CE->dz[k] * W->E->y[i][j][k] * jy(W, n, i, j, k);
		if (W->CE->z[i][j][k]) zsum += W->CH->dz[k] * W->E->z[i][j][k] * jz(W, n, i, j, k);
	}
	if (kMinJ == W->kMIN) {
		if (W->CE->x[i][j][W->kMIN]) xsum -= W->CE->dz[W->kMIN] * W->E->x[i][j][W->kMIN] * jx(W, n, i, j, W->kMIN) * 0.5;
		if (W->CE->y[i][j][W->kMIN]) ysum -= W->CE->dz[W->kMIN] * W->E->y[i][j][W->kMIN] * jy(W, n, i, j, W->kMIN) * 0.5;
		if (W->CE->z[i][j][W->kMIN]) zsum -= W->CH->dz[W->kMIN] * W->E->z[i][j][W->kMIN] * jz(W, n, i, j, W->kMIN);
	}
	if (kMaxJ == W->kMAX) {
		if (W->CE->x[i][j][W->kMAX]) xsum -= W->CE->dz[W->kMAX] * W->E->x[i][j][W->kMAX] * jx(W, n, i, j, W->kMAX) * 0.5;
		if (W->CE->y[i][j][W->kMAX]) ysum -= W->CE->dz[W->kMAX] * W->E->y[i][j][W->kMAX] * jy(W, n, i, j, W->kMAX) * 0.5;
	}
	Sum += x * xsum * W->CH->dx[i] * W->CE->dy[j];
	Sum += y * ysum * W->CE->dx[i] * W->CH->dy[j];
	Sum += z * zsum * W->CE->dx[i] * W->CE->dy[j];
	return Sum;
}


float objectAbsorption(world W, object O)
{
	float Sum=0;
	int n=0;
	for (int o=0; !objcmp(O, W->CE->O[o]); o++) if (auxiliaryFieldRequired(W->CE->M[o])) n++;
	int iMinJ = MAX(W->iMIN, W->CE->iMin[n]), iMaxJ = MIN(W->iMAX, W->CE->iMax[n]);
	int jMinJ = MAX(W->jMIN, W->CE->jMin[n]), jMaxJ = MIN(W->jMAX, W->CE->jMax[n]);

	#pragma omp parallel for reduction(+:Sum)
	for (int i=iMinJ; i<=iMaxJ; i++) {
		for (int j=jMinJ; j<=jMaxJ; j++) Sum += partialAbsorption(W, n, i, j, 1, 1, 1);
		if (jMinJ == W->jMIN) Sum += partialAbsorption(W, n, i, W->jMIN,-0.5,-1,-0.5);
		if (jMaxJ == W->jMAX) Sum += partialAbsorption(W, n, i, W->jMAX,-0.5, 0,-0.5);
	}
	if (iMinJ == W->iMIN) {
		#pragma omp parallel for reduction(+:Sum)
		for (int j=jMinJ; j<=jMaxJ; j++) Sum += partialAbsorption(W, n, W->iMIN, j,-1,-0.5,-0.5);
		if (jMinJ == W->jMIN) Sum += partialAbsorption(W, n, W->iMIN, W->jMIN, 0.5, 0.5, 0.25);
		if (jMaxJ == W->jMAX) Sum += partialAbsorption(W, n, W->iMIN, W->jMAX, 0.5, 0.0, 0.25);
	}
	if (iMaxJ == W->iMAX) {
		#pragma omp parallel for reduction(+:Sum)
		for (int j=jMinJ; j<=jMaxJ; j++) Sum += partialAbsorption(W, n, W->iMAX, j, 0,-0.5,-0.5);
		if (jMinJ == W->jMIN) Sum += partialAbsorption(W, n, W->iMAX, W->jMIN, 0.0, 0.5, 0.25);
		if (jMaxJ == W->jMAX) Sum += partialAbsorption(W, n, W->iMAX, W->jMAX, 0.0, 0.0, 0.25);
	}

	if (W->Dom.x.Min==W->Dom.x.Max) Sum /= W->CH->dx[W->iMax]+W->CH->dx[W->iMax-1];
	if (W->Dom.y.Min==W->Dom.y.Max) Sum /= W->CH->dy[W->jMax]+W->CH->dy[W->jMax-1];
	if (W->Dom.z.Min==W->Dom.z.Max) Sum /= W->CH->dz[W->kMax]+W->CH->dz[W->kMax-1];
	if (W->xMinSur==SYM) Sum *= 2;
	if (W->yMinSur==SYM) Sum *= 2;
	if (W->zMinSur==SYM) Sum *= 2;
	return Sum / W->dt;
}


void exec(char *format, ...)
{
	char str[1024];
	va_list ap;
	va_start(ap, format);
	vsprintf(str, format, ap);
	system(str);
}
