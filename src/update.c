/**********************************************************************
 * ALiS is Copyright (C) 2009-2019 by Ho-Seok Ee <hsee@kongju.ac.kr>. *
 * Redistribution and use with or without modification, are permitted *
 * under the terms of the Artistic License version 2.                 *
 **********************************************************************/

#include <sys/time.h>
#include "alis_c.h"
#define FOR3D(i,iMin,iMax,j,jMin,jMax,k,kMin,kMax) \
_Pragma("omp parallel for") for (int i=iMin; i<=iMax; i++) for (int j=jMin; j<=jMax; j++) \
_Pragma("GCC ivdep") for (int k=kMin; k<=kMax; k++)


int (timer)(long int n, long int N, char *str, ...)
{
	static long int initialStep, previousStep;
	static unsigned long long int initialTime, previousTime;
	unsigned long long int currentTime, elapsedTime, estimatedTime;
	struct timeval now;
	gettimeofday(&now, 0);
	currentTime = now.tv_sec * 1000000LL + now.tv_usec;

	if (!n || !initialStep) {
		initialStep = previousStep = n;
		initialTime = previousTime = currentTime;
		estimatedTime = elapsedTime = 0LL;
	}
	else if (n>N || currentTime/1000000LL > previousTime/1000000LL) {
		elapsedTime = currentTime - initialTime;
		estimatedTime = elapsedTime + (N-n+1) * (currentTime-previousTime) / (n-previousStep);
		if (estimatedTime - elapsedTime > 3600000000LL) estimatedTime -= estimatedTime%60000000LL;
		previousStep = n;
		previousTime = currentTime;
		if (str[0]) printf("\r%s ", str);
		else printf("\r");
		printf("[%02lld:%02lld:%02lld/%02lld:%02lld:%02lld] [%ld/%ld: %.1f%%]",
			elapsedTime/3600000000LL, elapsedTime%3600000000LL/60000000LL, elapsedTime%60000000LL/1000000LL,
			estimatedTime/3600000000LL, estimatedTime%3600000000LL/60000000LL, estimatedTime%60000000LL/1000000LL,
			(n-1), N, 100.0*(n-1)/N);
		fflush(stdout);
		if (n>N) printf("\n");
	}

	return n<=N;
}


static void updatePx(vfield F, vfield G, coeffs *C, float ***Py, float ***Pz, int Alt, int iMin, int iMax, int jMin, int jMax, int kMin, int kMax)
{
	if (C->y && C->z) FOR3D(j,jMin,jMax,i,iMin,iMax,k,kMin,kMax) {
		F.y[i][j][k] += (Py[i][j][k] = C->xP[i] * Py[i][j][k] + C->xPy[i] * (G.z[i][j][k] - G.z[i+Alt][j][k])) * C->y[i][j][k];
		F.z[i][j][k] += (Pz[i][j][k] = C->xP[i] * Pz[i][j][k] + C->xPz[i] * (G.y[i][j][k] - G.y[i+Alt][j][k])) * C->z[i][j][k];
	}
	else FOR3D(j,jMin,jMax,i,iMin,iMax,k,kMin,kMax) {
		F.y[i][j][k] += (Py[i][j][k] = C->xP[i] * Py[i][j][k] + C->xPy[i] * (G.z[i][j][k] - G.z[i+Alt][j][k]));
		F.z[i][j][k] += (Pz[i][j][k] = C->xP[i] * Pz[i][j][k] + C->xPz[i] * (G.y[i][j][k] - G.y[i+Alt][j][k]));
	}
}


static void updatePy(vfield F, vfield G, coeffs *C, float ***Pz, float ***Px, int Alt, int iMin, int iMax, int jMin, int jMax, int kMin, int kMax)
{
	if (C->z && C->x) FOR3D(i,iMin,iMax,j,jMin,jMax,k,kMin,kMax) {
		F.z[i][j][k] += (Pz[i][j][k] = C->yP[j] * Pz[i][j][k] + C->yPz[j] * (G.x[i][j][k] - G.x[i][j+Alt][k])) * C->z[i][j][k];
		F.x[i][j][k] += (Px[i][j][k] = C->yP[j] * Px[i][j][k] + C->yPx[j] * (G.z[i][j][k] - G.z[i][j+Alt][k])) * C->x[i][j][k];
	}
	else FOR3D(i,iMin,iMax,j,jMin,jMax,k,kMin,kMax) {
		F.z[i][j][k] += (Pz[i][j][k] = C->yP[j] * Pz[i][j][k] + C->yPz[j] * (G.x[i][j][k] - G.x[i][j+Alt][k]));
		F.x[i][j][k] += (Px[i][j][k] = C->yP[j] * Px[i][j][k] + C->yPx[j] * (G.z[i][j][k] - G.z[i][j+Alt][k]));
	}
}


static void updatePz(vfield F, vfield G, coeffs *C, float ***Px, float ***Py, int Alt, int iMin, int iMax, int jMin, int jMax, int kMin, int kMax)
{
	if (C->x && C->y) FOR3D(i,iMin,iMax,j,jMin,jMax,k,kMin,kMax) {
		F.x[i][j][k] += (Px[i][j][k] = C->zP[k] * Px[i][j][k] + C->zPx[k] * (G.y[i][j][k] - G.y[i][j][k+Alt])) * C->x[i][j][k];
		F.y[i][j][k] += (Py[i][j][k] = C->zP[k] * Py[i][j][k] + C->zPy[k] * (G.x[i][j][k] - G.x[i][j][k+Alt])) * C->y[i][j][k];
	}
	else FOR3D(i,iMin,iMax,j,jMin,jMax,k,kMin,kMax) {
		F.x[i][j][k] += (Px[i][j][k] = C->zP[k] * Px[i][j][k] + C->zPx[k] * (G.y[i][j][k] - G.y[i][j][k+Alt]));
		F.y[i][j][k] += (Py[i][j][k] = C->zP[k] * Py[i][j][k] + C->zPy[k] * (G.x[i][j][k] - G.x[i][j][k+Alt]));
	}
}


static void updateF(world W, vfield F, vfield G, coeffs *C, int Alt, int iMin, int iMax, int jMin, int jMax, int kMin, int kMax)
{
	for (int n=0; n<C->N; n++) {
		int iMinJ = MAX(iMin, C->iMin[n]), iMaxJ = MIN(iMax, C->iMax[n]);
		int jMinJ = MAX(jMin, C->jMin[n]), jMaxJ = MIN(jMax, C->jMax[n]);
		int kMinJ = MAX(kMin, C->kMin[n]), kMaxJ = MIN(kMax, C->kMax[n]);
		if (C->ax[n] || C->ay[n] || C->az[n]) FOR3D(i,iMinJ,iMaxJ,j,jMinJ,jMaxJ,k,kMinJ,kMaxJ) {
			F.x[i][j][k] *= (1 - C->ax[n] * C->Jx[n][i][j][k]) / (1 + C->ax[n] * C->Jx[n][i][j][k]);
			F.y[i][j][k] *= (1 - C->ay[n] * C->Jy[n][i][j][k]) / (1 + C->ay[n] * C->Jy[n][i][j][k]);
			F.z[i][j][k] *= (1 - C->az[n] * C->Jz[n][i][j][k]) / (1 + C->az[n] * C->Jz[n][i][j][k]);
		}
	}
	for (int n=0; n<C->N; n++) {
		int iMinJ = MAX(iMin, C->iMin[n]), iMaxJ = MIN(iMax, C->iMax[n]);
		int jMinJ = MAX(jMin, C->jMin[n]), jMaxJ = MIN(jMax, C->jMax[n]);
		int kMinJ = MAX(kMin, C->kMin[n]), kMaxJ = MIN(kMax, C->kMax[n]);
		for (int p=0; p<C->NK[n]; p++) FOR3D(i,iMinJ,iMaxJ,j,jMinJ,jMaxJ,k,kMinJ,kMaxJ) {
			F.x[i][j][k] -= F.Jx[n][p][i][j][k];
			F.y[i][j][k] -= F.Jy[n][p][i][j][k];
			F.z[i][j][k] -= F.Jz[n][p][i][j][k];
		}
	}
	if (C->x && C->y && C->z) FOR3D(i,iMin,iMax,j,jMin,jMax,k,kMin,kMax) {
		F.x[i][j][k] += (C->dtdz[k] * (G.y[i][j][k] - G.y[i][j][k+Alt]) + C->dtdy[j] * (G.z[i][j+Alt][k] - G.z[i][j][k])) * C->x[i][j][k];
		F.y[i][j][k] += (C->dtdx[i] * (G.z[i][j][k] - G.z[i+Alt][j][k]) + C->dtdz[k] * (G.x[i][j][k+Alt] - G.x[i][j][k])) * C->y[i][j][k];
		F.z[i][j][k] += (C->dtdy[j] * (G.x[i][j][k] - G.x[i][j+Alt][k]) + C->dtdx[i] * (G.y[i+Alt][j][k] - G.y[i][j][k])) * C->z[i][j][k];
	}
	else FOR3D(i,iMin,iMax,j,jMin,jMax,k,kMin,kMax) {
		F.x[i][j][k] += (C->dtdz[k] * (G.y[i][j][k] - G.y[i][j][k+Alt]) + C->dtdy[j] * (G.z[i][j+Alt][k] - G.z[i][j][k]));
		F.y[i][j][k] += (C->dtdx[i] * (G.z[i][j][k] - G.z[i+Alt][j][k]) + C->dtdz[k] * (G.x[i][j][k+Alt] - G.x[i][j][k]));
		F.z[i][j][k] += (C->dtdy[j] * (G.x[i][j][k] - G.x[i][j+Alt][k]) + C->dtdx[i] * (G.y[i+Alt][j][k] - G.y[i][j][k]));
	}
	if (W->iMIN < W->iMin) updatePx(F, G, C, F.xMinPy, F.xMinPz, Alt, W->iMIN-(Alt-1)/2, W->iMin-(Alt+3)/2, jMin, jMax, kMin, kMax);
	if (W->jMIN < W->jMin) updatePy(F, G, C, F.yMinPz, F.yMinPx, Alt, iMin, iMax, W->jMIN-(Alt-1)/2, W->jMin-(Alt+3)/2, kMin, kMax);
	if (W->kMIN < W->kMin) updatePz(F, G, C, F.zMinPx, F.zMinPy, Alt, iMin, iMax, jMin, jMax, W->kMIN-(Alt-1)/2, W->kMin-(Alt+3)/2);
	if (W->iMAX > W->iMax) updatePx(F, G, C, F.xMaxPy, F.xMaxPz, Alt, W->iMax+2, W->iMAX, jMin, jMax, kMin, kMax);
	if (W->jMAX > W->jMax) updatePy(F, G, C, F.yMaxPz, F.yMaxPx, Alt, iMin, iMax, W->jMax+2, W->jMAX, kMin, kMax);
	if (W->kMAX > W->kMax) updatePz(F, G, C, F.zMaxPx, F.zMaxPy, Alt, iMin, iMax, jMin, jMax, W->kMax+2, W->kMAX);
}


static void updateJ(world W, vfield F, coeffs *C)
{
	for (int n=0; n<C->N; n++) {
		int iMinJ = MAX(W->iMIN, C->iMin[n]), iMaxJ = MIN(W->iMAX, C->iMax[n]);
		int jMinJ = MAX(W->jMIN, C->jMin[n]), jMaxJ = MIN(W->jMAX, C->jMax[n]);
		int kMinJ = MAX(W->kMIN, C->kMin[n]), kMaxJ = MIN(W->kMAX, C->kMax[n]);
		for (int p=0; p<C->NJ[n]; p++) FOR3D(i,iMinJ,iMaxJ,j,jMinJ,jMaxJ,k,kMinJ,kMaxJ) {
			F.Jx[n][p][i][j][k] = C->JJx[n][p] * F.Jx[n][p][i][j][k] + C->JEx[n][p] * C->Jx[n][i][j][k] * F.x[i][j][k];
			F.Jy[n][p][i][j][k] = C->JJy[n][p] * F.Jy[n][p][i][j][k] + C->JEy[n][p] * C->Jy[n][i][j][k] * F.y[i][j][k];
			F.Jz[n][p][i][j][k] = C->JJz[n][p] * F.Jz[n][p][i][j][k] + C->JEz[n][p] * C->Jz[n][i][j][k] * F.z[i][j][k];
		}
		for (int p=C->NJ[n]; p<C->NK[n]; p++) FOR3D(i,iMinJ,iMaxJ,j,jMinJ,jMaxJ,k,kMinJ,kMaxJ) {
			F.Jx[n][p][i][j][k] = C->JJx[n][p] * F.Jx[n][p][i][j][k]
								+ C->JKx[n][p] * F.Kx[n][p][i][j][k] + C->JEx[n][p] * C->Jx[n][i][j][k] * F.x[i][j][k];
			F.Kx[n][p][i][j][k] = C->KJx[n][p] * F.Jx[n][p][i][j][k]
								+ C->KKx[n][p] * F.Kx[n][p][i][j][k] + C->KEx[n][p] * C->Jx[n][i][j][k] * F.x[i][j][k];
		}
		for (int p=C->NJ[n]; p<C->NK[n]; p++) FOR3D(i,iMinJ,iMaxJ,j,jMinJ,jMaxJ,k,kMinJ,kMaxJ) {
			F.Jy[n][p][i][j][k] = C->JJy[n][p] * F.Jy[n][p][i][j][k]
								+ C->JKy[n][p] * F.Ky[n][p][i][j][k] + C->JEy[n][p] * C->Jy[n][i][j][k] * F.y[i][j][k];
			F.Ky[n][p][i][j][k] = C->KJy[n][p] * F.Jy[n][p][i][j][k]
								+ C->KKy[n][p] * F.Ky[n][p][i][j][k] + C->KEy[n][p] * C->Jy[n][i][j][k] * F.y[i][j][k];
		}
		for (int p=C->NJ[n]; p<C->NK[n]; p++) FOR3D(i,iMinJ,iMaxJ,j,jMinJ,jMaxJ,k,kMinJ,kMaxJ) {
			F.Jz[n][p][i][j][k] = C->JJz[n][p] * F.Jz[n][p][i][j][k]
								+ C->JKz[n][p] * F.Kz[n][p][i][j][k] + C->JEz[n][p] * C->Jz[n][i][j][k] * F.z[i][j][k];
			F.Kz[n][p][i][j][k] = C->KJz[n][p] * F.Jz[n][p][i][j][k]
								+ C->KKz[n][p] * F.Kz[n][p][i][j][k] + C->KEz[n][p] * C->Jz[n][i][j][k] * F.z[i][j][k];
		}
	}
}


static void updateS(world W, source S, float t)
{
	if (W->n < S->N) {
		if (S->Phase) {
			if (S->Amplitude) FOR3D(i,S->iMin,S->iMax,j,S->jMin,S->jMax,k,S->kMin,S->kMax)
				S->F[i][j][k] += S->Amplitude[i][j][k] * S->waveForm(W, S, t, S->Phase[i][j][k]);
			else FOR3D(i,S->iMin,S->iMax,j,S->jMin,S->jMax,k,S->kMin,S->kMax)
				S->F[i][j][k] += S->amplitude * S->waveForm(W, S, t, S->Phase[i][j][k]);
		}
		else {
			float f = S->waveForm(W, S, t, S->phase);
			if (S->Amplitude) FOR3D(i,S->iMin,S->iMax,j,S->jMin,S->jMax,k,S->kMin,S->kMax)
				S->F[i][j][k] += S->Amplitude[i][j][k] * f;
			else FOR3D(i,S->iMin,S->iMax,j,S->jMin,S->jMax,k,S->kMin,S->kMax)
				S->F[i][j][k] += S->amplitude * f;
		}
	}
	if (S->S) updateS(W, S->S, t);
}


void (copyBoundary)(float ***F, float ***G, float cosk, float sink, int iMin, int iMax, int jMin, int jMax, int kMin, int kMax, int iAlt, int jAlt, int kAlt)
{
	if (G) FOR3D(i,iMin,iMax,j,jMin,jMax,k,kMin,kMax) {
		F[i][j][k] = cosk * F[i+iAlt][j+jAlt][k+kAlt] + sink * G[i+iAlt][j+jAlt][k+kAlt];
		G[i][j][k] = cosk * G[i+iAlt][j+jAlt][k+kAlt] - sink * F[i+iAlt][j+jAlt][k+kAlt];
	}
	else FOR3D(i,iMin,iMax,j,jMin,jMax,k,kMin,kMax) {
		F[i][j][k] = cosk * F[i+iAlt][j+jAlt][k+kAlt];
	}
}


static void update(float ***F, float ***G, float ***C, float dtds, int iMin, int iMax, int jMin, int jMax, int kMin, int kMax, int iAlt, int jAlt, int kAlt)
{
	if (C) FOR3D(i,iMin,iMax,j,jMin,jMax,k,kMin,kMax) {
		F[i][j][k] += dtds * G[i+iAlt][j+jAlt][k+kAlt] * C[i][j][k];
	}
	else FOR3D(i,iMin,iMax,j,jMin,jMax,k,kMin,kMax) {
		F[i][j][k] += dtds * G[i+iAlt][j+jAlt][k+kAlt];
	}
}


static void updateTFSFE(world W, tfsf S)
{
	updateE(S->W);

	for (int n=0; n<=W->complexField; n++) {
		if (S->iMin > W->iMIN) {
			update(W->E[n].z, S->W->H[n].y, W->CE->z, -W->CE->dtdx[S->iMin], S->iMin, S->iMin, S->jMin, S->jMax, S->kMIN, S->kMAX, 0, 0, 0);
			update(W->E[n].y, S->W->H[n].z, W->CE->y, +W->CE->dtdx[S->iMin], S->iMin, S->iMin, S->jMIN, S->jMAX, S->kMin, S->kMax, 0, 0, 0);
		}
		if (S->jMin > W->jMIN) {
			update(W->E[n].x, S->W->H[n].z, W->CE->x, -W->CE->dtdy[S->jMin], S->iMIN, S->iMAX, S->jMin, S->jMin, S->kMin, S->kMax, 0, 0, 0);
			update(W->E[n].z, S->W->H[n].x, W->CE->z, +W->CE->dtdy[S->jMin], S->iMin, S->iMax, S->jMin, S->jMin, S->kMIN, S->kMAX, 0, 0, 0);
		}
		if (S->kMin > W->kMIN) {
			update(W->E[n].y, S->W->H[n].x, W->CE->y, -W->CE->dtdz[S->kMin], S->iMin, S->iMax, S->jMIN, S->jMAX, S->kMin, S->kMin, 0, 0, 0);
			update(W->E[n].x, S->W->H[n].y, W->CE->x, +W->CE->dtdz[S->kMin], S->iMIN, S->iMAX, S->jMin, S->jMax, S->kMin, S->kMin, 0, 0, 0);
		}
		if (S->iMax < W->iMAX) {
			update(W->E[n].z, S->W->H[n].y, W->CE->z, +W->CE->dtdx[S->iMax], S->iMax, S->iMax, S->jMin, S->jMax, S->kMIN, S->kMAX, 1, 0, 0);
			update(W->E[n].y, S->W->H[n].z, W->CE->y, -W->CE->dtdx[S->iMax], S->iMax, S->iMax, S->jMIN, S->jMAX, S->kMin, S->kMax, 1, 0, 0);
		}
		if (S->jMax < W->jMAX) {
			update(W->E[n].x, S->W->H[n].z, W->CE->x, +W->CE->dtdy[S->jMax], S->iMIN, S->iMAX, S->jMax, S->jMax, S->kMin, S->kMax, 0, 1, 0);
			update(W->E[n].z, S->W->H[n].x, W->CE->z, -W->CE->dtdy[S->jMax], S->iMin, S->iMax, S->jMax, S->jMax, S->kMIN, S->kMAX, 0, 1, 0);
		}
		if (S->kMax < W->kMAX) {
			update(W->E[n].y, S->W->H[n].x, W->CE->y, +W->CE->dtdz[S->kMax], S->iMin, S->iMax, S->jMIN, S->jMAX, S->kMax, S->kMax, 0, 0, 1);
			update(W->E[n].x, S->W->H[n].y, W->CE->x, -W->CE->dtdz[S->kMax], S->iMIN, S->iMAX, S->jMin, S->jMax, S->kMax, S->kMax, 0, 0, 1);
		}
	}

	if (S->S) updateTFSFE(W, S->S);
}


static void updateTFSFH(world W, tfsf S)
{
	updateH(S->W);

	for (int n=0; n<=W->complexField; n++) {
		if (S->iMin > W->iMIN) {
			update(W->H[n].z, S->W->E[n].y, W->CH->z, +W->CH->dtdx[S->iMin], S->iMin, S->iMin, S->jMIN, S->jMAX, S->kMin, S->kMax, 0, 0, 0);
			update(W->H[n].y, S->W->E[n].z, W->CH->y, -W->CH->dtdx[S->iMin], S->iMin, S->iMin, S->jMin, S->jMax, S->kMIN, S->kMAX, 0, 0, 0);
		}
		if (S->jMin > W->jMIN) {
			update(W->H[n].x, S->W->E[n].z, W->CH->x, +W->CH->dtdy[S->jMin], S->iMin, S->iMax, S->jMin, S->jMin, S->kMIN, S->kMAX, 0, 0, 0);
			update(W->H[n].z, S->W->E[n].x, W->CH->z, -W->CH->dtdy[S->jMin], S->iMIN, S->iMAX, S->jMin, S->jMin, S->kMin, S->kMax, 0, 0, 0);
		}
		if (S->kMin > W->kMIN) {
			update(W->H[n].y, S->W->E[n].x, W->CH->y, +W->CH->dtdz[S->kMin], S->iMIN, S->iMAX, S->jMin, S->jMax, S->kMin, S->kMin, 0, 0, 0);
			update(W->H[n].x, S->W->E[n].y, W->CH->x, -W->CH->dtdz[S->kMin], S->iMin, S->iMax, S->jMIN, S->jMAX, S->kMin, S->kMin, 0, 0, 0);
		}
		if (S->iMax < W->iMAX) {
			update(W->H[n].z, S->W->E[n].y, W->CH->z, -W->CH->dtdx[S->iMax+1], S->iMax+1, S->iMax+1, S->jMIN, S->jMAX, S->kMin, S->kMax,-1, 0, 0);
			update(W->H[n].y, S->W->E[n].z, W->CH->y, +W->CH->dtdx[S->iMax+1], S->iMax+1, S->iMax+1, S->jMin, S->jMax, S->kMIN, S->kMAX,-1, 0, 0);
		}
		if (S->jMax < W->jMAX) {
			update(W->H[n].x, S->W->E[n].z, W->CH->x, -W->CH->dtdy[S->jMax+1], S->iMin, S->iMax, S->jMax+1, S->jMax+1, S->kMIN, S->kMAX, 0,-1, 0);
			update(W->H[n].z, S->W->E[n].x, W->CH->z, +W->CH->dtdy[S->jMax+1], S->iMIN, S->iMAX, S->jMax+1, S->jMax+1, S->kMin, S->kMax, 0,-1, 0);
		}
		if (S->kMax < W->kMAX) {
			update(W->H[n].y, S->W->E[n].x, W->CH->y, -W->CH->dtdz[S->kMax+1], S->iMIN, S->iMAX, S->jMin, S->jMax, S->kMax+1, S->kMax+1, 0, 0,-1);
			update(W->H[n].x, S->W->E[n].y, W->CH->x, +W->CH->dtdz[S->kMax+1], S->iMin, S->iMax, S->jMIN, S->jMAX, S->kMax+1, S->kMax+1, 0, 0,-1);
		}
	}

	if (S->S) updateTFSFH(W, S->S);
}


void updateE(world W)
{
	for (int m=0; m<=W->complexField; m++) {
		updateF(W, W->E[m], W->H[m], W->CE, 1, W->iMIN, W->iMAX, W->jMIN, W->jMAX, W->kMIN, W->kMAX);
		updateJ(W, W->H[m], W->CH);
	}
	if (W->S) updateTFSFE(W, W->S);
	if (W->SE) updateS(W, W->SE, W->H->t);
	int iAlt = (W->xMinSur==PBC || W->xMinSur==BBC ? W->iNUM-1 : 1);
	int jAlt = (W->yMinSur==PBC || W->yMinSur==BBC ? W->jNUM-1 : 1);
	int kAlt = (W->zMinSur==PBC || W->zMinSur==BBC ? W->kNUM-1 : 1);
	copyBoundary(W, E, x, W->coskx, -W->sinkx, W->iMAX+1, W->iMAX+1, W->jMIN, W->jMAX, W->kMIN, W->kMAX,-iAlt, 0, 0);
	copyBoundary(W, E, y, W->cosky, -W->sinky, W->iMIN, W->iMAX, W->jMAX+1, W->jMAX+1, W->kMIN, W->kMAX, 0,-jAlt, 0);
	copyBoundary(W, E, z, W->coskz, -W->sinkz, W->iMIN, W->iMAX, W->jMIN, W->jMAX, W->kMAX+1, W->kMAX+1, 0, 0,-kAlt);
	W->t = W->E->t = ++W->n * W->dt;
}


void updateH(world W)
{
	for (int m=0; m<=W->complexField; m++) {
		updateF(W, W->H[m], W->E[m], W->CH, -1, W->iMIN, W->iMAX, W->jMIN, W->jMAX, W->kMIN, W->kMAX);
		updateJ(W, W->E[m], W->CE);
	}
	if (W->S) updateTFSFH(W, W->S);
	if (W->SH) updateS(W, W->SH, W->E->t);
	int iAlt = (W->xMinSur==PBC || W->xMinSur==BBC ? W->iNUM-1 : 1);
	int jAlt = (W->yMinSur==PBC || W->yMinSur==BBC ? W->jNUM-1 : 1);
	int kAlt = (W->zMinSur==PBC || W->zMinSur==BBC ? W->kNUM-1 : 1);
	copyBoundary(W, H, y, W->xMaxSur==PEC?1:W->coskx, -W->sinkx, W->iMAX+1, W->iMAX+1, W->jMIN, W->jMAX+1, W->kMIN, W->kMAX+1,-iAlt, 0, 0);
	copyBoundary(W, H, z, W->xMaxSur==PEC?1:W->coskx, -W->sinkx, W->iMAX+1, W->iMAX+1, W->jMIN, W->jMAX+1, W->kMIN, W->kMAX+1,-iAlt, 0, 0);
	copyBoundary(W, H, z, W->yMaxSur==PEC?1:W->cosky, -W->sinky, W->iMIN, W->iMAX+1, W->jMAX+1, W->jMAX+1, W->kMIN, W->kMAX+1, 0,-jAlt, 0);
	copyBoundary(W, H, x, W->yMaxSur==PEC?1:W->cosky, -W->sinky, W->iMIN, W->iMAX+1, W->jMAX+1, W->jMAX+1, W->kMIN, W->kMAX+1, 0,-jAlt, 0);
	copyBoundary(W, H, x, W->zMaxSur==PEC?1:W->coskz, -W->sinkz, W->iMIN, W->iMAX+1, W->jMIN, W->jMAX+1, W->kMAX+1, W->kMAX+1, 0, 0,-kAlt);
	copyBoundary(W, H, y, W->zMaxSur==PEC?1:W->coskz, -W->sinkz, W->iMIN, W->iMAX+1, W->jMIN, W->jMAX+1, W->kMAX+1, W->kMAX+1, 0, 0,-kAlt);
	copyBoundary(W, H, y, W->coskx, W->sinkx, W->iMIN, W->iMIN, W->jMIN, W->jMAX+1, W->kMIN, W->kMAX+1, iAlt, 0, 0);
	copyBoundary(W, H, z, W->coskx, W->sinkx, W->iMIN, W->iMIN, W->jMIN, W->jMAX+1, W->kMIN, W->kMAX+1, iAlt, 0, 0);
	copyBoundary(W, H, z, W->cosky, W->sinky, W->iMIN, W->iMAX+1, W->jMIN, W->jMIN, W->kMIN, W->kMAX+1, 0, jAlt, 0);
	copyBoundary(W, H, x, W->cosky, W->sinky, W->iMIN, W->iMAX+1, W->jMIN, W->jMIN, W->kMIN, W->kMAX+1, 0, jAlt, 0);
	copyBoundary(W, H, x, W->coskz, W->sinkz, W->iMIN, W->iMAX+1, W->jMIN, W->jMAX+1, W->kMIN, W->kMIN, 0, 0, kAlt);
	copyBoundary(W, H, y, W->coskz, W->sinkz, W->iMIN, W->iMAX+1, W->jMIN, W->jMAX+1, W->kMIN, W->kMIN, 0, 0, kAlt);
	W->t = W->H->t = W->E->t + 0.5 * W->dt;
}


#define xx(G) (G.x[(i)][(j)][(k)])
#define yy(G) (G.y[(i)][(j)][(k)])
#define zz(G) (G.z[(i)][(j)][(k)])
#define yx(G) (0.25 * (G.y[(i)][(j)][(k)] + G.y[(i)-1][(j)][(k)] + G.y[(i)][(j)+1][(k)] + G.y[(i)-1][(j)+1][(k)]))
#define zx(G) (0.25 * (G.z[(i)][(j)][(k)] + G.z[(i)-1][(j)][(k)] + G.z[(i)][(j)][(k)+1] + G.z[(i)-1][(j)][(k)+1]))
#define zy(G) (0.25 * (G.z[(i)][(j)][(k)] + G.z[(i)][(j)-1][(k)] + G.z[(i)][(j)][(k)+1] + G.z[(i)][(j)-1][(k)+1]))
#define xy(G) (0.25 * (G.x[(i)][(j)][(k)] + G.x[(i)][(j)-1][(k)] + G.x[(i)+1][(j)][(k)] + G.x[(i)+1][(j)-1][(k)]))
#define xz(G) (0.25 * (G.x[(i)][(j)][(k)] + G.x[(i)][(j)][(k)-1] + G.x[(i)+1][(j)][(k)] + G.x[(i)+1][(j)][(k)-1]))
#define yz(G) (0.25 * (G.y[(i)][(j)][(k)] + G.y[(i)][(j)][(k)-1] + G.y[(i)][(j)+1][(k)] + G.y[(i)][(j)+1][(k)-1]))

#define twm(d,v,J,G) if (C->M[n].d) \
FOR3D(i,iMin,iMax,j,jMin,jMax,k,kMin,kMax) F.v[i][j][k] += 2 * Alt * C->M[n].d * C->v[i][j][k] * C->J[n][i][j][k] * (G);


static void coupleTWM(coeffs *C, vfield F, vfield G1, vfield G2, int n, float Alt)
{
	int iMin = C->iMin[n], iMax = C->iMax[n]-1;
	int jMin = C->jMin[n], jMax = C->jMax[n]-1;
	int kMin = C->kMin[n], kMax = C->kMax[n]-1;

	twm(d[0][0], x, Jx, xx(G1) * xx(G2));
	twm(d[0][1], x, Jx, yx(G1) * yx(G2));
	twm(d[0][2], x, Jx, zx(G1) * zx(G2));
	twm(d[0][3], x, Jx, yx(G1) * zx(G2));
	twm(d[0][3], x, Jx, zx(G1) * yx(G2));
	twm(d[0][4], x, Jx, zx(G1) * xx(G2));
	twm(d[0][4], x, Jx, xx(G1) * zx(G2));
	twm(d[0][5], x, Jx, xx(G1) * yx(G2));
	twm(d[0][5], x, Jx, yx(G1) * xx(G2));
	twm(d[1][0], y, Jy, xy(G1) * xy(G2));
	twm(d[1][1], y, Jy, yy(G1) * yy(G2));
	twm(d[1][2], y, Jy, zy(G1) * zy(G2));
	twm(d[1][3], y, Jy, yy(G1) * zy(G2));
	twm(d[1][3], y, Jy, zy(G1) * yy(G2));
	twm(d[1][4], y, Jy, zy(G1) * xy(G2));
	twm(d[1][4], y, Jy, xy(G1) * zy(G2));
	twm(d[1][5], y, Jy, xy(G1) * yy(G2));
	twm(d[1][5], y, Jy, yy(G1) * xy(G2));
	twm(d[2][0], z, Jz, xz(G1) * xz(G2));
	twm(d[2][1], z, Jz, yz(G1) * yz(G2));
	twm(d[2][2], z, Jz, zz(G1) * zz(G2));
	twm(d[2][3], z, Jz, yz(G1) * zz(G2));
	twm(d[2][3], z, Jz, zz(G1) * yz(G2));
	twm(d[2][4], z, Jz, zz(G1) * xz(G2));
	twm(d[2][4], z, Jz, xz(G1) * zz(G2));
	twm(d[2][5], z, Jz, xz(G1) * yz(G2));
	twm(d[2][5], z, Jz, yz(G1) * xz(G2));
}


void nonLinearSHG(world WS, world WP, float Alt)
{
	for (int n=0; n<WS->CE->N; n++) {
		for (int m=0; m<=WS->complexField; m++) {
			coupleTWM(WS->CE, WS->E[m], WP->E[m], WP->E[m], n, Alt);
		}
	}
}


void nonLinearTWM(world WS, world WP1, world WP2, float Alt)
{
	for (int n=0; n<WS->CE->N; n++) {
		for (int m=0; m<=WS->complexField; m++) {
			coupleTWM(WS->CE, WS->E[m], WP1->E[m], WP2->E[m], n, Alt);
		}
	}
}
