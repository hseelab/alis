/**********************************************************************
 * ALiS is Copyright (C) 2009-2019 by Ho-Seok Ee <hsee@kongju.ac.kr>. *
 * Redistribution and use with or without modification, are permitted *
 * under the terms of the Artistic License version 2.                 *
 **********************************************************************/

#include "alis_c.h"


static int auxiliaryFieldRequired(matter M)
{
	for (int i=0; i<3; i++) for (int j=0; j<6; j++) if (M.d[i][j]) return 2;
	if (M.P[0].w || M.P[0].r) return 1;
	return 0;
}

static float pole(world W, matter M, float w, int p)
{
	if (M.P[p].f || M.P[p].g)
		return 0.5
			*((M.P[p].f * M.P[p].w - M.P[p].g * M.P[p].r + I * M.P[p].g * M.P[p].w) / (M.P[p].w + w + I * M.P[p].r)
			+ (M.P[p].f * M.P[p].w - M.P[p].g * M.P[p].r - I * M.P[p].g * M.P[p].w) / (M.P[p].w - w - I * M.P[p].r));
	if (M.P[p].w) return - sq(M.P[p].w) / (sq(w) + I * w * M.P[p].r);
	return 0;
}

static float dpole(world W, matter M, float w, int p)
{
	if (M.P[p].f || M.P[p].g)
		return 0.5
			*((M.P[p].f * M.P[p].w - M.P[p].g * M.P[p].r + I * M.P[p].g * M.P[p].w) / sq(M.P[p].w - w - I * M.P[p].r)
			- (M.P[p].f * M.P[p].w - M.P[p].g * M.P[p].r - I * M.P[p].g * M.P[p].w) / sq(M.P[p].w + w + I * M.P[p].r));
	if (M.P[p].w) return (2 * w + I * M.P[p].r) * sq(M.P[p].w) / sq(sq(w) + I * w * M.P[p].r);
	return 0;
}


/* Refractive index */

float raweEffx(world W, int i, int j, int k)
{
	if (!W->CE->x) return 1;
	if (i<0 && W->xMinSur==SYM) i=1-i;
	if (j<0 && W->yMinSur==SYM) j=-j;
	if (k<0 && W->zMinSur==SYM) k=-k;
	float f = 1 / W->CE->x[i][j][k];
	for (int n=0, o=0; n<W->CE->N; n++, o++) {
		while (!auxiliaryFieldRequired(W->CE->M[o])) o++;
		if (W->CE->iMin[n]<=i && i<=W->CE->iMax[n] && W->CE->jMin[n]<=j && j<=W->CE->jMax[n] && W->CE->kMin[n]<=k && k<=W->CE->kMax[n]) {
			int m = W->CE->M[o].e.y || W->CE->M[o].e.z ? 1 : 0;
			f -= W->CE->ax[n] * (W->CE->M[o].e.x ? W->CE->M[o].e.x : W->CE->M[o].e.x) * W->CE->Jx[n][i][j][k];
			for (int p=0; p<MAXPOLES; p+=1+2*m)
				f += (pole(W, W->CE->M[o], W->f*W->eV, p+0*m) + W->f*W->eV * dpole(W, W->CE->M[o], W->f*W->eV, p+0*m)) * W->CE->Jx[n][i][j][k];
		}
	}
	return f;
}

float raweEffy(world W, int i, int j, int k)
{
	if (!W->CE->y) return 1;
	if (i<0 && W->xMinSur==SYM) i=-i;
	if (j<0 && W->yMinSur==SYM) j=1-j;
	if (k<0 && W->zMinSur==SYM) k=-k;
	float f = 1 / W->CE->y[i][j][k];
	for (int n=0, o=0; n<W->CE->N; n++, o++) {
		while (!auxiliaryFieldRequired(W->CE->M[o])) o++;
		if (W->CE->iMin[n]<=i && i<=W->CE->iMax[n] && W->CE->jMin[n]<=j && j<=W->CE->jMax[n] && W->CE->kMin[n]<=k && k<=W->CE->kMax[n]) {
			int m = W->CE->M[o].e.y || W->CE->M[o].e.z ? 1 : 0;
			f -= W->CE->ay[n] * (W->CE->M[o].e.y ? W->CE->M[o].e.y : W->CE->M[o].e.x) * W->CE->Jy[n][i][j][k];
			for (int p=0; p<MAXPOLES; p+=1+2*m)
				f += (pole(W, W->CE->M[o], W->f*W->eV, p+1*m) + W->f*W->eV * dpole(W, W->CE->M[o], W->f*W->eV, p+1*m)) * W->CE->Jy[n][i][j][k];
		}
	}
	return f;
}

float raweEffz(world W, int i, int j, int k)
{
	if (!W->CE->z) return 1;
	if (i<0 && W->xMinSur==SYM) i=-i;
	if (j<0 && W->yMinSur==SYM) j=-j;
	if (k<0 && W->zMinSur==SYM) k=1-k;
	float f = 1 / W->CE->z[i][j][k];
	for (int n=0, o=0; n<W->CE->N; n++, o++) {
		while (!auxiliaryFieldRequired(W->CE->M[o])) o++;
		if (W->CE->iMin[n]<=i && i<=W->CE->iMax[n] && W->CE->jMin[n]<=j && j<=W->CE->jMax[n] && W->CE->kMin[n]<=k && k<=W->CE->kMax[n]) {
			int m = W->CE->M[o].e.y || W->CE->M[o].e.z ? 1 : 0;
			f -= W->CE->az[n] * (W->CE->M[o].e.z ? W->CE->M[o].e.z : W->CE->M[o].e.x) * W->CE->Jz[n][i][j][k];
			for (int p=0; p<MAXPOLES; p+=1+2*m)
				f += (pole(W, W->CE->M[o], W->f*W->eV, p+2*m) + W->f*W->eV * dpole(W, W->CE->M[o], W->f*W->eV, p+2*m)) * W->CE->Jz[n][i][j][k];
		}
	}
	return f;
}

float eEff(world W, int i, int j, int k) {
	if (i<0 && W->xMinSur==SYM) i=-i;
	if (j<0 && W->yMinSur==SYM) j=-j;
	if (k<0 && W->zMinSur==SYM) k=-k;
	return (raweEffx(W,i,j,k) + raweEffx(W,i+1,j,k) + raweEffy(W,i,j,k) + raweEffy(W,i,j+1,k) + raweEffz(W,i,j,k) + raweEffz(W,i,j,k+1)) / 6;
}

float RI(world W, int i, int j, int k) {
	if (i<0 && W->xMinSur==SYM) i=-i;
	if (j<0 && W->yMinSur==SYM) j=-j;
	if (k<0 && W->zMinSur==SYM) k=-k;
	return sqrtf(fabsf(eEff(W,i,j,k)));
}

float ContourRI(world W, int i, int j, int k)
{
	float x=0, y=0, z=0;
	i = W->xtoi(W, i*W->dx), j = W->ytoj(W, j*W->dy), k = W->ztok(W, k*W->dz);
	if (i<0 && W->xMinSur==SYM) i=-i;
	if (j<0 && W->yMinSur==SYM) j=-j;
	if (k<0 && W->zMinSur==SYM) k=-k;
	x = log10f(1+sqrtf(fabsf(raweEffx(W, i+1, j, k)))) - log10f(1+sqrtf(fabsf(raweEffx(W, i, j, k))));
	y = log10f(1+sqrtf(fabsf(raweEffy(W, i, j+1, k)))) - log10f(1+sqrtf(fabsf(raweEffy(W, i, j, k))));
	z = log10f(1+sqrtf(fabsf(raweEffz(W, i, j, k+1)))) - log10f(1+sqrtf(fabsf(raweEffz(W, i, j, k))));
	return sqrt(x*x+y*y+z*z);
}

float LogRI(world W, int i, int j, int k) { return log10f(RI(W, i, j, k)); }
float Contour(world W, int i, int j, int k) { return ContourRI(W, i, j, k) ? 1 : 0; }


/* Electric and magnetic fields */

float rawEx(world W, int i, int j, int k) {
	float p=1;
	if (i<0 && W->xMinSur==SYM) i=-i, p*= W->coskx;
	if (j<0 && W->yMinSur==SYM) j=-j, p*=-W->cosky;
	if (k<0 && W->zMinSur==SYM) k=-k, p*=-W->coskz;
	return p * W->E->x[i][j][k];
}

float rawEy(world W, int i, int j, int k) {
	float p=1;
	if (i<0 && W->xMinSur==SYM) i=-i, p*=-W->coskx;
	if (j<0 && W->yMinSur==SYM) j=-j, p*= W->cosky;
	if (k<0 && W->zMinSur==SYM) k=-k, p*=-W->coskz;
	return p * W->E->y[i][j][k];
}

float rawEz(world W, int i, int j, int k) {
	float p=1;
	if (i<0 && W->xMinSur==SYM) i=-i, p*=-W->coskx;
	if (j<0 && W->yMinSur==SYM) j=-j, p*=-W->cosky;
	if (k<0 && W->zMinSur==SYM) k=-k, p*= W->coskz;
	return p * W->E->z[i][j][k];
}

float rawHx(world W, int i, int j, int k) {
	float p=1;
	if (i<0 && W->xMinSur==SYM) i=-i, p*=-W->coskx;
	if (j<0 && W->yMinSur==SYM) j=-j, p*= W->cosky;
	if (k<0 && W->zMinSur==SYM) k=-k, p*= W->coskz;
	return p * W->H->x[i][j][k];
}

float rawHy(world W, int i, int j, int k) {
	float p=1;
	if (i<0 && W->xMinSur==SYM) i=-i, p*= W->coskx;
	if (j<0 && W->yMinSur==SYM) j=-j, p*=-W->cosky;
	if (k<0 && W->zMinSur==SYM) k=-k, p*= W->coskz;
	return p * W->H->y[i][j][k];
}

float rawHz(world W, int i, int j, int k) {
	float p=1;
	if (i<0 && W->xMinSur==SYM) i=-i, p*= W->coskx;
	if (j<0 && W->yMinSur==SYM) j=-j, p*= W->cosky;
	if (k<0 && W->zMinSur==SYM) k=-k, p*=-W->coskz;
	return p * W->H->z[i][j][k];
}


float Ex(world W, int i, int j, int k)
{
	float p=0.5;
	if (i<0 && W->xMinSur==SYM) i=-i, p*= W->coskx;
	if (j<0 && W->yMinSur==SYM) j=-j, p*=-W->cosky;
	if (k<0 && W->zMinSur==SYM) k=-k, p*=-W->coskz;
	return p * (W->E->x[i][j][k] + W->E->x[i+1][j][k]);
}

float Ey(world W, int i, int j, int k)
{
	float p=0.5;
	if (i<0 && W->xMinSur==SYM) i=-i, p*=-W->coskx;
	if (j<0 && W->yMinSur==SYM) j=-j, p*= W->cosky;
	if (k<0 && W->zMinSur==SYM) k=-k, p*=-W->coskz;
	return p * (W->E->y[i][j][k] + W->E->y[i][j+1][k]);
}

float Ez(world W, int i, int j, int k)
{
	float p=0.5;
	if (i<0 && W->xMinSur==SYM) i=-i, p*=-W->coskx;
	if (j<0 && W->yMinSur==SYM) j=-j, p*=-W->cosky;
	if (k<0 && W->zMinSur==SYM) k=-k, p*= W->coskz;
	return p * (W->E->z[i][j][k] + W->E->z[i][j][k+1]);
}

float iEx(world W, int i, int j, int k)
{
	float p=0.5;
	if (!W->complexField) return 0;
	if (i<0 && W->xMinSur==SYM) i=-i, p*= W->coskx;
	if (j<0 && W->yMinSur==SYM) j=-j, p*=-W->cosky;
	if (k<0 && W->zMinSur==SYM) k=-k, p*=-W->coskz;
	return p * (W->E[1].x[i][j][k] + W->E[1].x[i+1][j][k]);
}

float iEy(world W, int i, int j, int k)
{
	float p=0.5;
	if (!W->complexField) return 0;
	if (i<0 && W->xMinSur==SYM) i=-i, p*=-W->coskx;
	if (j<0 && W->yMinSur==SYM) j=-j, p*= W->cosky;
	if (k<0 && W->zMinSur==SYM) k=-k, p*=-W->coskz;
	return p * (W->E[1].y[i][j][k] + W->E[1].y[i][j+1][k]);
}

float iEz(world W, int i, int j, int k)
{
	float p=0.5;
	if (!W->complexField) return 0;
	if (i<0 && W->xMinSur==SYM) i=-i, p*=-W->coskx;
	if (j<0 && W->yMinSur==SYM) j=-j, p*=-W->cosky;
	if (k<0 && W->zMinSur==SYM) k=-k, p*= W->coskz;
	return p * (W->E[1].z[i][j][k] + W->E[1].z[i][j][k+1]);
}

float Hx(world W, int i, int j, int k)
{
	float p=0.25;
	if (i<0 && W->xMinSur==SYM) i=-i, p*=-W->coskx;
	if (j<0 && W->yMinSur==SYM) j=-j, p*= W->cosky;
	if (k<0 && W->zMinSur==SYM) k=-k, p*= W->coskz;
	return p * (W->H->x[i][j][k] + W->H->x[i][j+1][k] + W->H->x[i][j][k+1] + W->H->x[i][j+1][k+1]);
}

float Hy(world W, int i, int j, int k)
{
	float p=0.25;
	if (i<0 && W->xMinSur==SYM) i=-i, p*= W->coskx;
	if (j<0 && W->yMinSur==SYM) j=-j, p*=-W->cosky;
	if (k<0 && W->zMinSur==SYM) k=-k, p*= W->coskz;
	return p * (W->H->y[i][j][k] + W->H->y[i][j][k+1] + W->H->y[i+1][j][k] + W->H->y[i+1][j][k+1]);
}

float Hz(world W, int i, int j, int k)
{
	float p=0.25;
	if (i<0 && W->xMinSur==SYM) i=-i, p*= W->coskx;
	if (j<0 && W->yMinSur==SYM) j=-j, p*= W->cosky;
	if (k<0 && W->zMinSur==SYM) k=-k, p*=-W->coskz;
	return p * (W->H->z[i][j][k] + W->H->z[i+1][j][k] + W->H->z[i][j+1][k] + W->H->z[i+1][j+1][k]);
}

float iHx(world W, int i, int j, int k)
{
	float p=0.25;
	if (!W->complexField) return 0;
	if (i<0 && W->xMinSur==SYM) i=-i, p*=-W->coskx;
	if (j<0 && W->yMinSur==SYM) j=-j, p*= W->cosky;
	if (k<0 && W->zMinSur==SYM) k=-k, p*= W->coskz;
	return p * (W->H[1].x[i][j][k] + W->H[1].x[i][j+1][k] + W->H[1].x[i][j][k+1] + W->H[1].x[i][j+1][k+1]);
}

float iHy(world W, int i, int j, int k)
{
	float p=0.25;
	if (!W->complexField) return 0;
	if (i<0 && W->xMinSur==SYM) i=-i, p*= W->coskx;
	if (j<0 && W->yMinSur==SYM) j=-j, p*=-W->cosky;
	if (k<0 && W->zMinSur==SYM) k=-k, p*= W->coskz;
	return p * (W->H[1].y[i][j][k] + W->H[1].y[i][j][k+1] + W->H[1].y[i+1][j][k] + W->H[1].y[i+1][j][k+1]);
}

float iHz(world W, int i, int j, int k)
{
	float p=0.25;
	if (!W->complexField) return 0;
	if (i<0 && W->xMinSur==SYM) i=-i, p*= W->coskx;
	if (j<0 && W->yMinSur==SYM) j=-j, p*= W->cosky;
	if (k<0 && W->zMinSur==SYM) k=-k, p*=-W->coskz;
	return p * (W->H[1].z[i][j][k] + W->H[1].z[i+1][j][k] + W->H[1].z[i][j+1][k] + W->H[1].z[i+1][j+1][k]);
}


/* Current density */

float rawJx(world W, int i, int j, int k)
{
	float f = 0, p = 1;
	source S = W->SE;
	if (i<0 && W->xMinSur==SYM) i=1-i, p*= W->coskx;
	if (j<0 && W->yMinSur==SYM) j= -j, p*=-W->cosky;
	if (k<0 && W->zMinSur==SYM) k= -k, p*=-W->coskz;

	if (W->CE->x) {
		f += (W->CE->dtdz[k] * (W->H->y[i][j][k] - W->H->y[i][j][k+1])
			+W->CE->dtdy[j] * (W->H->z[i][j+1][k] - W->H->z[i][j][k])) * (1 - W->CE->x[i][j][k]);
		for (int n=0; n<W->CE->N; n++) {
			if (W->CE->iMin[n]<=i && i<=W->CE->iMax[n] && W->CE->jMin[n]<=j && j<W->CE->jMax[n] && W->CE->kMin[n]<=k && k<W->CE->kMax[n]) {
				if (W->CE->ax[n]) f += 2 * W->CE->Jx[n][i][j][k] * W->E->x[i][j][k] / (1 + 1 / W->CE->ax[n]);
				for (int p=0; p<W->CE->NK[n]; p++) if (W->E->Jx[n][p]) f += W->E->Jx[n][p][i][j][k];
			}
		}
	}
	while (S) {
		if (S->F == W->E->x && W->n < S->N && S->iMin <= i && i <= S->iMax && S->jMin <= j && j <= S->jMax && S->kMin <= k && k <= S->kMax)
			f -= (S->Amplitude ? S->Amplitude[i][j][k] : S->amplitude) * S->waveForm(W, S, W->H->t, S->Phase ? S->Phase[i][j][k] : S->phase);
		S = S->S;
	}
	return p * f / W->dt;
}

float rawJy(world W, int i, int j, int k)
{
	float f = 0, p = 1;
	source S = W->SE;
	if (i<0 && W->xMinSur==SYM) i= -i, p*=-W->coskx;
	if (j<0 && W->yMinSur==SYM) j=1-j, p*= W->cosky;
	if (k<0 && W->zMinSur==SYM) k= -k, p*=-W->coskz;

	if (W->CE->y) {
		f += (W->CE->dtdx[i] * (W->H->z[i][j][k] - W->H->z[i+1][j][k])
			+W->CE->dtdz[k] * (W->H->x[i][j][k+1] - W->H->x[i][j][k])) * (1 - W->CE->y[i][j][k]);
		for (int n=0; n<W->CE->N; n++) {
			if (W->CE->iMin[n]<=i && i<W->CE->iMax[n] && W->CE->jMin[n]<=j && j<=W->CE->jMax[n] && W->CE->kMin[n]<=k && k<W->CE->kMax[n]) {
				if (W->CE->ay[n]) f += 2 * W->CE->Jy[n][i][j][k] * W->E->y[i][j][k] / (1 + 1 / W->CE->ay[n]);
				for (int p=0; p<W->CE->NK[n]; p++) if (W->E->Jy[n][p]) f += W->E->Jy[n][p][i][j][k];
			}
		}
	}
	while (S) {
		if (S->F == W->E->y && W->n < S->N && S->iMin <= i && i <= S->iMax && S->jMin <= j && j <= S->jMax && S->kMin <= k && k <= S->kMax)
			f -= (S->Amplitude ? S->Amplitude[i][j][k] : S->amplitude) * S->waveForm(W, S, W->H->t, S->Phase ? S->Phase[i][j][k] : S->phase);
		S = S->S;
	}
	return p * f / W->dt;
}

float rawJz(world W, int i, int j, int k)
{
	float f = 0, p = 1;
	source S = W->SE;
	if (i<0 && W->xMinSur==SYM) i= -i, p*=-W->coskx;
	if (j<0 && W->yMinSur==SYM) j= -j, p*=-W->cosky;
	if (k<0 && W->zMinSur==SYM) k=1-k, p*= W->coskz;

	if (W->CE->z) {
		f += (W->CE->dtdy[j] * (W->H->x[i][j][k] - W->H->x[i][j+1][k])
			+W->CE->dtdx[i] * (W->H->y[i+1][j][k] - W->H->y[i][j][k])) * (1 - W->CE->z[i][j][k]);
		for (int n=0; n<W->CE->N; n++) {
			if (W->CE->iMin[n]<=i && i<W->CE->iMax[n] && W->CE->jMin[n]<=j && j<W->CE->jMax[n] && W->CE->kMin[n]<=k && k<=W->CE->kMax[n]) {
				if (W->CE->az[n]) f += 2 * W->CE->Jz[n][i][j][k] * W->E->z[i][j][k] / (1 + 1 / W->CE->az[n]);
				for (int p=0; p<W->CE->NK[n]; p++) if (W->E->Jz[n][p]) f += W->E->Jz[n][p][i][j][k];
			}
		}
	}
	while (S) {
		if (S->F == W->E->z && W->n < S->N && S->iMin <= i && i <= S->iMax && S->jMin <= j && j <= S->jMax && S->kMin <= k && k <= S->kMax)
			f -= (S->Amplitude ? S->Amplitude[i][j][k] : S->amplitude) * S->waveForm(W, S, W->H->t, S->Phase ? S->Phase[i][j][k] : S->phase);
		S = S->S;
	}
	return p * f / W->dt;
}

float Jx(world W, int i, int j, int k) { return 0.5 * (rawJx(W,i,j,k) + rawJx(W,i+1,j,k)); }
float Jy(world W, int i, int j, int k) { return 0.5 * (rawJy(W,i,j,k) + rawJy(W,i,j+1,k)); }
float Jz(world W, int i, int j, int k) { return 0.5 * (rawJz(W,i,j,k) + rawJz(W,i,j,k+1)); }


/* Fields related to energy */

float ExHy(world W, int i, int j, int k)
{
	float p=0.25;
	if (i<0 && W->xMinSur==SYM) i=-i, p*= 1;
	if (j<0 && W->yMinSur==SYM) j=-j, p*= 1;
	if (k<0 && W->zMinSur==SYM) k=-k, p*=-1;
	return p * (W->E->x[i][j][k] * (W->H->y[i][j][k] + W->H->y[i][j][k+1]) + W->E->x[i+1][j][k] * (W->H->y[i+1][j][k] + W->H->y[i+1][j][k+1]));
}

float ExHz(world W, int i, int j, int k)
{
	float p=0.25;
	if (i<0 && W->xMinSur==SYM) i=-i, p*= 1;
	if (j<0 && W->yMinSur==SYM) j=-j, p*=-1;
	if (k<0 && W->zMinSur==SYM) k=-k, p*= 1;
	return p * (W->E->x[i][j][k] * (W->H->z[i][j][k] + W->H->z[i][j+1][k]) + W->E->x[i+1][j][k] * (W->H->z[i+1][j][k] + W->H->z[i+1][j+1][k]));
}

float EyHz(world W, int i, int j, int k)
{
	float p=0.25;
	if (i<0 && W->xMinSur==SYM) i=-i, p*=-1;
	if (j<0 && W->yMinSur==SYM) j=-j, p*= 1;
	if (k<0 && W->zMinSur==SYM) k=-k, p*= 1;
	return p * (W->E->y[i][j][k] * (W->H->z[i][j][k] + W->H->z[i+1][j][k]) + W->E->y[i][j+1][k] * (W->H->z[i][j+1][k] + W->H->z[i+1][j+1][k]));
}

float EyHx(world W, int i, int j, int k)
{
	float p=0.25;
	if (i<0 && W->xMinSur==SYM) i=-i, p*= 1;
	if (j<0 && W->yMinSur==SYM) j=-j, p*= 1;
	if (k<0 && W->zMinSur==SYM) k=-k, p*=-1;
	return p * (W->E->y[i][j][k] * (W->H->x[i][j][k] + W->H->x[i][j][k+1]) + W->E->y[i][j+1][k] * (W->H->x[i][j+1][k] + W->H->x[i][j+1][k+1]));
}

float EzHx(world W, int i, int j, int k)
{
	float p=0.25;
	if (i<0 && W->xMinSur==SYM) i=-i, p*= 1;
	if (j<0 && W->yMinSur==SYM) j=-j, p*=-1;
	if (k<0 && W->zMinSur==SYM) k=-k, p*= 1;
	return p * (W->E->z[i][j][k] * (W->H->x[i][j][k] + W->H->x[i][j+1][k]) + W->E->z[i][j][k+1] * (W->H->x[i][j][k+1] + W->H->x[i][j+1][k+1]));
}

float EzHy(world W, int i, int j, int k)
{
	float p=0.25;
	if (i<0 && W->xMinSur==SYM) i=-i, p*=-1;
	if (j<0 && W->yMinSur==SYM) j=-j, p*= 1;
	if (k<0 && W->zMinSur==SYM) k=-k, p*= 1;
	return p * (W->E->z[i][j][k] * (W->H->y[i][j][k] + W->H->y[i+1][j][k]) + W->E->z[i][j][k+1] * (W->H->y[i][j][k+1] + W->H->y[i+1][j][k+1]));
}

float JE(world W, int i, int j, int k)
{
	float f=0;
	if (i<0 && W->xMinSur==SYM) i=-i;
	if (j<0 && W->yMinSur==SYM) j=-j;
	if (k<0 && W->zMinSur==SYM) k=-k;
	f += W->E->x[i][j][k] * rawJx(W,i,j,k) + W->E->x[i+1][j][k] * rawJx(W,i+1,j,k);
	f += W->E->y[i][j][k] * rawJy(W,i,j,k) + W->E->y[i][j+1][k] * rawJy(W,i,j+1,k);
	f += W->E->z[i][j][k] * rawJz(W,i,j,k) + W->E->z[i][j][k+1] * rawJz(W,i,j,k+1);
	return 0.5 * f;
}

float EE(world W, int i, int j, int k)
{
	float f=0;
	if (i<0 && W->xMinSur==SYM) i=-i;
	if (j<0 && W->yMinSur==SYM) j=-j;
	if (k<0 && W->zMinSur==SYM) k=-k;
	f += sq(W->E->x[i][j][k]) + sq(W->E->x[i+1][j][k]);
	f += sq(W->E->y[i][j][k]) + sq(W->E->y[i][j+1][k]);
	f += sq(W->E->z[i][j][k]) + sq(W->E->z[i][j][k+1]);
	return 0.5 * f;
}

float UE(world W, int i, int j, int k)
{
	float f=0;
	if (i<0 && W->xMinSur==SYM) i=-i;
	if (j<0 && W->yMinSur==SYM) j=-j;
	if (k<0 && W->zMinSur==SYM) k=-k;
	f += sq(W->E->x[i][j][k]) * raweEffx(W, i, j, k) + sq(W->E->x[i+1][j][k]) * raweEffx(W, i+1, j, k);
	f += sq(W->E->y[i][j][k]) * raweEffy(W, i, j, k) + sq(W->E->y[i][j+1][k]) * raweEffy(W, i, j+1, k);
	f += sq(W->E->z[i][j][k]) * raweEffz(W, i, j, k) + sq(W->E->z[i][j][k+1]) * raweEffz(W, i, j, k+1);
	return 0.25 * f;
}

float HH(world W, int i, int j, int k)
{
	float f=0;
	if (i<0 && W->xMinSur==SYM) i=-i;
	if (j<0 && W->yMinSur==SYM) j=-j;
	if (k<0 && W->zMinSur==SYM) k=-k;
	f += sq(W->H->x[i][j][k]) + sq(W->H->x[i][j+1][k]) + sq(W->H->x[i][j][k+1]) + sq(W->H->x[i][j+1][k+1]);
	f += sq(W->H->y[i][j][k]) + sq(W->H->y[i][j][k+1]) + sq(W->H->y[i+1][j][k]) + sq(W->H->y[i+1][j][k+1]);
	f += sq(W->H->z[i][j][k]) + sq(W->H->z[i+1][j][k]) + sq(W->H->z[i][j+1][k]) + sq(W->H->z[i+1][j+1][k]);
	return 0.25 * f;
}

float Sx(world W, int i, int j, int k) { return EyHz(W,i,j,k) - EzHy(W,i,j,k); }
float Sy(world W, int i, int j, int k) { return EzHx(W,i,j,k) - ExHz(W,i,j,k); }
float Sz(world W, int i, int j, int k) { return ExHy(W,i,j,k) - EyHx(W,i,j,k); }
float UH(world W, int i, int j, int k) { return 0.5*HH(W, i, j, k); }
float U(world W, int i, int j, int k) { return UE(W, i, j, k)+UH(W, i, j, k); }
float LogJE(world W, int i, int j, int k) { return log10f(JE(W, i, j, k)); }
float LogEE(world W, int i, int j, int k) { return log10f(EE(W, i, j, k)); }
float LogUE(world W, int i, int j, int k) { return log10f(UE(W, i, j, k)); }
float LogHH(world W, int i, int j, int k) { return log10f(HH(W, i, j, k)); }
float LogUH(world W, int i, int j, int k) { return log10f(HH(W, i, j, k)); }
float LogU(world W, int i, int j, int k) { return log10f(U(W, i, j, k)); }


/* Fields related to optical force */

float JxHy(world W, int i, int j, int k)
{
	float p=0.25;
	if (i<0 && W->xMinSur==SYM) i=-i, p*= 1;
	if (j<0 && W->yMinSur==SYM) j=-j, p*= 1;
	if (k<0 && W->zMinSur==SYM) k=-k, p*=-1;
	return p * (rawJx(W,i,j,k) * (W->H->y[i][j][k] + W->H->y[i][j][k+1]) + rawJx(W,i+1,j,k) * (W->H->y[i+1][j][k] + W->H->y[i+1][j][k+1]));
}

float JxHz(world W, int i, int j, int k)
{
	float p=0.25;
	if (i<0 && W->xMinSur==SYM) i=-i, p*= 1;
	if (j<0 && W->yMinSur==SYM) j=-j, p*=-1;
	if (k<0 && W->zMinSur==SYM) k=-k, p*= 1;
	return p * (rawJx(W,i,j,k) * (W->H->z[i][j][k] + W->H->z[i][j+1][k]) + rawJx(W,i+1,j,k) * (W->H->z[i+1][j][k] + W->H->z[i+1][j+1][k]));
}

float JyHz(world W, int i, int j, int k)
{
	float p=0.25;
	if (i<0 && W->xMinSur==SYM) i=-i, p*=-1;
	if (j<0 && W->yMinSur==SYM) j=-j, p*= 1;
	if (k<0 && W->zMinSur==SYM) k=-k, p*= 1;
	return p * (rawJy(W,i,j,k) * (W->H->z[i][j][k] + W->H->z[i+1][j][k]) + rawJy(W,i,j+1,k) * (W->H->z[i][j+1][k] + W->H->z[i+1][j+1][k]));
}

float JyHx(world W, int i, int j, int k)
{
	float p=0.25;
	if (i<0 && W->xMinSur==SYM) i=-i, p*= 1;
	if (j<0 && W->yMinSur==SYM) j=-j, p*= 1;
	if (k<0 && W->zMinSur==SYM) k=-k, p*=-1;
	return p * (rawJy(W,i,j,k) * (W->H->x[i][j][k] + W->H->x[i][j][k+1]) + rawJy(W,i,j+1,k) * (W->H->x[i][j+1][k] + W->H->x[i][j+1][k+1]));
}

float JzHx(world W, int i, int j, int k)
{
	float p=0.25;
	if (i<0 && W->xMinSur==SYM) i=-i, p*= 1;
	if (j<0 && W->yMinSur==SYM) j=-j, p*=-1;
	if (k<0 && W->zMinSur==SYM) k=-k, p*= 1;
	return p * (rawJz(W,i,j,k) * (W->H->x[i][j][k] + W->H->x[i][j+1][k]) + rawJz(W,i,j,k+1) * (W->H->x[i][j][k+1] + W->H->x[i][j+1][k+1]));
}

float JzHy(world W, int i, int j, int k)
{
	float p=0.25;
	if (i<0 && W->xMinSur==SYM) i=-i, p*=-1;
	if (j<0 && W->yMinSur==SYM) j=-j, p*= 1;
	if (k<0 && W->zMinSur==SYM) k=-k, p*= 1;
	return p * (rawJz(W,i,j,k) * (W->H->y[i][j][k] + W->H->y[i+1][j][k]) + rawJz(W,i,j,k+1) * (W->H->y[i][j][k+1] + W->H->y[i+1][j][k+1]));
}

float DivE(world W, int i, int j, int k)
{
	float p=1;
	if (i<0 && W->xMinSur==SYM) i=-i, p*=-W->coskx;
	if (j<0 && W->yMinSur==SYM) j=-j, p*=-W->cosky;
	if (k<0 && W->zMinSur==SYM) k=-k, p*=-W->coskz;
	return ((W->E->x[i+1][j][k] - W->E->x[i][j][k]) * W->CE->dtdx[i] +
		   	(W->E->y[i][j+1][k] - W->E->y[i][j][k]) * W->CE->dtdy[j] +
		   	(W->E->z[i][j][k+1] - W->E->z[i][j][k]) * W->CE->dtdz[k]) * p / W->dt;
}

float Txx(world W, int i, int j, int k)
{
	float f=0;
	if (i<0 && W->xMinSur==SYM) i=-i;
	if (j<0 && W->yMinSur==SYM) j=-j;
	if (k<0 && W->zMinSur==SYM) k=-k;
	f += sq(W->H->x[i][j][k]) + sq(W->H->x[i][j+1][k]) + sq(W->H->x[i][j][k+1]) + sq(W->H->x[i][j+1][k+1]);
	f -= sq(W->H->y[i][j][k]) + sq(W->H->y[i][j][k+1]) + sq(W->H->y[i+1][j][k]) + sq(W->H->y[i+1][j][k+1]);
	f -= sq(W->H->z[i][j][k]) + sq(W->H->z[i+1][j][k]) + sq(W->H->z[i][j+1][k]) + sq(W->H->z[i+1][j+1][k]);
	f *= 0.5;
	f += sq(W->E->x[i][j][k]) + sq(W->E->x[i+1][j][k]);
	f -= sq(W->E->y[i][j][k]) + sq(W->E->y[i][j+1][k]);
	f -= sq(W->E->z[i][j][k]) + sq(W->E->z[i][j][k+1]);
	return 0.25 * f;
}

float Tyy(world W, int i, int j, int k)
{
	float f=0;
	if (i<0 && W->xMinSur==SYM) i=-i;
	if (j<0 && W->yMinSur==SYM) j=-j;
	if (k<0 && W->zMinSur==SYM) k=-k;
	f -= sq(W->H->x[i][j][k]) + sq(W->H->x[i][j+1][k]) + sq(W->H->x[i][j][k+1]) + sq(W->H->x[i][j+1][k+1]);
	f += sq(W->H->y[i][j][k]) + sq(W->H->y[i][j][k+1]) + sq(W->H->y[i+1][j][k]) + sq(W->H->y[i+1][j][k+1]);
	f -= sq(W->H->z[i][j][k]) + sq(W->H->z[i+1][j][k]) + sq(W->H->z[i][j+1][k]) + sq(W->H->z[i+1][j+1][k]);
	f *= 0.5;
	f -= sq(W->E->x[i][j][k]) + sq(W->E->x[i+1][j][k]);
	f += sq(W->E->y[i][j][k]) + sq(W->E->y[i][j+1][k]);
	f -= sq(W->E->z[i][j][k]) + sq(W->E->z[i][j][k+1]);
	return 0.25 * f;
}

float Tzz(world W, int i, int j, int k)
{
	float f=0;
	if (i<0 && W->xMinSur==SYM) i=-i;
	if (j<0 && W->yMinSur==SYM) j=-j;
	if (k<0 && W->zMinSur==SYM) k=-k;
	f -= sq(W->H->x[i][j][k]) + sq(W->H->x[i][j+1][k]) + sq(W->H->x[i][j][k+1]) + sq(W->H->x[i][j+1][k+1]);
	f -= sq(W->H->y[i][j][k]) + sq(W->H->y[i][j][k+1]) + sq(W->H->y[i+1][j][k]) + sq(W->H->y[i+1][j][k+1]);
	f += sq(W->H->z[i][j][k]) + sq(W->H->z[i+1][j][k]) + sq(W->H->z[i][j+1][k]) + sq(W->H->z[i+1][j+1][k]);
	f *= 0.5;
	f -= sq(W->E->x[i][j][k]) + sq(W->E->x[i+1][j][k]);
	f -= sq(W->E->y[i][j][k]) + sq(W->E->y[i][j+1][k]);
	f += sq(W->E->z[i][j][k]) + sq(W->E->z[i][j][k+1]);
	return 0.25 * f;
}

float Txy(world W, int i, int j, int k)
{
	float f=0, p=0.25;
	if (i<0 && W->xMinSur==SYM) i=-i, p*=-1;
	if (j<0 && W->yMinSur==SYM) j=-j, p*=-1;
	if (k<0 && W->zMinSur==SYM) k=-k, p*= 1;
	f += (W->E->x[i][j][k] + W->E->x[i+1][j][k]) * (W->E->y[i][j][k] + W->E->y[i][j+1][k]);
	f += (W->H->x[i][j][k] + W->H->x[i][j+1][k] + W->H->x[i][j][k+1] + W->H->x[i][j+1][k+1])
	     * (W->H->y[i][j][k] + W->H->y[i][j][k+1] + W->H->y[i+1][j][k] + W->H->y[i+1][j][k+1]) * 0.25;
	return p * f;
}

float Txz(world W, int i, int j, int k)
{
	float f=0, p=0.25;
	if (i<0 && W->xMinSur==SYM) i=-i, p*=-1;
	if (j<0 && W->yMinSur==SYM) j=-j, p*= 1;
	if (k<0 && W->zMinSur==SYM) k=-k, p*=-1;
	f += (W->E->z[i][j][k] + W->E->z[i][j][k+1]) * (W->E->x[i][j][k] + W->E->x[i+1][j][k]);
	f += (W->H->z[i][j][k] + W->H->z[i+1][j][k] + W->H->z[i][j+1][k] + W->H->z[i+1][j+1][k])
	     * (W->H->x[i][j][k] + W->H->x[i][j+1][k] + W->H->x[i][j][k+1] + W->H->x[i][j+1][k+1]) * 0.25;
	return p * f;
}

float Tyx(world W, int i, int j, int k)
{
	float f=0, p=0.25;
	if (i<0 && W->xMinSur==SYM) i=-i, p*=-1;
	if (j<0 && W->yMinSur==SYM) j=-j, p*=-1;
	if (k<0 && W->zMinSur==SYM) k=-k, p*= 1;
	f += (W->E->x[i][j][k] + W->E->x[i+1][j][k]) * (W->E->y[i][j][k] + W->E->y[i][j+1][k]);
	f += (W->H->x[i][j][k] + W->H->x[i][j+1][k] + W->H->x[i][j][k+1] + W->H->x[i][j+1][k+1])
	     * (W->H->y[i][j][k] + W->H->y[i][j][k+1] + W->H->y[i+1][j][k] + W->H->y[i+1][j][k+1]) * 0.25;
	return p * f;
}

float Tyz(world W, int i, int j, int k)
{
	float f=0, p=0.25;
	if (i<0 && W->xMinSur==SYM) i=-i, p*= 1;
	if (j<0 && W->yMinSur==SYM) j=-j, p*=-1;
	if (k<0 && W->zMinSur==SYM) k=-k, p*=-1;
	f += (W->E->y[i][j][k] + W->E->y[i][j+1][k]) * (W->E->z[i][j][k] + W->E->z[i][j][k+1]);
	f += (W->H->y[i][j][k] + W->H->y[i][j][k+1] + W->H->y[i+1][j][k] + W->H->y[i+1][j][k+1])
	     * (W->H->z[i][j][k] + W->H->z[i+1][j][k] + W->H->z[i][j+1][k] + W->H->z[i+1][j+1][k]) * 0.25;
	return p * f;
}

float Tzx(world W, int i, int j, int k)
{
	float f=0, p=0.25;
	if (i<0 && W->xMinSur==SYM) i=-i, p*=-1;
	if (j<0 && W->yMinSur==SYM) j=-j, p*= 1;
	if (k<0 && W->zMinSur==SYM) k=-k, p*=-1;
	f += (W->E->z[i][j][k] + W->E->z[i][j][k+1]) * (W->E->x[i][j][k] + W->E->x[i+1][j][k]);
	f += (W->H->z[i][j][k] + W->H->z[i+1][j][k] + W->H->z[i][j+1][k] + W->H->z[i+1][j+1][k])
	     * (W->H->x[i][j][k] + W->H->x[i][j+1][k] + W->H->x[i][j][k+1] + W->H->x[i][j+1][k+1]) * 0.25;
	return p * f;
}

float Tzy(world W, int i, int j, int k)
{
	float f=0, p=0.25;
	if (i<0 && W->xMinSur==SYM) i=-i, p*= 1;
	if (j<0 && W->yMinSur==SYM) j=-j, p*=-1;
	if (k<0 && W->zMinSur==SYM) k=-k, p*=-1;
	f += (W->E->y[i][j][k] + W->E->y[i][j+1][k]) * (W->E->z[i][j][k] + W->E->z[i][j][k+1]);
	f += (W->H->y[i][j][k] + W->H->y[i][j][k+1] + W->H->y[i+1][j][k] + W->H->y[i+1][j][k+1])
	     * (W->H->z[i][j][k] + W->H->z[i+1][j][k] + W->H->z[i][j+1][k] + W->H->z[i+1][j+1][k]) * 0.25;
	return p * f;
}

float Fx(world W, int i, int j, int k) { return DivE(W, i, j, k) * Ex(W, i, j, k) + (JyHz(W, i, j, k) - JzHy(W, i, j, k)); }
float Fy(world W, int i, int j, int k) { return DivE(W, i, j, k) * Ey(W, i, j, k) + (JzHx(W, i, j, k) - JxHz(W, i, j, k)); }
float Fz(world W, int i, int j, int k) { return DivE(W, i, j, k) * Ez(W, i, j, k) + (JxHy(W, i, j, k) - JyHx(W, i, j, k)); }


/* Source electric fields */

float SoEx(world W, int i, int j, int k)
{
	tfsf S = W->S;
	float f=0, p=0.5;
	if (i<0 && W->xMinSur==SYM) i=-i, p*= W->coskx;
	if (j<0 && W->yMinSur==SYM) j=-j, p*=-W->cosky;
	if (k<0 && W->zMinSur==SYM) k=-k, p*=-W->coskz;
	while (S) {
		f += S->W->E->x[i][j][k] + S->W->E->x[i+1][j][k];
		S = S->S;
	}
	return p * f;
}

float SoEy(world W, int i, int j, int k)
{
	tfsf S = W->S;
	float f=0, p=0.5;
	if (i<0 && W->xMinSur==SYM) i=-i, p*=-W->coskx;
	if (j<0 && W->yMinSur==SYM) j=-j, p*= W->cosky;
	if (k<0 && W->zMinSur==SYM) k=-k, p*=-W->coskz;
	while (S) {
		f += S->W->E->y[i][j][k] + S->W->E->y[i][j+1][k];
		S = S->S;
	}
	return p * f;
}

float SoEz(world W, int i, int j, int k)
{
	tfsf S = W->S;
	float f=0, p=0.5;
	if (i<0 && W->xMinSur==SYM) i=-i, p*=-W->coskx;
	if (j<0 && W->yMinSur==SYM) j=-j, p*=-W->cosky;
	if (k<0 && W->zMinSur==SYM) k=-k, p*= W->coskz;
	while (S) {
		f += S->W->E->z[i][j][k] + S->W->E->z[i][j][k+1];
		S = S->S;
	}
	return p * f;
}

float SoEE(world W, int i, int j, int k) { return sq(SoEx(W,i,j,k))+sq(SoEy(W,i,j,k))+sq(SoEz(W,i,j,k)); }


/* Scattered electric fields */

float ScEx(world W, int i, int j, int k)
{
	tfsf S = W->S;
	float f, p=0.5;
	if (i<0 && W->xMinSur==SYM) i=-i, p*= W->coskx;
	if (j<0 && W->yMinSur==SYM) j=-j, p*=-W->cosky;
	if (k<0 && W->zMinSur==SYM) k=-k, p*=-W->coskz;
	f = W->E->x[i][j][k] + W->E->x[i+1][j][k];
	while (S) {
		if (S->jMin <= j && j <= S->jMax && S->kMin <= k && k <= S->kMax) {
			if (S->iMIN <= i && i <= S->iMAX) f -= S->W->E->x[i][j][k];
			if (S->iMIN <= i+1 && i+1 <= S->iMAX) f -= S->W->E->x[i+1][j][k];
		}
		S = S->S;
	}
	return p * f;
}

float ScEy(world W, int i, int j, int k)
{
	tfsf S = W->S;
	float f, p=0.5;
	if (i<0 && W->xMinSur==SYM) i=-i, p*=-W->coskx;
	if (j<0 && W->yMinSur==SYM) j=-j, p*= W->cosky;
	if (k<0 && W->zMinSur==SYM) k=-k, p*=-W->coskz;
	f = W->E->y[i][j][k] + W->E->y[i][j+1][k];
	while (S) {
		if (S->iMin <= i && i <= S->iMax && S->kMin <= k && k <= S->kMax) {
			if (S->jMIN <= j && j <= S->jMAX) f -= S->W->E->y[i][j][k];
			if (S->jMIN <= j+1 && j+1 <= S->jMAX) f -= S->W->E->y[i][j+1][k];
		}
		S = S->S;
	}
	return p * f;
}

float ScEz(world W, int i, int j, int k)
{
	tfsf S = W->S;
	float f, p=0.5;
	if (i<0 && W->xMinSur==SYM) i=-i, p*=-W->coskx;
	if (j<0 && W->yMinSur==SYM) j=-j, p*=-W->cosky;
	if (k<0 && W->zMinSur==SYM) k=-k, p*= W->coskz;
	f = W->E->z[i][j][k] + W->E->z[i][j][k+1];
	while (S) {
		if (S->iMin <= i && i <= S->iMax && S->jMin <= j && j <= S->jMax) {
			if (S->kMIN <= k && k <= S->kMAX) f -= S->W->E->z[i][j][k];
			if (S->kMIN <= k+1 && k+1 <= S->kMAX) f -= S->W->E->z[i][j][k+1];
		}
		S = S->S;
	}
	return p * f;
}

float ScEE(world W, int i, int j, int k) { return sq(ScEx(W,i,j,k))+sq(ScEy(W,i,j,k))+sq(ScEz(W,i,j,k)); }
