/**********************************************************************
 * ALiS is Copyright (C) 2009-2019 by Ho-Seok Ee <hsee@kongju.ac.kr>. *
 * Redistribution and use with or without modification, are permitted *
 * under the terms of the Artistic License version 2.                 *
 **********************************************************************/

#include <string.h>
#include "alis_c.h"


float Sine(world W, source S, float t, float phase)
{
	if (!S->N) {
		S->frequency = 2 * PI / S->wavelength;
		S->peaktime = W->dt * floorf(S->wavelength * S->supplement / W->dt);
		S->deviation = S->peaktime / PI;
		S->N = -1;
	}
	if (t < S->peaktime)
		return sinf(S->frequency * (t-S->peaktime) - phase) * (0.5 * (1 - cosf(t / S->deviation)));
	else
		return sinf(S->frequency * (t-S->peaktime) - phase);
}


float Pulse(world W, source S, float t, float phase)
{
	if (!S->N) {
		S->frequency = 2 * PI / S->wavelength;
		S->deviation = 2 * sq(S->wavelength) / (PI * S->supplement);
		S->peaktime = W->dt * floorf(PI * S->deviation / W->dt);
		S->N = 2 * floorf(PI * S->deviation / W->dt);
	}
	return sinf(S->frequency * (t-S->peaktime) - phase) * expf(-sq((t-S->peaktime) / S->deviation));
}


float Band(world W, source S, float t, float phase)
{
	if (!S->N) {
		S->resolution = S->supplement > S->wavelength ? ceilf(1.25 * S->supplement / S->wavelength + 0.5) : 2;
		S->frequency = 2 * PI / S->wavelength;
		S->deviation = 2 * PI / S->wavelength / (S->resolution - 0.5);
		S->peaktime = floorf(S->wavelength * (S->resolution - 0.5) / W->dt) * W->dt;
		S->N = 2 * S->wavelength * (S->resolution - 0.5) / W->dt;
	}
	if (t == S->peaktime) return 0;
	return 2 * (cosf((t-S->peaktime) * S->deviation * S->resolution) - cosf((t-S->peaktime) * S->deviation)) / ((t-S->peaktime) * S->deviation);
}


static source createSource(world W, field F, waveform waveForm, float wavelength, float supplement, float resolution)
{
	source S;
	MALLOC(S, 0, 1);
	*S = (struct source) { waveForm, wavelength, supplement, resolution };
	waveForm(W, S, 0, 0);

	if (F==Ex || F==Ey || F==Ez || F==iEx || F==iEy || F==iEz) {
		if (W->SE) { source s = W->SE; while (s->S) s = s->S; s->S = S; }
		else W->SE = S;
		S->F = F==Ex ? W->E->x : F==Ey ? W->E->y : F==Ez ? W->E->z : F==iEx ? W->E[1].x : F==iEy ? W->E[1].y : F==iEz ? W->E[1].z : 0;
	}

	if (F==Hx || F==Hy || F==Hz || F==iHx || F==iHy || F==iHz) {
		if (W->SH) { source s = W->SH; while (s->S) s = s->S; s->S = S; }
		else W->SH = S;
		S->F = F==Hx ? W->H->x : F==Hy ? W->H->y : F==Hz ? W->H->z : F==iHx ? W->H[1].x : F==iHy ? W->H[1].y : F==iHz ? W->H[1].z : 0;
	}

	if (W->N < S->N && S->N != -1) W->N = S->N;
	if (wavelength && W->T != wavelength/W->dt) W->f = 1/wavelength, W->T = wavelength/W->dt;
	if (!W->sourceType) W->sourceType = 1;
	return S;
}


static void setAmplitude(source S, float *dtdx, float *dtdy, float *dtdz, float ***C, float amplitude)
{
	#pragma omp parallel for
	for (int i=S->iMin; i<=S->iMax; i++) for (int j=S->jMin; j<=S->jMax; j++) for (int k=S->kMin; k<=S->kMax; k++) {
		S->Amplitude[i][j][k] = amplitude * dtdx[i] * dtdy[j] * dtdz[k] * (C ? C[i][j][k] : 1);
	}
}


void (pointDipole)(world W, field F, float x, float y, float z, waveform waveForm, float wavelength, float supplement, float amplitude, float phase, ...)
{
	int i = W->xtoi(W, x), j = W->ytoj(W, y), k = W->ztok(W, z);

	if (W->xMinSur==SYM && !i) W->coskx = (F==Ex || F==Hy || F==Hz || F==iEx || F==iHy || F==iHz) ? 1 : -1;
	if (W->yMinSur==SYM && !j) W->cosky = (F==Ey || F==Hz || F==Hx || F==iEy || F==iHz || F==iHx) ? 1 : -1;
	if (W->zMinSur==SYM && !k) W->coskz = (F==Ez || F==Hx || F==Hy || F==iEz || F==iHx || F==iHy) ? 1 : -1;
	if (W->xMinSur==SYM && !W->coskx) W->coskx = (F==Ex || F==Hy || F==Hz || F==iEx || F==iHy || F==iHz) ? -1 : 1;
	if (W->yMinSur==SYM && !W->cosky) W->cosky = (F==Ey || F==Hz || F==Hx || F==iEy || F==iHz || F==iHx) ? -1 : 1;
	if (W->zMinSur==SYM && !W->coskz) W->coskz = (F==Ez || F==Hx || F==Hy || F==iEz || F==iHx || F==iHy) ? -1 : 1;

	if (waveForm) {
		if (!W->sourceType)
			W->sourceType = (W->Dom.x.Min==W->Dom.x.Max || W->Dom.y.Min==W->Dom.y.Max || W->Dom.z.Min==W->Dom.z.Max) ? 2 : 3;
		source S = createSource(W, F, waveForm, wavelength, supplement, amplitude);

		S->iMin = S->iMax = i, S->jMin = S->jMax = j, S->kMin = S->kMax = k;
		if (F==Ex || F==Hy || F==Hz || F==iEx || F==iHy || F==iHz) S->iMax++;
		if (F==Ey || F==Hz || F==Hx || F==iEy || F==iHz || F==iHx) S->jMax++;
		if (F==Ez || F==Hx || F==Hy || F==iEz || F==iHx || F==iHy) S->kMax++;

		if (!amplitude) amplitude = 1;
		if      (W->Dom.x.Min==W->Dom.x.Max) amplitude *= W->dx * sqrtf(2.5) * sqrtf(wavelength) / sq(W->dt);
		else if (W->Dom.y.Min==W->Dom.y.Max) amplitude *= W->dy * sqrtf(2.5) * sqrtf(wavelength) / sq(W->dt);
		else if (W->Dom.z.Min==W->Dom.z.Max) amplitude *= W->dz * sqrtf(2.5) * sqrtf(wavelength) / sq(W->dt);
		else amplitude *= sqrtf(1.5 / PI) * wavelength / sq(W->dt);

		S->phase = phase * PI / 180;
		S->Amplitude = makeField(S->iMin, S->iMax, S->jMin, S->jMax, S->kMin, S->kMax);
		if (F==Ex || F==iEx) setAmplitude(S, W->CH->dtdx, W->CE->dtdy, W->CE->dtdz, W->CE->x, amplitude*0.5);
		if (F==Ey || F==iEy) setAmplitude(S, W->CE->dtdx, W->CH->dtdy, W->CE->dtdz, W->CE->y, amplitude*0.5);
		if (F==Ez || F==iEz) setAmplitude(S, W->CE->dtdx, W->CE->dtdy, W->CH->dtdz, W->CE->z, amplitude*0.5);
		if (F==Hx || F==iHx) setAmplitude(S, W->CE->dtdx, W->CH->dtdy, W->CH->dtdz, W->CH->x, amplitude*0.25);
		if (F==Hy || F==iHy) setAmplitude(S, W->CH->dtdx, W->CE->dtdy, W->CH->dtdz, W->CH->y, amplitude*0.25);
		if (F==Hz || F==iHz) setAmplitude(S, W->CH->dtdx, W->CH->dtdy, W->CE->dtdz, W->CH->z, amplitude*0.25);
	}
	if (wavelength && W->T != wavelength/W->dt) W->f = 1/wavelength, W->T = wavelength/W->dt;
}


void (guidedWaveX)(world W, char *file, field F, float x, waveform waveForm, float wavelength, float supplement, float amplitude, float phase, ...)
{
	int i = W->xtoi(W, x);

	if (W->yMinSur==SYM) W->cosky = (F==Ey || F==Hz || F==Hx || F==iEy || F==iHz || F==iHx) ? 1 : -1;
	if (W->zMinSur==SYM) W->coskz = (F==Ez || F==Hx || F==Hy || F==iEz || F==iHx || F==iHy) ? 1 : -1;
	if (W->xMinSur==SYM && !i) W->coskx = (F==Ex || F==Hy || F==Hz || F==iEx || F==iHy || F==iHz) ? 1 : -1;
	if (W->xMinSur==SYM && !W->coskx) W->coskx = (F==Ex || F==Hy || F==Hz || F==iEx || F==iHy || F==iHz) ? -1 : 1;

	if (waveForm) {
		int d[2];
		float **f = h5read(d, file);
		if (!amplitude) amplitude = 1;
		source S = createSource(W, F, waveForm, wavelength, supplement, amplitude);

		if (d[0] != 1+W->jMax+(W->yMinSur==SYM?W->jMax:-W->jMin)) printf("The dimension of guided wave does not match!\n"), exit(0);
		if (d[1] != 1+W->kMax+(W->zMinSur==SYM?W->kMax:-W->kMin)) printf("The dimension of guided wave does not match!\n"), exit(0);
		S->phase = phase * PI / 180;
		S->iMin = S->iMax = i;
		S->jMin = W->jMin, S->jMax = W->jMax;
		S->kMin = W->kMin, S->kMax = W->kMax;
		S->Amplitude = makeField(S->iMin, S->iMax, S->jMin, S->jMax, S->kMin, S->kMax);

		#pragma omp parallel for
		for (int j=S->jMin; j<=S->jMax; j++) for (int k=S->kMin; k<=S->kMax; k++) {
			S->Amplitude[i][j][k] = amplitude * f[j-S->jMax+d[0]-1][k-S->kMax+d[1]-1];
		}
	}
	if (wavelength && W->T != wavelength/W->dt) W->f = 1/wavelength, W->T = wavelength/W->dt;
}


void (guidedWaveY)(world W, char *file, field F, float y, waveform waveForm, float wavelength, float supplement, float amplitude, float phase, ...)
{
	int j = W->ytoj(W, y);

	if (W->xMinSur==SYM) W->coskx = (F==Ex || F==Hy || F==Hz || F==iEx || F==iHy || F==iHz) ? 1 : -1;
	if (W->zMinSur==SYM) W->coskz = (F==Ez || F==Hx || F==Hy || F==iEz || F==iHx || F==iHy) ? 1 : -1;
	if (W->yMinSur==SYM && !j) W->cosky = (F==Ey || F==Hz || F==Hx || F==iEy || F==iHz || F==iHx) ? 1 : -1;
	if (W->yMinSur==SYM && !W->cosky) W->cosky = (F==Ey || F==Hz || F==Hx || F==iEy || F==iHz || F==iHx) ? -1 : 1;

	if (waveForm) {
		int d[2];
		float **f = h5read(d, file);
		if (!amplitude) amplitude = 1;
		source S = createSource(W, F, waveForm, wavelength, supplement, amplitude);

		if (d[0] != 1+W->iMax+(W->xMinSur==SYM?W->iMax:-W->iMin)) printf("The dimension of guided wave does not match!\n"), exit(0);
		if (d[1] != 1+W->kMax+(W->zMinSur==SYM?W->kMax:-W->kMin)) printf("The dimension of guided wave does not match!\n"), exit(0);
		S->phase = phase * PI / 180;
		S->jMin = S->jMax = j;
		S->iMin = W->iMin, S->iMax = W->iMax;
		S->kMin = W->kMin, S->kMax = W->kMax;
		S->Amplitude = makeField(S->iMin, S->iMax, S->jMin, S->jMax, S->kMin, S->kMax);

		#pragma omp parallel for
		for (int i=S->iMin; i<=S->iMax; i++) for (int k=S->kMin; k<=S->kMax; k++) {
			S->Amplitude[i][j][k] = amplitude * f[i-S->iMax+d[0]-1][k-S->kMax+d[1]-1];
		}
	}
	if (wavelength && W->T != wavelength/W->dt) W->f = 1/wavelength, W->T = wavelength/W->dt;
}


void (guidedWaveZ)(world W, char *file, field F, float z, waveform waveForm, float wavelength, float supplement, float amplitude, float phase, ...)
{
	int k = W->ztok(W, z);

	if (W->xMinSur==SYM) W->coskx = (F==Ex || F==Hy || F==Hz || F==iEx || F==iHy || F==iHz) ? 1 : -1;
	if (W->yMinSur==SYM) W->cosky = (F==Ey || F==Hz || F==Hx || F==iEy || F==iHz || F==iHx) ? 1 : -1;
	if (W->zMinSur==SYM && !k) W->coskz = (F==Ez || F==Hx || F==Hy || F==iEz || F==iHx || F==iHy) ? 1 : -1;
	if (W->zMinSur==SYM && !W->coskz) W->coskz = (F==Ez || F==Hx || F==Hy || F==iEz || F==iHx || F==iHy) ? -1 : 1;

	if (waveForm) {
		int d[2];
		float **f = h5read(d, file);
		if (!amplitude) amplitude = 1;
		source S = createSource(W, F, waveForm, wavelength, supplement, amplitude);

		if (d[0] != 1+W->iMax+(W->xMinSur==SYM?W->iMax:-W->iMin)) printf("The dimension of guided wave does not match!\n"), exit(0);
		if (d[1] != 1+W->jMax+(W->yMinSur==SYM?W->jMax:-W->jMin)) printf("The dimension of guided wave does not match!\n"), exit(0);
		S->phase = phase * PI / 180;
		S->kMin = S->kMax = k;
		S->iMin = W->iMin, S->iMax = W->iMax;
		S->jMin = W->jMin, S->jMax = W->jMax;
		S->Amplitude = makeField(S->iMin, S->iMax, S->jMin, S->jMax, S->kMin, S->kMax);

		#pragma omp parallel for
		for (int i=S->iMin; i<=S->iMax; i++) for (int j=S->jMin; j<=S->jMax; j++) {
			S->Amplitude[i][j][k] = amplitude * f[i-S->iMax+d[0]-1][j-S->jMax+d[1]-1];
		}
	}
	if (wavelength && W->T != wavelength/W->dt) W->f = 1/wavelength, W->T = wavelength/W->dt;
}


static float ***expandFieldX(world W, world SW, float ***g, int Alt)
{
	float ***f;
	MALLOC(f, W->iMIN-1, W->iNUM+2);
	for (int i=SW->iMIN-1; i<=SW->iMAX; i++) f[i] = g[i];
	for (int i=SW->iMAX+1; i<=W->iMAX+1; i++) f[i] = g[SW->iMAX+Alt];
	free(&g[SW->iMIN-1]);
	return f;
}


static float ***expandFieldY(world W, world SW, float ***g, int Alt)
{
	float ***f;
	MALLOC(f, SW->iMIN-1, SW->iNUM+2);
	MALLOC(f[SW->iMIN-1], W->jMIN-1, (SW->iNUM+2)*(W->jNUM+2));

	for (int i=SW->iMIN-1; i<=SW->iMAX+1; i++) {
		f[i] = f[SW->iMIN-1]+(W->jNUM+2)*(i-SW->iMIN+1);
		for (int j=SW->jMIN-1; j<=SW->jMAX; j++) f[i][j] = g[i][j];
		for (int j=SW->jMAX+1; j<=W->jMAX+1; j++) f[i][j] = g[i][SW->jMAX+Alt];
	}
	free(&g[SW->iMIN-1][SW->jMIN-1]);
	free(&g[SW->iMIN-1]);
	return f;
}


static void putObjectArrayTFSF(world W, tfsf S, matter *M, object *O)
{
	int N = 0;
	while (O[N].S) N++;

	for (int i=W->iMIN; i<=W->iMAX; i++) for (int j=W->jMIN; j<=W->jMAX; j++) for (int k=W->kMIN; k<=W->kMAX; k++) {
		if (i < S->iMin || S->iMax < i || j < S->jMin || S->jMax < j || k < S->kMin || S->kMax < k) {
			int n = 0;
			while (O[n].S && !(O[n].S(O[n], W->itox(W, i), W->jtoy(W, j), W->ktoz(W, k)))) n++;
			if (n < N) N = n;
		}
	}

	if (O[N].S || M[N].e.x != 1 || M[N].e.y && M[N].e.y != 1 || M[N].e.z && M[N].e.z != 1 || M[N].P->omega || M[N].P->gamma)
		putObjectArray(S->W, M+N, O+N);
}


tfsf (createTFSF)(world W, dom DOM, sur SUR)
{
	tfsf S;
	MALLOC(S, 0, 1);
	if (W->S) { tfsf s = W->S; while (s->S) s = s->S; s->S = S; }
	else W->S = S;

	dom Dom = W->Dom;
	sur Sur = SUR;

	if (SUR.x.MinSur==NIL) {
		Dom.x.Min = (W->xMinSur==SYM) ? 0 : W->xMIN;
		Dom.x.Max = (W->xMinSur==SYM) ? 0 : W->itox(W, 1+W->xtoi(W, W->xMIN));
		Sur.x.MinSur = (W->xMinSur==SYM) ? SYM : (W->xMinSur==BBC) ? BBC : PBC;
	}
	if (SUR.y.MinSur==NIL) {
		Dom.y.Min = (W->yMinSur==SYM) ? 0 : W->yMIN;
		Dom.y.Max = (W->yMinSur==SYM) ? 0 : W->jtoy(W, 1+W->ytoj(W, W->yMIN));
		Sur.y.MinSur = (W->yMinSur==SYM) ? SYM : (W->yMinSur==BBC) ? BBC : PBC;
	}
	if (SUR.z.MinSur==NIL) {
		Dom.z.Min = W->ktoz(W, W->ztok(W, W->zMin)-1);
		Dom.z.Max = W->ktoz(W, W->ztok(W, W->zMax)+1);
		Sur.z.MinSur = Sur.z.MaxSur = PML;
	}

	S->W = createWorld(Dom, W->Res, Sur, "");
	strcpy(S->W->ID, W->ID);

	if (SUR.y.MinSur==NIL) {
		for (int n=0; n<=S->W->complexField; n++) {
			S->W->E[n].x = expandFieldY(W, S->W, S->W->E[n].x, 0);
			S->W->E[n].y = expandFieldY(W, S->W, S->W->E[n].y, 1);
			S->W->E[n].z = expandFieldY(W, S->W, S->W->E[n].z, 0);
			S->W->H[n].x = expandFieldY(W, S->W, S->W->H[n].x, 1);
			S->W->H[n].y = expandFieldY(W, S->W, S->W->H[n].y, 0);
			S->W->H[n].z = expandFieldY(W, S->W, S->W->H[n].z, 1);
		}
	}
	if (SUR.x.MinSur==NIL) {
		for (int n=0; n<=S->W->complexField; n++) {
			S->W->E[n].x = expandFieldX(W, S->W, S->W->E[n].x, 1);
			S->W->E[n].y = expandFieldX(W, S->W, S->W->E[n].y, 0);
			S->W->E[n].z = expandFieldX(W, S->W, S->W->E[n].z, 0);
			S->W->H[n].x = expandFieldX(W, S->W, S->W->H[n].x, 0);
			S->W->H[n].y = expandFieldX(W, S->W, S->W->H[n].y, 1);
			S->W->H[n].z = expandFieldX(W, S->W, S->W->H[n].z, 1);
		}
	}

	if (!DOM.x.Max && DOM.x.Min>0 || SUR.x.MinSur==SYM) DOM.x.Min = -(DOM.x.Max = DOM.x.Max ? DOM.x.Max : DOM.x.Min);
	if (!DOM.y.Max && DOM.y.Min>0 || SUR.y.MinSur==SYM) DOM.y.Min = -(DOM.y.Max = DOM.y.Max ? DOM.y.Max : DOM.y.Min);
	if (!DOM.z.Max && DOM.z.Min>0 || SUR.z.MinSur==SYM) DOM.z.Min = -(DOM.z.Max = DOM.z.Max ? DOM.z.Max : DOM.z.Min);
	S->iMin = W->xMinSur==SYM ? 0 : DOM.x.Min == DOM.x.Max ? W->iMin : W->iMin > W->xtoi(W, DOM.x.Min) ? W->iMIN : W->xtoi(W, DOM.x.Min);
	S->jMin = W->yMinSur==SYM ? 0 : DOM.y.Min == DOM.y.Max ? W->jMin : W->jMin > W->ytoj(W, DOM.y.Min) ? W->jMIN : W->ytoj(W, DOM.y.Min);
	S->kMin = W->zMinSur==SYM ? 0 : DOM.z.Min == DOM.z.Max ? W->kMin : W->kMin > W->ztok(W, DOM.z.Min) ? W->kMIN : W->ztok(W, DOM.z.Min);
	S->iMax = DOM.x.Min == DOM.x.Max ? W->iMax : W->iMax < W->xtoi(W, DOM.x.Max) ? W->iMAX : W->xtoi(W, DOM.x.Max);
	S->jMax = DOM.y.Min == DOM.y.Max ? W->jMax : W->jMax < W->ytoj(W, DOM.y.Max) ? W->jMAX : W->ytoj(W, DOM.y.Max);
	S->kMax = DOM.z.Min == DOM.z.Max ? W->kMax : W->kMax < W->ztok(W, DOM.z.Max) ? W->kMAX : W->ztok(W, DOM.z.Max);
	S->iMIN = (S->iMin == W->iMIN) ? S->iMin : S->iMin + 1, S->iMAX = (S->iMax == W->iMAX) ? S->iMax + 1 : S->iMax;
	S->jMIN = (S->jMin == W->jMIN) ? S->jMin : S->jMin + 1, S->jMAX = (S->jMax == W->jMAX) ? S->jMax + 1 : S->jMax;
	S->kMIN = (S->kMin == W->kMIN) ? S->kMin : S->kMin + 1, S->kMAX = (S->kMax == W->kMAX) ? S->kMax + 1 : S->kMax;

	if (W->CE->M && W->CE->O) putObjectArrayTFSF(W, S, W->CE->M, W->CE->O);
	return S;
}


static void sourcePlane(world W, field F, int k, float kx, float ky, waveform waveForm, float wavelength, float supplement, float resolution, float amplitude, float phase)
{
	source S = createSource(W, F, waveForm, wavelength, supplement, resolution);

	S->amplitude = amplitude;
	S->phase = phase;
	S->iMin = W->iMIN, S->iMax = W->iMAX, S->jMin = W->jMIN, S->jMax = W->jMAX, S->kMin = k, S->kMax = k;

	if (W->complexField) {
		float di = (F==Ex||F==Hy||F==iEx||F==iHy) ? 0.5 : 0;
		float dj = (F==Ey||F==Hx||F==iEy||F==iHx) ? 0.5 : 0;

		S->Phase = makeField(S->iMin, S->iMax, S->jMin, S->jMax, S->kMin, S->kMax);
		#pragma omp parallel for
		for (int i=S->iMin; i<=S->iMax; i++) for (int j=S->jMin; j<=S->jMax; j++) {
			S->Phase[i][j][k] = phase + kx * W->itox(W, i-di) + ky * W->jtoy(W, j-dj);
		}
	}
}


static void sourcePlanes(world W, tfield pol(phaser,float,float), float theta, float phi, waveform waveForm, float wavelength, float supplement, float resolution, float amplitude, float phase)
{
	int kE = (90<theta && theta<270) ? W->kMin : (theta<90 || 270<theta) ? W->kMax : 0;
	int kH = (90<theta && theta<270) ? W->kMin : (theta<90 || 270<theta) ? W->kMax+1 : 0;
	float sinth = sinf(theta*PI/180) * sqrtf(W->CE->z?1/W->CE->z[W->iMin][W->jMin][kH]:1);
	float costh = cosf(theta*PI/180) * sqrtf(W->CE->z?1/W->CE->z[W->iMin][W->jMin][kH]:1);
	float sinph = sinf(phi*PI/180);
	float cosph = cosf(phi*PI/180);
	float kx = -2 * PI * sinth * cosph / wavelength;
	float ky = -2 * PI * sinth * sinph / wavelength;

	if (W->xMinSur == BBC) W->coskx = cosf(kx*(W->xMax-W->xMin)), W->sinkx = sinf(kx*(W->xMax-W->xMin));
	if (W->yMinSur == BBC) W->cosky = cosf(ky*(W->yMax-W->yMin)), W->sinky = sinf(ky*(W->yMax-W->yMin));
	if ((W->xMinSur == SYM) && (phi==90 || phi==270)) W->coskx = pol==Ppol ? -1 : 1;
	if ((W->xMinSur == SYM) && (phi==0  || phi==180) && (theta==0 || theta==180)) W->coskx = pol==Ppol ? 1 : -1;
	if ((W->yMinSur == SYM) && (phi==0  || phi==180)) W->cosky = pol==Ppol ? -1 : 1;
	if ((W->yMinSur == SYM) && (phi==90 || phi==270) && (theta==0 || theta==180)) W->cosky = pol==Ppol ? 1 : -1;

	if (waveForm) {
		float n = (W->CE->z ? 1/sqrtf(W->CE->z[W->iMin][W->jMin][kH]) : 1);
		float kzdz = W->CE->dz[kE] * PI * fabsf(costh) / wavelength;
		float ampE = W->CE->dtdz[kE] * (amplitude ? amplitude : 1) * 2 / ((1+n) * sqrtf(n) * n);
		float ampH = W->CH->dtdz[kH] * (amplitude ? amplitude : 1) * 2 / ((1+n) * sqrtf(n));

		if (pol == Ppol) {
			if (phi != 90 && phi != 270) sourcePlane(W, Ex, kE, kx, ky, waveForm, wavelength, supplement, resolution, ampE*cosph*costh, (phase+90)*PI/180);
			if (phi !=  0 && phi != 180) sourcePlane(W, Ey, kE, kx, ky, waveForm, wavelength, supplement, resolution, ampE*sinph*costh, (phase+90)*PI/180);
			if (phi !=  0 && phi != 180) sourcePlane(W, Hx, kH, kx, ky, waveForm, wavelength, supplement, resolution, ampH*sinph, (phase+90)*PI/180+kzdz);
			if (phi != 90 && phi != 270) sourcePlane(W, Hy, kH, kx, ky, waveForm, wavelength, supplement, resolution,-ampH*cosph, (phase+90)*PI/180+kzdz);

			if (W->complexField) {
				if (phi != 90 && phi != 270) sourcePlane(W, iEx, kE, kx, ky, waveForm, wavelength, supplement, resolution, ampE*cosph*costh, phase*PI/180);
				if (phi !=  0 && phi != 180) sourcePlane(W, iEy, kE, kx, ky, waveForm, wavelength, supplement, resolution, ampE*sinph*costh, phase*PI/180);
				if (phi !=  0 && phi != 180) sourcePlane(W, iHx, kH, kx, ky, waveForm, wavelength, supplement, resolution, ampH*sinph, phase*PI/180+kzdz);
				if (phi != 90 && phi != 270) sourcePlane(W, iHy, kH, kx, ky, waveForm, wavelength, supplement, resolution,-ampH*cosph, phase*PI/180+kzdz);
			}
		}

		if (pol == Spol) {
			if (phi !=  0 && phi != 180) sourcePlane(W, Ex, kE, kx, ky, waveForm, wavelength, supplement, resolution,-ampE*sinph, (phase+90)*PI/180);
			if (phi != 90 && phi != 270) sourcePlane(W, Ey, kE, kx, ky, waveForm, wavelength, supplement, resolution, ampE*cosph, (phase+90)*PI/180);
			if (phi != 90 && phi != 270) sourcePlane(W, Hx, kH, kx, ky, waveForm, wavelength, supplement, resolution, ampH*cosph*costh, (phase+90)*PI/180+kzdz);
			if (phi !=  0 && phi != 180) sourcePlane(W, Hy, kH, kx, ky, waveForm, wavelength, supplement, resolution, ampH*sinph*costh, (phase+90)*PI/180+kzdz);

			if (W->complexField) {
				if (phi !=  0 && phi != 180) sourcePlane(W, iEx, kE, kx, ky, waveForm, wavelength, supplement, resolution,-ampE*sinph, phase*PI/180);
				if (phi != 90 && phi != 270) sourcePlane(W, iEy, kE, kx, ky, waveForm, wavelength, supplement, resolution, ampE*cosph, phase*PI/180);
				if (phi != 90 && phi != 270) sourcePlane(W, iHx, kH, kx, ky, waveForm, wavelength, supplement, resolution, ampH*cosph*costh, phase*PI/180+kzdz);
				if (phi !=  0 && phi != 180) sourcePlane(W, iHy, kH, kx, ky, waveForm, wavelength, supplement, resolution, ampH*sinph*costh, phase*PI/180+kzdz);
			}
		}
	}
	if (wavelength && W->T != wavelength/W->dt) W->f = 1/wavelength, W->T = wavelength/W->dt;
}


void (planewave)(world W, dom Dom, tfield pol(phaser,float,float), float theta, float phi, waveform waveForm, float wavelength, float supplement, float amplitude, float phase, ...)
{
	while (theta < 0) theta += 360;
	while (theta >= 360) theta -= 360;
	while (phi < 0) phi += 360;
	while (phi >= 360) phi -= 360;

	sur Sur = {{NIL}, {NIL}, {NIL}};
	if ((theta !=0 && theta != 180) && (phi != 90 && phi != 270)) Sur.x.MinSur = BBC;
	if ((theta !=0 && theta != 180) && (phi !=  0 && phi != 180)) Sur.y.MinSur = BBC;

	float resolution = amplitude;
	if (!amplitude) amplitude = 1;

	tfsf S = createTFSF(W, Dom, Sur);
	sourcePlanes(W, pol, theta, phi, 0, wavelength, supplement, resolution, 0, 0);
	sourcePlanes(S->W, pol, theta, phi, waveForm, wavelength, supplement, resolution, amplitude, phase);
	W->N = S->W->N;
}


void beam(world W, field F, float x, float y, float z, float theta, waveform waveForm, float wavelength, float supplement, float spotsize, float amplitude, float phase, ...)
{
	while (theta < 0) theta += 360;
	while (theta >= 360) theta -= 360;

	if (W->xMinSur==SYM) W->coskx = (F==Ex || F==iEx) ? 1 : -1;
	if (W->yMinSur==SYM) W->cosky = (F==Ey || F==iEy) ? 1 : -1;

	if (waveForm) {
		int k = (90<theta && theta<270) ? W->kMin : (theta<90 || 270<theta) ? W->kMax : 0;
		float n  = W->CE->z ? 1/sqrtf(W->CE->z[W->xtoi(W,MIN(x,W->xMax))][W->ytoj(W,MIN(y,W->yMax))][k]) : 1;
		source S = createSource(W, F, waveForm, wavelength, supplement, amplitude);
		S->iMin = W->iMIN, S->iMax = W->iMAX, S->jMin = W->jMIN, S->jMax = W->jMAX, S->kMin = k, S->kMax = k;
		if (!amplitude) amplitude = 1;
		amplitude /= sqrtf(n*n*n);

		if (spotsize) {
			float kz = 2 * PI * n / wavelength;
			float zz = fabsf(z - W->ktoz(W,k));
			float w0 = wavelength * (spotsize ? spotsize : 1) / (n * sqrtf(-2.0 * logf(0.5)));
			float z0 = PI * n * sq(w0) / wavelength;
			float wz = w0 * sqrtf(1 + sq(zz/z0));
			float Rz = zz * (1 + sq(z0/zz));
			float nz = atanf(zz/z0);

			S->Amplitude = makeField(S->iMin, S->iMax, S->jMin, S->jMax, S->kMin, S->kMax);
			S->Phase = makeField(S->iMin, S->iMax, S->jMin, S->jMax, S->kMin, S->kMax);

			for (int i=S->iMin; i<=S->iMax; i++) for (int j=S->jMin; j<=S->jMax; j++) {
				float r2 = 0;
				if (x < INF) r2 += sq(x - W->itox(W,F==Ex?i-0.5:i));
				if (y < INF) r2 += sq(y - W->jtoy(W,F==Ey?j-0.5:j));
				S->Amplitude[i][j][k] = amplitude * w0/wz * expf(-r2/sq(wz));
				S->Phase[i][j][k] = phase * PI / 180 - kz * (0.5 * r2 / Rz) + nz;
			}
		}
		else {
			S->amplitude = amplitude;
			S->phase = phase * PI / 180;
		}
	}
	if (wavelength && W->T != wavelength/W->dt) W->f = 1/wavelength, W->T = wavelength/W->dt;
}


void (focusedBeam)(world W, dom Dom, field F, float x, float y, float z, float theta, waveform waveForm, float wavelength, float supplement, float spotsize, float amplitude, float phase, ...)
{
	if (Dom.x.Min==W->Dom.x.Min && Dom.x.Max==W->Dom.x.Max &&
		Dom.y.Min==W->Dom.y.Min && Dom.y.Max==W->Dom.y.Max &&
		Dom.z.Min==W->Dom.z.Min && Dom.z.Max==W->Dom.z.Max) {
		beam(W, F, x, y, z, theta, waveForm, wavelength, supplement, spotsize, amplitude, phase);
	}
	else {
		sur Sur = {{NIL}, {NIL}, {NIL}};
		if (x < INF) Sur.x.MinSur = W->Sur.x.MinSur, Sur.x.MaxSur = W->Sur.x.MaxSur;
		if (y < INF) Sur.y.MinSur = W->Sur.y.MinSur, Sur.y.MaxSur = W->Sur.y.MaxSur;
		tfsf S = createTFSF(W, Dom, Sur);
		beam(W, F, x, y, z, theta, 0, wavelength, supplement, 0, 0, 0);
		beam(S->W, F, x, y, z, theta, waveForm, wavelength, supplement, spotsize, amplitude, phase);
		W->N = S->W->N;
	}
}


static void removeSource(source S)
{
	if (S->S) removeSource(S->S);
	if (S->Amplitude) removeField(S->Amplitude, S->iMin, S->jMin, S->kMin);
	if (S->Phase) removeField(S->Phase, S->iMin, S->jMin, S->kMin);
	free(S);
}


static void removeTFSF(tfsf S)
{
	if (S->S) removeTFSF(S->S);
	deleteWorld(S->W);
	free(S);
}


void removeSources(world W)
{
	if (W->SE) removeSource(W->SE), W->SE = 0;
	if (W->SH) removeSource(W->SH), W->SH = 0;
	if (W->S) removeTFSF(W->S), W->S = 0;
	W->f = 0;
	W->T = 0;
	W->N = 0;
}
