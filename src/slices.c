/**********************************************************************
 * ALiS is Copyright (C) 2009-2019 by Ho-Seok Ee <hsee@kongju.ac.kr>. *
 * Redistribution and use with or without modification, are permitted *
 * under the terms of the Artistic License version 2.                 *
 **********************************************************************/

#include <sys/stat.h>
#include <png.h>
#include <string.h>
#include "alis_c.h"


typedef struct queue *queue;
struct queue {
	world W;
	char *ID;
	FILE *File;
	float max, Max, ***F, ***G;
	long int N;
	method M;
	queue Q;
};


static float ***bilinearInterpolation(world W, slice S, float ***F)
{
	float ***T, ***G = makeField(S->iMIN, S->iMAX, S->jMIN, S->jMAX, S->kMIN, S->kMAX);

	if (S->iMin != S->iMIN || S->iMax != S->iMAX) {
		T = F, F = G, G = T;
		#pragma omp parallel for
		for (int i=S->iMin; i<S->iMax; i++) {
			float x0 = W->itox(W, i), x1 = W->itox(W, i+1);
			int n0 = x0/W->dx+(i>0?1:0), n1 = x1/W->dx-(i<0?1:0);
			for (int j=S->jMin; j<=S->jMax; j++) for (int n = n0; n <= n1; n++) for (int k=S->kMin; k<=S->kMax; k++)
				F[n][j][k] = (G[i][j][k] * (x1 - n*W->dx) + G[i+1][j][k] * (n*W->dx - x0)) / (x1 - x0);
		}
	}

	if (S->jMin != S->jMIN || S->jMax != S->jMAX) {
		T = F, F = G, G = T;
		#pragma omp parallel for
		for (int j=S->jMin; j<S->jMax; j++) {
			float y0 = W->jtoy(W, j), y1 = W->jtoy(W, j+1);
			int n0 = y0/W->dy+(j>0?1:0), n1 = y1/W->dy-(j<0?1:0);
			for (int i=S->iMIN; i<=S->iMAX; i++) for (int n = n0; n <= n1; n++) for (int k=S->kMin; k<=S->kMax; k++)
				F[i][n][k] = (G[i][j][k] * (y1 - n*W->dy) + G[i][j+1][k] * (n*W->dy - y0)) / (y1 - y0);
		}
	}

	if (S->kMin != S->kMIN || S->kMax != S->kMAX) {
		T = F, F = G, G = T;
		#pragma omp parallel for
		for (int k=S->kMin; k<S->kMax; k++) {
			float z0 = W->ktoz(W, k), z1 = W->ktoz(W, k+1);
			int n0 = z0/W->dz+(k>0?1:0), n1 = z1/W->dz-(k<0?1:0);
			for (int i=S->iMIN; i<=S->iMAX; i++) for (int j=S->jMIN; j<=S->jMAX; j++) for (int n = n0; n <= n1; n++)
				F[i][j][n] = (G[i][j][k] * (z1 - n*W->dz) + G[i][j][k+1] * (n*W->dz - z0)) / (z1 - z0);
		}
	}

	removeField(G, S->iMIN, S->jMIN, S->kMIN);
	return F;
}


slice createSlice(world W, float xMin, float xMax, float yMin, float yMax, float zMin, float zMax)
{
	slice S;
	MALLOC(S, 0, 1);
	*S = (struct slice) {
		W->xtoi(W, xMin), W->xtoi(W, xMax),
		W->ytoj(W, yMin), W->ytoj(W, yMax),
		W->ztok(W, zMin), W->ztok(W, zMax),
	};
	if (W->Res.nug.x.m>1 || W->Res.nug.y.m>1 || W->Res.nug.z.m>1) S->interpolate = bilinearInterpolation;

	if (S->iMax > W->iMAX) S->iMax = W->iMAX;
	if (S->jMax > W->jMAX) S->jMax = W->jMAX;
	if (S->kMax > W->kMAX) S->kMax = W->kMAX;
	if (S->iMin < W->iMIN && W->xMinSur!=SYM) S->iMin = W->iMIN;
	if (S->jMin < W->jMIN && W->yMinSur!=SYM) S->jMin = W->jMIN;
	if (S->kMin < W->kMIN && W->zMinSur!=SYM) S->kMin = W->kMIN;
	if (S->iMin <-W->iMAX && W->xMinSur==SYM) S->iMin =-W->iMAX;
	if (S->jMin <-W->jMAX && W->yMinSur==SYM) S->jMin =-W->jMAX;
	if (S->kMin <-W->kMAX && W->zMinSur==SYM) S->kMin =-W->kMAX;
	S->iMIN = (S->iMin==S->iMax || !S->interpolate) ? S->iMin : MIN(S->iMin, W->itox(W, S->iMin)/W->dx);
	S->jMIN = (S->jMin==S->jMax || !S->interpolate) ? S->jMin : MIN(S->jMin, W->jtoy(W, S->jMin)/W->dy);
	S->kMIN = (S->kMin==S->kMax || !S->interpolate) ? S->kMin : MIN(S->kMin, W->ktoz(W, S->kMin)/W->dz);
	S->iMAX = (S->iMin==S->iMax || !S->interpolate) ? S->iMax : MAX(S->iMax, W->itox(W, S->iMax)/W->dx);
	S->jMAX = (S->jMin==S->jMax || !S->interpolate) ? S->jMax : MAX(S->jMax, W->jtoy(W, S->jMax)/W->dy);
	S->kMAX = (S->kMin==S->kMax || !S->interpolate) ? S->kMax : MAX(S->kMax, W->ktoz(W, S->kMax)/W->dz);
	S->iNUM = 1-S->iMIN+S->iMAX;
	S->jNUM = 1-S->jMIN+S->jMAX;
	S->kNUM = 1-S->kMIN+S->kMAX;
	S->F = makeField(S->iMIN, S->iMAX, S->jMIN, S->jMAX, S->kMIN, S->kMAX);

	if (!W->CE->x || !W->CE->y || !W->CE->z || (S->iMin!=S->iMax && S->jMin!=S->jMax && S->kMin!=S->kMax)) return S;

	float Max=0;
	S->C = makeField(S->iMIN, S->iMAX, S->jMIN, S->jMAX, S->kMIN, S->kMAX);
	#pragma omp parallel
	{
		float max=0;
		#pragma omp for
		for (int i=S->iMIN; i<=S->iMAX; i++) for (int j=S->jMIN; j<=S->jMAX; j++) for (int k=S->kMIN; k<=S->kMAX; k++)
			if (max < (S->C[i][j][k] = ContourRI(W, i, j, k))) max = S->C[i][j][k];
		#pragma omp critical
		if (Max < max) Max = max;
	}
	if (Max) {
		#pragma omp parallel for
		for (int i=S->iMIN; i<=S->iMAX; i++) for (int j=S->jMIN; j<=S->jMAX; j++) for (int k=S->kMIN; k<=S->kMAX; k++)
			S->C[i][j][k] = 0.5*S->C[i][j][k]/Max;
	}
	return S;
}


void deleteSlice(slice S)
{
	if (S->F) removeField(S->F, S->iMIN, S->jMIN, S->kMIN);
	if (S->C) removeField(S->C, S->iMIN, S->jMIN, S->kMIN);
	free(S);
}


void (h5)(world W, slice S, float Max, char *F, char *filename, unsigned long int* cmap)
{
	int dim[3] = {S->iMin==S->iMax?0:S->iNUM, S->jMin==S->jMax?0:S->jNUM, S->kMin==S->kMax?0:S->kNUM};

	if (S->interpolate) S->F = S->interpolate(W, S, S->F);
	h5write(dim, &S->F[S->iMIN][S->jMIN][S->kMIN], filename);
}


void (txt)(world W, slice S, float Max, char *F, char *filename, unsigned long int* cmap)
{
	int d=0, dim[3]={0};
	FILE *File = fopen(strcat(filename, ".txt"), "w");

	if (S->iMin != S->iMax) dim[d++] = S->iNUM;
	if (S->jMin != S->jMax) dim[d++] = S->jNUM;
	if (S->kMin != S->kMax) dim[d++] = S->kNUM;
	if (d==3) printf("Can't write 3D slice to txt file.\n"), exit(0);

	if (d==1) {
		if (S->iMin != S->iMax) for (int i=S->iMin; i<=S->iMax; i++)
			fprintf(File, "%f\t%e\n", W->itox(W, i), S->F[i][S->jMin][S->kMin]);
		if (S->jMin != S->jMax) for (int j=S->jMin; j<=S->jMax; j++)
			fprintf(File, "%f\t%e\n", W->jtoy(W, j), S->F[S->iMin][j][S->kMin]);
		if (S->kMin != S->kMax) for (int k=S->kMin; k<=S->kMax; k++)
			fprintf(File, "%f\t%e\n", W->ktoz(W, k), S->F[S->iMin][S->jMin][k]);
	}
	else {
		if (S->interpolate) S->F = S->interpolate(W, S, S->F);
		for (int j=dim[1]-1; j>=0; j--) {
			fprintf(File, "%e", S->F[S->iMIN][S->jMIN][S->kMIN+j]);
			for (int i=1; i<dim[0]; i++) fprintf(File, "\t%e", S->F[S->iMIN][S->jMIN][S->kMIN+dim[1]*i+j]);
			fprintf(File, "\n");
		}
	}
	fclose(File);
}


void (png)(world W, slice S, float Max, char *F, char *filename, unsigned long int* cmap)
{
	int d=0, dim[3]={1,1,1};
	FILE *File = fopen(strcat(filename, ".png"), "wb");
	float Min = 0;
	if (cmap[0] > LONG_MAX) Min = -Max;
	if (F[0]=='L' && F[1]=='o' && F[2]=='g' && F[3]!='R' && F[4]!='I') Min = Max-3;

	if (S->iMin != S->iMax) dim[d++] = S->iNUM;
	if (S->jMin != S->jMax) dim[d++] = S->jNUM;
	if (S->kMin != S->kMax) dim[d++] = S->kNUM;
	if (d==1) dim[2] = 64;
	if (d==3) printf("Can't write 3D slice to png file.\n"), exit(0);

	if (S->interpolate) S->F = S->interpolate(W, S, S->F);
	png_structp png = png_create_write_struct(PNG_LIBPNG_VER_STRING, 0, 0, 0);
	png_infop info = png_create_info_struct(png);
	png_init_io(png, File);
	png_set_compression_level(png, PNG_Z_DEFAULT_COMPRESSION);
	png_set_IHDR(png, info, dim[0], dim[1]*dim[2], 8,
		PNG_COLOR_TYPE_RGB, PNG_INTERLACE_NONE, PNG_COMPRESSION_TYPE_DEFAULT, PNG_FILTER_TYPE_DEFAULT);
//	png_set_text(png, info, (png_text[]){PNG_TEXT_COMPRESSION_NONE, "Title", "Value"}, 1);
	png_write_info(png, info);
	png_byte *row;
	MALLOC(row, 0, 3*dim[0]);

	for (int j=dim[1]-1; j>=0; j--) {
		#pragma omp parallel for
		for (int i=0; i<dim[0]; i++) {
			float C = !S->C ? 1 : 1-S->C[S->iMIN][S->jMIN][S->kMIN+dim[1]*i+j];
			float f = (Max-Min)?MAX(0,MIN(abs((long)cmap[0])-1,(abs((long)cmap[0])-1)*(S->F[S->iMIN][S->jMIN][S->kMIN+dim[1]*i+j]-Min)/(Max-Min))):0;
			int n = f;
			if (F[0]!='E' && F[0]!='H' && F[0]!='J' && F[0]!='U' && F[0]!='S' && F[0]!='i') C = 1;
			if (n > abs((long)cmap[0])-2) n = abs((long)cmap[0])-2;
			row[3*i]   = C * ((1-f+n) * (cmap[n+1]>>0x10&0xFF) + (f-n) * (cmap[n+2]>>0x10&0xFF));
			row[3*i+1] = C * ((1-f+n) * (cmap[n+1]>>0x08&0xFF) + (f-n) * (cmap[n+2]>>0x08&0xFF));
			row[3*i+2] = C * ((1-f+n) * (cmap[n+1]>>0x00&0xFF) + (f-n) * (cmap[n+2]>>0x00&0xFF));
		}
		for (int i=0; i<dim[2]; i++) png_write_row(png, row);
	}
	png_write_end(png, info);
	png_destroy_write_struct(&png, 0);
	fclose(File);
	free(row);
}


static queue searchQueue(world W, char *ID, method M, queue *staticQ)
{
	if (!*staticQ) {
		MALLOC(*staticQ, 0, 1);
		**staticQ = (struct queue) {0};
	}
	queue Q = *staticQ;

	while (Q->Q && (W != Q->W || M != Q->M || strcmp(ID, Q->ID))) Q = Q->Q;
	if (Q->Q) {
		strcat(strcpy(ID, W->ID), Q->ID);
	}
	else {
		Q->W = W;
		Q->M = M;
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


void (sliceSnap)(world W, field F, slice S, char *f, char *s, ...)
{
	static queue staticQ;
	char *format, string[1024], filename[1024];
	va_list ap;
	va_start(ap, s);
	int N = (F==RI || F==LogRI || F==ContourRI || F==Contour || F==raweEffx || F==raweEffy || F==raweEffz || F==eEff) ? 0 : va_arg(ap, int);
	method M = va_arg(ap, method);
	unsigned long int *cmap = (M == png) ? va_arg(ap, unsigned long int *) : 0;
	float Max = (M == png) ? va_arg(ap, double) : 1;
	int step = (!N || N>W->T) ? 1 : W->T/N;
	int skip = Max<0 ? ((int)(-Max*W->T/step))*step : 0;
	format = va_arg(ap, char*);
	vsprintf(filename, format, ap);
	if (format = strchr(filename,'%')) strcpy(string, format+1), strcat(strcat(strcat(strcpy(format, f), "-"), s), string);
	queue Q = searchQueue(W, filename, M, &staticQ);
	if (N) sprintf(string, "%05ld", (Q->N-skip)/step), strcat(filename, string);

	if (Q->N % step == 0) {
		float max=0;
		#pragma omp parallel
		{
			float m=0;
			#pragma omp for
			for (int i=S->iMin; i<=S->iMax; i++) for (int j=S->jMin; j<=S->jMax; j++) for (int k=S->kMin; k<=S->kMax; k++)
				if (m < fabsf(S->F[i][j][k] = F(W, i, j, k))) m = fabsf(S->F[i][j][k]);
			#pragma omp critical
			if (max < m) max = m;
		}
		if (!Max) Max = max;
		else if (Max<0) {
			Max = Q->Max;
			if (Q->max < max) Q->max = max;
			if (Q->N % skip == 0) Q->Max = Q->max, Q->max = 0;
		}
		if (Q->N > skip) M(W, S, Max?Max:1, f, filename, cmap);
	}
}


void (sliceTimeAvg)(world W, field F, slice S, char *f, char *s, ...)
{
	static queue staticQ;
	char *format, string[1024], filename[1024];
	va_list ap;
	va_start(ap, s);
	int N = (F==RI || F==LogRI || F==ContourRI || F==Contour || F==raweEffx || F==raweEffy || F==raweEffz || F==eEff) ? 0 : va_arg(ap, int);
	method M = va_arg(ap, method);
	unsigned long int *cmap = (M == png) ? va_arg(ap, unsigned long int *) : 0;
	float Max = (M == png) ? va_arg(ap, double) : 1;
	format = va_arg(ap, char*);
	vsprintf(filename, format, ap);
	if (format = strchr(filename,'%')) strcpy(string, format+1), strcat(strcat(strcat(strcpy(format, f), "-"), s), string);
	queue Q = searchQueue(W, filename, M, &staticQ);
	if (N) sprintf(string, "%05ld", Q->N / N), strcat(filename, string);

	if (Q->N == 1) {
		Q->F = makeField(S->iMin, S->iMax, S->jMin, S->jMax, S->kMin, S->kMax);
	}
	#pragma omp parallel for
	for (int i=S->iMin; i<=S->iMax; i++) for (int j=S->jMin; j<=S->jMax; j++) for (int k=S->kMin; k<=S->kMax; k++)
		Q->F[i][j][k] += F(W, i, j, k);

	if (Q->N % N == 0) {
		float max=0;
		#pragma omp parallel
		{
			float m=0;
			#pragma omp for
			for (int i=S->iMin; i<=S->iMax; i++) for (int j=S->jMin; j<=S->jMax; j++) for (int k=S->kMin; k<=S->kMax; k++)
				if (m < fabsf(S->F[i][j][k] = Q->F[i][j][k] / N)) m = fabsf(S->F[i][j][k]);
			#pragma omp critical
			if (max < m) max = m;
		}
		if (Max<=0) Max = max;
		M(W, S, Max?Max:1, f, filename, cmap);

		#pragma omp parallel for
		for (int i=S->iMin; i<=S->iMax; i++) for (int j=S->jMin; j<=S->jMax; j++) for (int k=S->kMin; k<=S->kMax; k++)
			Q->F[i][j][k] = 0;
	}
}


void (sliceFreqDom)(world W, field F, slice S, char *f, char *s, int N, float L, ...)
{
	int n=0;
	field *G[30] = {Ex, Ey, Ez, Hx, Hy, Hz, iEx, iEy, iEz, iHx, iHy, iHz, rawEx, rawEy, rawEz, rawHx, rawHy, rawHz, Jx, Jy, Jz, rawJx, rawJy, rawJz, ScEx, ScEy, ScEz, SoEx, SoEy, SoEz};
	while (n<30 && *F!=G[n]) n++;
	if (n==30) printf("\nThe %s field is not supported in sliceFreqDom().\n", f), exit(0);

	static queue staticQ;
	char *format, string[1024], filename[1024];
	va_list ap;
	va_start(ap, L);
	method M = va_arg(ap, method);
	unsigned long int *cmap = (M == png) ? va_arg(ap, unsigned long int *) : 0;
	float Max = (M == png) ? va_arg(ap, double) : 1;
	format = va_arg(ap, char*);
	vsprintf(filename, format, ap);
	if (format = strchr(filename,'%')) strcpy(string, format+1), strcat(strcat(strcat(strcpy(format, f), "-"), s), string);
	queue Q = searchQueue(W, filename, M, &staticQ);

	if (Q->N == 1) {
		Q->F = makeField(S->iMin, S->iMax, S->jMin, S->jMax, S->kMin, S->kMax);
		Q->G = makeField(S->iMin, S->iMax, S->jMin, S->jMax, S->kMin, S->kMax);
	}

	float real = cosf(2 * PI * Q->N * W->dt / L);
	float imag = sinf(2 * PI * Q->N * W->dt / L);
	#pragma omp parallel for
	for (int i=S->iMin; i<=S->iMax; i++) for (int j=S->jMin; j<=S->jMax; j++) for (int k=S->kMin; k<=S->kMax; k++) {
		Q->F[i][j][k] += F(W, i, j, k) * real;
		Q->G[i][j][k] += F(W, i, j, k) * imag;
	}

	if (Q->N % N == 0) {
		if (Max <= 0) {
			#pragma omp parallel
			{
				float max=0;
				#pragma omp for
				for (int i=S->iMin; i<=S->iMax; i++) for (int j=S->jMin; j<=S->jMax; j++) for (int k=S->kMin; k<=S->kMax; k++)
					if (max < sq(Q->F[i][j][k])+sq(Q->G[i][j][k])) max = sq(Q->F[i][j][k])+sq(Q->G[i][j][k]);
				#pragma omp critical
				if (Max < sqrtf(max) * 2 / Q->N) Max = sqrtf(max) * 2 / Q->N;
			}
		}

		for (int n=0; n<12; n++) {
			#pragma omp parallel for
			for (int i=S->iMin; i<=S->iMax; i++) for (int j=S->jMin; j<=S->jMax; j++) for (int k=S->kMin; k<=S->kMax; k++)
				S->F[i][j][k] = (Q->F[i][j][k] * cosf(n * PI / 6) + Q->G[i][j][k] * sinf(n * PI / 6)) * 2 / Q->N;
			sprintf(string, "%s%02d", filename, n+1);
			M(W, S, Max, f, string, cmap);
		}
	}
}


float sliceMax(world W, field F, slice S)
{
	float Max=0;
	#pragma omp parallel
	{
		float f, max=0;
		#pragma omp for
		for (int i=S->iMin; i<=S->iMax; i++) for (int j=S->jMin; j<=S->jMax; j++) for (int k=S->kMin; k<=S->kMax; k++)
			if (max < (f=fabsf(F(W, i, j, k)))) max = f;
		#pragma omp critical
		if (Max < max) Max = max;
	}
	return Max;
}


static float sliceSumZ(world W, field F, slice S, int i, int j)
{
	float Sum=0;
	for (int k=S->kMin+1; k<S->kMax; k++) Sum += W->CE->dz[k] * F(W, i, j, k);
	Sum += 0.5 * W->CE->dz[S->kMin] * F(W, i, j, S->kMin);
	Sum += 0.5 * W->CE->dz[S->kMax] * F(W, i, j, S->kMax);
	return Sum;
}

static float sliceSumY(world W, field F, slice S, int i)
{
	float Sum=0;
	for (int j=S->jMin+1; j<S->jMax; j++) Sum += W->CE->dy[j] * sliceSumZ(W, F, S, i, j);
	Sum += 0.5 * W->CE->dy[S->jMin] * sliceSumZ(W, F, S, i, S->jMin);
	Sum += 0.5 * W->CE->dy[S->jMax] * sliceSumZ(W, F, S, i, S->jMax);
	return Sum;
}

float sliceSum(world W, field F, slice S)
{
	float Sum=0;
	#pragma omp parallel for reduction(+:Sum)
	for (int i=S->iMin+1; i<S->iMax; i++) Sum += W->CE->dx[i] * sliceSumY(W, F, S, i);
	Sum += 0.5 * W->CE->dx[S->iMin] * sliceSumY(W, F, S, S->iMin);
	Sum += 0.5 * W->CE->dx[S->iMax] * sliceSumY(W, F, S, S->iMax);
	if (S->iMin==S->iMax) Sum /= W->CE->dx[S->iMin];
	if (S->jMin==S->jMax) Sum /= W->CE->dy[S->jMin];
	if (S->kMin==S->kMax) Sum /= W->CE->dz[S->kMin];
	return Sum;
}
