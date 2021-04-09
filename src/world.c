/**********************************************************************
 * ALiS is Copyright (C) 2009-2019 by Ho-Seok Ee <hsee@kongju.ac.kr>. *
 * Redistribution and use with or without modification, are permitted *
 * under the terms of the Artistic License version 2.                 *
 **********************************************************************/

#include <string.h>
#include <hdf5.h>
#include "alis_c.h"


static int xyztoijk(world W, float ijktoxyz(world, float), float x)
{
	int i=0, j;
	float y;
	for (j=15; j>=0; j--) {
		y = ijktoxyz(W, i);
		if (x > y) i+=1u<<j;
		if (x < y) i-=1u<<j;
	}
	y = fabsf(x-ijktoxyz(W, i));
	if (fabsf(x-ijktoxyz(W, i+1)) < y) i++;
	if (fabsf(x-ijktoxyz(W, i-1)) < y) i--;
	return i;
}

static float nonUniform(int m, float min, float max, float ds, float i)
{
	float d = 10;
	float iMin = min/ds-0.5, iMax = max/ds+0.5;
	float r = d/(d-1), rr;
	int n = floor((d*m-0.5*m-0.5)*logf(m)/(m-1)+0.5);

	do rr = r, r = rr - (powf(rr,n)-(d*m-0.5*m-0.5)*logf(rr)-1)/(n*powf(rr,n-1)-(d*m-0.5*m-0.5)/rr);
	while (fabsf(r-rr) > 1E-6);

	if (i>=iMax+n) return ds * (iMax + (d*m-0.5*m-0.5) + (i-iMax-n) * m);
	if (i<=iMin-n) return ds * (iMin - (d*m-0.5*m-0.5) + (i-iMin+n) * m);
	if (i>iMax) return ds * (iMax + (powf(r,i-iMax)-1) / logf(r));
	if (i<iMin) return ds * (iMin - (powf(r,iMin-i)-1) / logf(r));
	return ds * i;
}

static int xtoi(world W, float x) { return xyztoijk(W, W->itox, x); }
static int ytoj(world W, float y) { return xyztoijk(W, W->jtoy, y); }
static int ztok(world W, float z) { return xyztoijk(W, W->ktoz, z); }

static float itoxUniform(world W, float i) { return W->dx*i; }
static float jtoyUniform(world W, float j) { return W->dy*j; }
static float ktozUniform(world W, float k) { return W->dz*k; }

static float itoxNonUniform(world W, float i) { return nonUniform(W->Res.nug.x.m, W->Res.nug.x.Min, W->Res.nug.x.Max, W->dx, i); }
static float jtoyNonUniform(world W, float j) { return nonUniform(W->Res.nug.y.m, W->Res.nug.y.Min, W->Res.nug.y.Max, W->dy, j); }
static float ktozNonUniform(world W, float k) { return nonUniform(W->Res.nug.z.m, W->Res.nug.z.Min, W->Res.nug.z.Max, W->dz, k); }


#define removeField(F,i,j,k) free(&F[i][j][k]),free(&F[i][j]),free(&F[i]),F=0
float ***makeField(int iMin, int iMax, int jMin, int jMax, int kMin, int kMax)
{
	float ***f;
	MALLOC(f, iMin, 1-iMin+iMax);
	MALLOC(f[iMin], jMin, (1-iMin+iMax)*(1-jMin+jMax));
	MALLOC(f[iMin][jMin], kMin, (1-iMin+iMax)*(1-jMin+jMax)*(1-kMin+kMax));
	#pragma omp parallel for
	for (int i=iMin; i<=iMax; i++) {
		f[i] = f[iMin]+(1-jMin+jMax)*(i-iMin);
		for (int j=jMin; j<=jMax; j++) f[i][j] = f[iMin][jMin]+(1-kMin+kMax)*((1-jMin+jMax)*(i-iMin)+(j-jMin));
	}
	return f;
}


float complex ***makeComplexField(int iMin, int iMax, int jMin, int jMax, int kMin, int kMax)
{
	float complex ***f;
	MALLOC(f, iMin, 1-iMin+iMax);
	MALLOC(f[iMin], jMin, (1-iMin+iMax)*(1-jMin+jMax));
	MALLOC(f[iMin][jMin], kMin, (1-iMin+iMax)*(1-jMin+jMax)*(1-kMin+kMax));
	#pragma omp parallel for
	for (int i=iMin; i<=iMax; i++) {
		f[i] = f[iMin]+(1-jMin+jMax)*(i-iMin);
		for (int j=jMin; j<=jMax; j++) f[i][j] = f[iMin][jMin]+(1-kMin+kMax)*((1-jMin+jMax)*(i-iMin)+(j-jMin));
	}
	return f;
}


static void resetField(float ***f, int iMin, int iMax, int jMin, int jMax, int kMin, int kMax)
{
	#pragma omp parallel for
	for (int i=iMin; i<=iMax; i++) for (int j=jMin; j<=jMax; j++) {
		#pragma GCC ivdep
		for (int k=kMin; k<=kMax; k++) f[i][j][k] = 0;
	}
}


static void makeCpmlCoeff(sur S, int Offset, float **C0, float **C1, float **C2, float *dtds, float dt, int MinSur, int MaxSur, int Min, int Max, int MIN, int NUM)
{
	if (!S.pml.SigmaMax) S.pml.SigmaMax = 1;
	if (!S.pml.KappaMax) S.pml.KappaMax = S.pml.Num < 24 ? 1 : 10;
	if (!S.pml.AlphaMax) S.pml.AlphaMax = S.pml.Num < 24 ? 0 : 0.01;
	if (!S.pml.SigmaPower) S.pml.SigmaPower = 3.5;
	if (!S.pml.AlphaPower) S.pml.AlphaPower = 1;
	if (MinSur==PML || MaxSur==PML) {
		MALLOC(*C0, MIN, NUM);
		MALLOC(*C1, MIN, NUM);
		MALLOC(*C2, MIN, NUM);

		for (int n=1; n<=S.pml.Num; n++) {
			float Sigma = 0.8*S.pml.SigmaMax*(S.pml.SigmaPower+1)*powf((n-0.5*Offset)/S.pml.Num, S.pml.SigmaPower);
			float Kappa = 1+(S.pml.KappaMax-1)*powf((n-0.5*Offset)/S.pml.Num, S.pml.SigmaPower);
			float Alpha = S.pml.AlphaMax*powf((S.pml.Num-(n-0.5*Offset))/S.pml.Num, S.pml.AlphaPower);
			if (MinSur==PML) {
				dtds[Min-n-1+Offset] /= Kappa;
				float Ca = (1-0.5*(Sigma*dtds[Min-n-1+Offset]+Alpha*dt))/(1+0.5*(Sigma*dtds[Min-n-1+Offset]+Alpha*dt));
				float Cb = Sigma*dtds[Min-n-1+Offset]/(1+0.5*(Sigma*dtds[Min-n-1+Offset]+Alpha*dt));
				(*C0)[Min-n-1+Offset] = Ca;
				(*C1)[Min-n-1+Offset] = -dtds[Min-n-1+Offset]*Cb;
				(*C2)[Min-n-1+Offset] = +dtds[Min-n-1+Offset]*Cb;
			}
			if (MaxSur==PML) {
				dtds[Max+n+1] /= Kappa;
				float Ca = (1-0.5*(Sigma*dtds[Max+n+1]+Alpha*dt))/(1+0.5*(Sigma*dtds[Max+n+1]+Alpha*dt));
				float Cb = Sigma*dtds[Max+n+1]/(1+0.5*(Sigma*dtds[Max+n+1]+Alpha*dt));
				(*C0)[Max+n+1] = Ca;
				(*C1)[Max+n+1] = -dtds[Max+n+1]*Cb;
				(*C2)[Max+n+1] = +dtds[Max+n+1]*Cb;
			}
		}
	}
}


static coeffs *makeVFieldCoeff(world W, sur S, int Offset)
{
	coeffs *C;
	MALLOC(C, 0, 1);
	S.pml.Num = S.pml.Num ? S.pml.Num : 12;

	MALLOC(C->dx, W->xMinSur!=SYM ? W->iMIN : -W->iMAX, W->xMinSur!=SYM ? W->iNUM : 2*W->iNUM-1);
	MALLOC(C->dy, W->yMinSur!=SYM ? W->jMIN : -W->jMAX, W->yMinSur!=SYM ? W->jNUM : 2*W->jNUM-1);
	MALLOC(C->dz, W->zMinSur!=SYM ? W->kMIN : -W->kMAX, W->zMinSur!=SYM ? W->kNUM : 2*W->kNUM-1);
	MALLOC(C->dtdx, W->iMIN, W->iNUM);
	MALLOC(C->dtdy, W->jMIN, W->jNUM);
	MALLOC(C->dtdz, W->kMIN, W->kNUM);

	for (int i=W->iMIN; i<=W->iMAX; i++) C->dtdx[i] = W->dt/(C->dx[i]=(W->itox(W, i-0.5*(Offset-1))-W->itox(W, i-0.5*(Offset+1))));
	for (int j=W->jMIN; j<=W->jMAX; j++) C->dtdy[j] = W->dt/(C->dy[j]=(W->jtoy(W, j-0.5*(Offset-1))-W->jtoy(W, j-0.5*(Offset+1))));
	for (int k=W->kMIN; k<=W->kMAX; k++) C->dtdz[k] = W->dt/(C->dz[k]=(W->ktoz(W, k-0.5*(Offset-1))-W->ktoz(W, k-0.5*(Offset+1))));
	if (W->xMinSur==SYM) for (int i=1; i<=W->iMAX; i++) C->dx[-i] = C->dx[i];
	if (W->yMinSur==SYM) for (int j=1; j<=W->jMAX; j++) C->dy[-j] = C->dy[j];
	if (W->zMinSur==SYM) for (int k=1; k<=W->kMAX; k++) C->dz[-k] = C->dz[k];

	makeCpmlCoeff(S, Offset, &C->xP, &C->xPy, &C->xPz, C->dtdx, W->dt, W->xMinSur, W->xMaxSur, W->iMin, W->iMax, W->iMIN, W->iNUM);
	makeCpmlCoeff(S, Offset, &C->yP, &C->yPz, &C->yPx, C->dtdy, W->dt, W->yMinSur, W->yMaxSur, W->jMin, W->jMax, W->jMIN, W->jNUM);
	makeCpmlCoeff(S, Offset, &C->zP, &C->zPx, &C->zPy, C->dtdz, W->dt, W->zMinSur, W->zMaxSur, W->kMin, W->kMax, W->kMIN, W->kNUM);

	C->N = 0;
	return C;
}


static vfield *makeVectorField(world W, int Offset)
{
	vfield *F;
	MALLOC(F, 0, 1+W->complexField);

	for (int m=0; m<=W->complexField; m++) {
		F[m].x = makeField(W->iMIN-1, W->iMAX+1, W->jMIN-1, W->jMAX+1, W->kMIN-1, W->kMAX+1);
		F[m].y = makeField(W->iMIN-1, W->iMAX+1, W->jMIN-1, W->jMAX+1, W->kMIN-1, W->kMAX+1);
		F[m].z = makeField(W->iMIN-1, W->iMAX+1, W->jMIN-1, W->jMAX+1, W->kMIN-1, W->kMAX+1);

		if (W->xMinSur==PML) F[m].xMinPy = makeField(W->iMIN+Offset, W->iMin+Offset-2, W->jMIN, W->jMAX, W->kMIN, W->kMAX);
		if (W->xMinSur==PML) F[m].xMinPz = makeField(W->iMIN+Offset, W->iMin+Offset-2, W->jMIN, W->jMAX, W->kMIN, W->kMAX);
		if (W->yMinSur==PML) F[m].yMinPz = makeField(W->iMIN, W->iMAX, W->jMIN+Offset, W->jMin+Offset-2, W->kMIN, W->kMAX);
		if (W->yMinSur==PML) F[m].yMinPx = makeField(W->iMIN, W->iMAX, W->jMIN+Offset, W->jMin+Offset-2, W->kMIN, W->kMAX);
		if (W->zMinSur==PML) F[m].zMinPx = makeField(W->iMIN, W->iMAX, W->jMIN, W->jMAX, W->kMIN+Offset, W->kMin+Offset-2);
		if (W->zMinSur==PML) F[m].zMinPy = makeField(W->iMIN, W->iMAX, W->jMIN, W->jMAX, W->kMIN+Offset, W->kMin+Offset-2);

		if (W->xMaxSur==PML) F[m].xMaxPy = makeField(W->iMax+2, W->iMAX, W->jMIN, W->jMAX, W->kMIN, W->kMAX);
		if (W->xMaxSur==PML) F[m].xMaxPz = makeField(W->iMax+2, W->iMAX, W->jMIN, W->jMAX, W->kMIN, W->kMAX);
		if (W->yMaxSur==PML) F[m].yMaxPz = makeField(W->iMIN, W->iMAX, W->jMax+2, W->jMAX, W->kMIN, W->kMAX);
		if (W->yMaxSur==PML) F[m].yMaxPx = makeField(W->iMIN, W->iMAX, W->jMax+2, W->jMAX, W->kMIN, W->kMAX);
		if (W->zMaxSur==PML) F[m].zMaxPx = makeField(W->iMIN, W->iMAX, W->jMIN, W->jMAX, W->kMax+2, W->kMAX);
		if (W->zMaxSur==PML) F[m].zMaxPy = makeField(W->iMIN, W->iMAX, W->jMIN, W->jMAX, W->kMax+2, W->kMAX);
	}
	return F;
}


static world createEmptyWorld(res R, dom D, sur S)
{
	world W;
	MALLOC(W, 0, 1);

	if (R.nug.x.Min>0) R.nug.x.Max = R.nug.x.Min, R.nug.x.Min = -R.nug.x.Max;
	if (R.nug.y.Min>0) R.nug.y.Max = R.nug.y.Min, R.nug.y.Min = -R.nug.y.Max;
	if (R.nug.z.Min>0) R.nug.z.Max = R.nug.z.Min, R.nug.z.Min = -R.nug.z.Max;
	if (S.x.MinSur==PBC && S.x.MaxSur) S.x.MinSur=BBC;
	if (S.y.MinSur==PBC && S.y.MaxSur) S.y.MinSur=BBC;
	if (S.z.MinSur==PBC && S.z.MaxSur) S.z.MinSur=BBC;
	if (!S.pml.Num) S.pml.Num = 12;

	*W = (struct world) { D, R, S, 0, R.dt, R.d.x, R.d.y?R.d.y:R.d.x, R.d.z?R.d.z:R.d.x, R.spec.eV?R.spec.eV:1, INF };
	if (!D.x.Max && D.x.Min>0 || S.x.MinSur==SYM) D.x.Min = -(D.x.Max = D.x.Max ? D.x.Max : D.x.Min);
	if (!D.y.Max && D.y.Min>0 || S.y.MinSur==SYM) D.y.Min = -(D.y.Max = D.y.Max ? D.y.Max : D.y.Min);
	if (!D.z.Max && D.z.Min>0 || S.z.MinSur==SYM) D.z.Min = -(D.z.Max = D.z.Max ? D.z.Max : D.z.Min);

	W->itox = R.gridFunction.itox ? R.gridFunction.itox : R.nug.x.m > 1 ? itoxNonUniform : itoxUniform;
	W->jtoy = R.gridFunction.jtoy ? R.gridFunction.jtoy : R.nug.y.m > 1 ? jtoyNonUniform : jtoyUniform;
	W->ktoz = R.gridFunction.ktoz ? R.gridFunction.ktoz : R.nug.z.m > 1 ? ktozNonUniform : ktozUniform;
	W->xtoi = xtoi, W->ytoj = ytoj, W->ztok = ztok;

	W->xMinSur = S.x.MinSur ? S.x.MinSur : S.x.MinSur;
	W->yMinSur = S.y.MinSur ? S.y.MinSur : S.x.MinSur;
	W->zMinSur = S.z.MinSur ? S.z.MinSur : S.x.MinSur;
	W->xMaxSur = S.x.MinSur!=BBC && S.x.MaxSur ? S.x.MaxSur : S.x.MinSur ? S.x.MinSur : S.x.MinSur;
	W->yMaxSur = S.y.MinSur!=BBC && S.y.MaxSur ? S.y.MaxSur : S.y.MinSur ? S.y.MinSur : S.x.MinSur;
	W->zMaxSur = S.z.MinSur!=BBC && S.z.MaxSur ? S.z.MaxSur : S.z.MinSur ? S.z.MinSur : S.x.MinSur;

	W->iMin = W->xMinSur==SYM ? 0 : D.x.Min || D.x.Max ? W->xtoi(W, D.x.Min) : -1;
	W->jMin = W->yMinSur==SYM ? 0 : D.y.Min || D.y.Max ? W->ytoj(W, D.y.Min) : -1;
	W->kMin = W->zMinSur==SYM ? 0 : D.z.Min || D.z.Max ? W->ztok(W, D.z.Min) : -1;
	W->iMax = D.x.Min || D.x.Max ? W->xtoi(W, D.x.Max) : 1, W->iNum = 1-W->iMin+W->iMax;
	W->jMax = D.y.Min || D.y.Max ? W->ytoj(W, D.y.Max) : 1, W->jNum = 1-W->jMin+W->jMax;
	W->kMax = D.z.Min || D.z.Max ? W->ztok(W, D.z.Max) : 1, W->kNum = 1-W->kMin+W->kMax;

	W->iMIN = W->iMin-(W->xMinSur==PML ? S.pml.Num+1 : 0);
	W->jMIN = W->jMin-(W->yMinSur==PML ? S.pml.Num+1 : 0);
	W->kMIN = W->kMin-(W->zMinSur==PML ? S.pml.Num+1 : 0);
	W->iMAX = W->iMax+(W->xMaxSur==PML ? S.pml.Num+1 : 0), W->iNUM = 1-W->iMIN+W->iMAX;
	W->jMAX = W->jMax+(W->yMaxSur==PML ? S.pml.Num+1 : 0), W->jNUM = 1-W->jMIN+W->jMAX;
	W->kMAX = W->kMax+(W->zMaxSur==PML ? S.pml.Num+1 : 0), W->kNUM = 1-W->kMIN+W->kMAX;

	W->xMax = W->itox(W, W->iMax), W->xMin = W->xMinSur==SYM ? -W->xMax : W->itox(W, W->iMin), W->xSize = W->xMax - W->xMin,
	W->yMax = W->jtoy(W, W->jMax), W->yMin = W->yMinSur==SYM ? -W->yMax : W->jtoy(W, W->jMin), W->ySize = W->yMax - W->yMin;
	W->zMax = W->ktoz(W, W->kMax), W->zMin = W->zMinSur==SYM ? -W->zMax : W->ktoz(W, W->kMin), W->zSize = W->zMax - W->zMin;
	W->xMAX = W->itox(W, W->iMAX), W->xMIN = W->xMinSur==SYM ? -W->xMAX : W->itox(W, W->iMIN);
	W->yMAX = W->jtoy(W, W->jMAX), W->yMIN = W->yMinSur==SYM ? -W->yMAX : W->jtoy(W, W->jMIN);
	W->zMAX = W->ktoz(W, W->kMAX), W->zMIN = W->zMinSur==SYM ? -W->zMAX : W->ktoz(W, W->kMIN);

	W->coskx = S.x.MinSur==BBC && S.x.MaxSur ? cosf(2*PI*(W->xMax-W->xMin)/S.x.MaxSur) : S.x.MinSur==SYM ? 0 : 1;
	W->cosky = S.y.MinSur==BBC && S.y.MaxSur ? cosf(2*PI*(W->yMax-W->yMin)/S.y.MaxSur) : S.y.MinSur==SYM ? 0 : 1;
	W->coskz = S.z.MinSur==BBC && S.z.MaxSur ? cosf(2*PI*(W->zMax-W->zMin)/S.z.MaxSur) : S.z.MinSur==SYM ? 0 : 1;
	W->sinkx = S.x.MinSur==BBC && S.x.MaxSur ? sinf(2*PI*(W->xMax-W->xMin)/S.x.MaxSur) : 0;
	W->sinky = S.y.MinSur==BBC && S.y.MaxSur ? sinf(2*PI*(W->yMax-W->yMin)/S.y.MaxSur) : 0;
	W->sinkz = S.z.MinSur==BBC && S.z.MaxSur ? sinf(2*PI*(W->zMax-W->zMin)/S.z.MaxSur) : 0;

	W->coskyp = S.y.MinSur==HBC && S.y.MaxSur ? cosf(PI*(2*(W->yMax-W->yMin)/S.y.MaxSur+(W->xMax-W->xMin)/S.x.MaxSur)) : 1;
	W->coskym = S.y.MinSur==HBC && S.y.MaxSur ? cosf(PI*(2*(W->yMax-W->yMin)/S.y.MaxSur-(W->xMax-W->xMin)/S.x.MaxSur)) : 1;
	W->sinkyp = S.y.MinSur==HBC && S.y.MaxSur ? sinf(PI*(2*(W->yMax-W->yMin)/S.y.MaxSur+(W->xMax-W->xMin)/S.x.MaxSur)) : 0;
	W->sinkym = S.y.MinSur==HBC && S.y.MaxSur ? sinf(PI*(2*(W->yMax-W->yMin)/S.y.MaxSur-(W->xMax-W->xMin)/S.x.MaxSur)) : 0;

	W->xyArea = W->xSize * W->ySize, W->xzArea = W->xSize * W->zSize, W->yzArea = W->ySize * W->zSize;
	if (S.x.MinSur==BBC || S.y.MinSur==BBC || S.z.MinSur==BBC) W->complexField = 1;
	return W;
}


world createWorld(dom D, res R, sur S, char *format, ...)
{
	int iNum=0, jNum=0;
	va_list ap;
	world W = createEmptyWorld(R, D, S);

	W->E = makeVectorField(W, 0);
	W->H = makeVectorField(W, 1);
	W->CE = makeVFieldCoeff(W, S, 0);
	W->CH = makeVFieldCoeff(W, S, 1);

	va_start(ap, format);
	vsprintf(W->ID, format, ap);
	va_end(ap);
	if (W->ID[0]) {
		printf("Deleting previous files.");
		exec("rm -fr %s %s.c?* %s[#--0-~]* %s.[#-bd-~]* 2> /dev/null", W->ID, W->ID, W->ID, W->ID);
	}
	return W;
}


static void makeAuxiliaryFields(world W, vfield *F, coeffs *C)
{
	for (int m=0; m<=W->complexField; m++) {
		MALLOC(F[m].Jx, 0, C->N);
		MALLOC(F[m].Jy, 0, C->N);
		MALLOC(F[m].Jz, 0, C->N);
		MALLOC(F[m].Kx, 0, C->N);
		MALLOC(F[m].Ky, 0, C->N);
		MALLOC(F[m].Kz, 0, C->N);
	}

	for (int n=0; n<C->N; n++) {
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
}


world cloneWorld(world SW)
{
	va_list ap;
	world W;
	MALLOC(W, 0, 1);
	*W = *SW;

	W->E = makeVectorField(W, 0);
	W->H = makeVectorField(W, 1);
	makeAuxiliaryFields(W, W->E, W->CE);
	makeAuxiliaryFields(W, W->H, W->CH);

	return W;
}


static void resetVectorField(world W, vfield *F, coeffs *C, int Offset)
{
	for (int m=0; m<=W->complexField; m++) {
		resetField(F[m].x, W->iMIN-1, W->iMAX+1, W->jMIN-1, W->jMAX+1, W->kMIN-1, W->kMAX+1);
		resetField(F[m].y, W->iMIN-1, W->iMAX+1, W->jMIN-1, W->jMAX+1, W->kMIN-1, W->kMAX+1);
		resetField(F[m].z, W->iMIN-1, W->iMAX+1, W->jMIN-1, W->jMAX+1, W->kMIN-1, W->kMAX+1);

		if (F[m].xMinPy) resetField(F[m].xMinPy, W->iMIN+Offset, W->iMin+Offset-2, W->jMIN, W->jMAX, W->kMIN, W->kMAX);
		if (F[m].xMinPz) resetField(F[m].xMinPz, W->iMIN+Offset, W->iMin+Offset-2, W->jMIN, W->jMAX, W->kMIN, W->kMAX);
		if (F[m].yMinPz) resetField(F[m].yMinPz, W->iMIN, W->iMAX, W->jMIN+Offset, W->jMin+Offset-2, W->kMIN, W->kMAX);
		if (F[m].yMinPx) resetField(F[m].yMinPx, W->iMIN, W->iMAX, W->jMIN+Offset, W->jMin+Offset-2, W->kMIN, W->kMAX);
		if (F[m].zMinPx) resetField(F[m].zMinPx, W->iMIN, W->iMAX, W->jMIN, W->jMAX, W->kMIN+Offset, W->kMin+Offset-2);
		if (F[m].zMinPy) resetField(F[m].zMinPy, W->iMIN, W->iMAX, W->jMIN, W->jMAX, W->kMIN+Offset, W->kMin+Offset-2);

		if (F[m].xMaxPy) resetField(F[m].xMaxPy, W->iMax+2, W->iMAX, W->jMIN, W->jMAX, W->kMIN, W->kMAX);
		if (F[m].xMaxPz) resetField(F[m].xMaxPz, W->iMax+2, W->iMAX, W->jMIN, W->jMAX, W->kMIN, W->kMAX);
		if (F[m].yMaxPz) resetField(F[m].yMaxPz, W->iMIN, W->iMAX, W->jMax+2, W->jMAX, W->kMIN, W->kMAX);
		if (F[m].yMaxPx) resetField(F[m].yMaxPx, W->iMIN, W->iMAX, W->jMax+2, W->jMAX, W->kMIN, W->kMAX);
		if (F[m].zMaxPx) resetField(F[m].zMaxPx, W->iMIN, W->iMAX, W->jMIN, W->jMAX, W->kMax+2, W->kMAX);
		if (F[m].zMaxPy) resetField(F[m].zMaxPy, W->iMIN, W->iMAX, W->jMIN, W->jMAX, W->kMax+2, W->kMAX);

		for (int n=0; n<C->N; n++) {
			for (int p=0; p<C->NK[n]; p++) {
				resetField(F[m].Jx[n][p], C->iMin[n], C->iMax[n], C->jMin[n], C->jMax[n], C->kMin[n], C->kMax[n]);
				resetField(F[m].Jy[n][p], C->iMin[n], C->iMax[n], C->jMin[n], C->jMax[n], C->kMin[n], C->kMax[n]);
				resetField(F[m].Jz[n][p], C->iMin[n], C->iMax[n], C->jMin[n], C->jMax[n], C->kMin[n], C->kMax[n]);
			}
			for (int p=C->NJ[n]; p<C->NK[n]; p++) {
				resetField(F[m].Kx[n][p], C->iMin[n], C->iMax[n], C->jMin[n], C->jMax[n], C->kMin[n], C->kMax[n]);
				resetField(F[m].Ky[n][p], C->iMin[n], C->iMax[n], C->jMin[n], C->jMax[n], C->kMin[n], C->kMax[n]);
				resetField(F[m].Kz[n][p], C->iMin[n], C->iMax[n], C->jMin[n], C->jMax[n], C->kMin[n], C->kMax[n]);
			}
		}
	}
}


void rebootWorld(world W)
{
	removeSources(W);
	resetVectorField(W, W->E, W->CE, 0);
	resetVectorField(W, W->H, W->CH, 1);
	W->t = W->E->t = 0;
	W->H->t = 0.5;
	W->n = 0;
	printf("\n");
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


static void removeVectorField(world W, vfield *F, int Offset)
{
	for (int m=0; m<=W->complexField; m++) {
		removeField(F[m].x, W->iMIN-1, W->jMIN-1, W->kMIN-1);
		removeField(F[m].y, W->iMIN-1, W->jMIN-1, W->kMIN-1);
		removeField(F[m].z, W->iMIN-1, W->jMIN-1, W->kMIN-1);

		if (F[m].xMinPy) removeField(F[m].xMinPy, W->iMIN+Offset, W->jMIN, W->kMIN);
		if (F[m].xMinPz) removeField(F[m].xMinPz, W->iMIN+Offset, W->jMIN, W->kMIN);
		if (F[m].yMinPz) removeField(F[m].yMinPz, W->iMIN, W->jMIN+Offset, W->kMIN);
		if (F[m].yMinPx) removeField(F[m].yMinPx, W->iMIN, W->jMIN+Offset, W->kMIN);
		if (F[m].zMinPx) removeField(F[m].zMinPx, W->iMIN, W->jMIN, W->kMIN+Offset);
		if (F[m].zMinPy) removeField(F[m].zMinPy, W->iMIN, W->jMIN, W->kMIN+Offset);

		if (F[m].xMaxPy) removeField(F[m].xMaxPy, W->iMax+2, W->jMIN, W->kMIN);
		if (F[m].xMaxPz) removeField(F[m].xMaxPz, W->iMax+2, W->jMIN, W->kMIN);
		if (F[m].yMaxPz) removeField(F[m].yMaxPz, W->iMIN, W->jMax+2, W->kMIN);
		if (F[m].yMaxPx) removeField(F[m].yMaxPx, W->iMIN, W->jMax+2, W->kMIN);
		if (F[m].zMaxPx) removeField(F[m].zMaxPx, W->iMIN, W->jMIN, W->kMax+2);
		if (F[m].zMaxPy) removeField(F[m].zMaxPy, W->iMIN, W->jMIN, W->kMax+2);
	}
	free(F);
}


static void removeVFieldCoeff(world W, coeffs *C)
{
	if (C->xP) free(&C->xP[W->iMIN]), free(&C->xPy[W->iMIN]), free(&C->xPz[W->iMIN]);
	if (C->yP) free(&C->yP[W->jMIN]), free(&C->yPz[W->jMIN]), free(&C->yPx[W->jMIN]);
	if (C->zP) free(&C->zP[W->kMIN]), free(&C->zPx[W->kMIN]), free(&C->zPy[W->kMIN]);
	free(&C->dtdx[W->iMIN]), free(&C->dx[W->xMinSur!=SYM ? W->iMIN :-W->iMAX]);
	free(&C->dtdy[W->jMIN]), free(&C->dy[W->yMinSur!=SYM ? W->jMIN :-W->jMAX]);
	free(&C->dtdz[W->kMIN]), free(&C->dz[W->zMinSur!=SYM ? W->kMIN :-W->kMAX]);
	free(C);
}


void deleteWorld(world W)
{
	removeSources(W);
	removeObjects(W);
	removeVFieldCoeff(W, W->CE);
	removeVFieldCoeff(W, W->CH);
	removeVectorField(W, W->E, 0);
	removeVectorField(W, W->H, 1);
	free(W);
}


void deleteCloneWorld(world W)
{
	removeSources(W);
	removeObjects(W);
	removeAuxiliaryFields(W, W->E, W->CE);
	removeAuxiliaryFields(W, W->H, W->CH);
	removeVectorField(W, W->E, 0);
	removeVectorField(W, W->H, 1);
	free(W);
}


void *h5read(int *dim, char *format, ...)
{
	void* p;
	va_list ap;
	char name[1024];
	hsize_t dims[3], maxdims[3];
	va_start(ap, format);
	vsprintf(name, format, ap);

	hid_t file = H5Fopen(strcat(name, ".h5"), H5F_ACC_RDONLY, H5P_DEFAULT);
	H5Gget_objname_by_idx(file, 0, name, 64);
	hid_t dataset = H5Dopen(file, name, H5P_DEFAULT);
	hid_t space = H5Dget_space(dataset);
	int N = H5Sget_simple_extent_ndims(space);
	H5Sget_simple_extent_dims(space, dims, maxdims);
	if (dim) for (int n=0; n<N; n++) dim[n] = dims[n];

	if (N == 1) {
		float *f;
		f = calloc(dims[0], sizeof(float));
		H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, f);
		p = f;
	}

	if (N == 2) {
		float **f;
		f = calloc(dims[0], sizeof(float*));
		f[0] = calloc(dims[0]*dims[1], sizeof(float));
		for (int i=1; i<dims[0]; i++) f[i] = f[0]+dims[1]*i;
		H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, f[0]);
		p = f;
	}

	if (N == 3) {
		float ***f;
		f = calloc(dims[0], sizeof(float**));
		f[0] = calloc(dims[0]*dims[1], sizeof(float*));
		f[0][0] = calloc(dims[0]*dims[1]*dims[2], sizeof(float));
		for (int i=1; i<dims[0]; i++) f[i] = f[0]+dims[1]*i;
		for (int i=0; i<dims[0]; i++) for (int j=1; j<dims[1]; j++) f[i][j] = f[i][0]+dims[2]*j;
		H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, f[0][0]);
		p = f;
	}

	H5Dclose(dataset);
	H5Sclose(space);
	H5Fclose(file);
	return p;
}


void h5write(int *dim, float *p, char *format, ...)
{
	va_list ap;
	char filename[1024];
	int d = 0;
	hsize_t dims[3] = {0};
	va_start(ap, format);
	vsprintf(filename, format, ap);

	for (int n=0; n<3; n++) if (dim[n]) dims[d++] = dim[n];
	hid_t file = H5Fcreate(strcat(filename, ".h5"), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
	hid_t space = H5Screate_simple(d, dims, 0);
	hid_t dataset = H5Dcreate(file, "data", H5T_NATIVE_FLOAT, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	H5Dwrite(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, p);
	H5Dclose(dataset);
	H5Sclose(space);
	H5Fclose(file);
}
