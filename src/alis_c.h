/**********************************************************************
 * ALiS is Copyright (C) 2009-2019 by Ho-Seok Ee <hsee@kongju.ac.kr>. *
 * Redistribution and use with or without modification, are permitted *
 * under the terms of the Artistic License version 2.                 *
 **********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <math.h>
#include <values.h>
#include <complex.h>
#include <ieee754.h>

#define NIL 6
#define BBC 5
#define PBC 4
#define PML 3
#define PEC 2
#define SYM 1
#define MAXPOLES 6

#define INF FLT_MAX
#define PI 3.14159265358979323846f

#define sq(x) ((x)*(x))
#define n(n) ((matter){(n)*(n)})
#define nk(n,k,f) ((matter){{(n)*(n)-(k)*(k)},{0,2*n*k*((float)f)}})
#define png(cmap,max) png,cmap,(double)max
#define objects (object[])
#define MIN(x,y) ((x)<(y)?(x):(y))
#define MAX(x,y) ((x)>(y)?(x):(y))
#define MALLOC(p,m,n) if ((p)=calloc((n),sizeof(*(p)))) (p)-=(m); else printf("Insufficient memory!\n"), exit(1)
#define INT2NAN(i) ((union ieee754_double){.ieee.exponent=-1,.ieee.mantissa1=i}.d)
#define NAN2INT(n) ((union ieee754_double){.ieee.exponent=-1}.ieee.exponent==(union ieee754_double){.d=n}.ieee.exponent&&!(union ieee754_double){.d=n}.ieee.mantissa0?(union ieee754_double){.d=n}.ieee.mantissa1:0)


/* world.c */

typedef struct world *world;
typedef struct matter matter;
typedef struct object object;
typedef struct source *source;
typedef struct tfsf *tfsf;
typedef struct slice *slice;
typedef struct phaser *phaser;
typedef float field(world, int, int, int);

typedef struct {
	float dt;
	struct {
		float x, y, z;
	} d;
	struct {
		float eV;
		int subN;
	} spec;
	struct {
		struct {
			int m;
			float Min, Max;
		} x, y, z;
	} nug;
	struct {
		float (*itox)(world, float);
		float (*jtoy)(world, float);
		float (*ktoz)(world, float);
	} gridFunction;
} res;

typedef struct {
	struct {
		float Min, Max;
	} x, y, z;
} dom;

typedef struct {
	struct {
		float MinSur, MaxSur;
	} x, y, z;
	struct {
		int Num;
		float SigmaMax, KappaMax, AlphaMax, SigmaPower, AlphaPower;
	} pml;
} sur;

typedef struct {
	float ***x, ***y, ***z;
	float ***yMinPx, ***zMinPx, ***xMinPy, ***zMinPy, ***xMinPz, ***yMinPz;
	float ***yMaxPx, ***zMaxPx, ***xMaxPy, ***zMaxPy, ***xMaxPz, ***yMaxPz;
	float *****Jx, *****Jy, *****Jz, *****Kx, *****Ky, *****Kz;
	float t;
} vfield;

typedef struct {
	float *dx, *dy, *dz, *dtdx, *dtdy, *dtdz;
	float *xP, *yP, *zP, *xPy, *xPz, *yPz, *yPx, *zPx, *zPy;
	float ***x, ***y, ***z;
	float ****Jx, ****Jy, ****Jz;
	float *ax, *ay, *az;
	float **JJx, **JJy, **JJz, **JKx, **JKy, **JKz, **JEx, **JEy, **JEz;
	float **KJx, **KJy, **KJz, **KKx, **KKy, **KKz, **KEx, **KEy, **KEz;
	int N, *NJ, *NK, *iMin, *iMax, *jMin, *jMax, *kMin, *kMax;
	object *O;
	matter *M;
} coeffs;

struct world {
	dom Dom;
	res Res;
	sur Sur;
	float t, dt, dx, dy, dz, eV, f;
	float xMin, xMax, yMin, yMax, zMin, zMax;
	float xMIN, xMAX, yMIN, yMAX, zMIN, zMAX;
	float coskx, sinkx, cosky, sinky, coskz, sinkz;
	float xSize, ySize, zSize, xyArea, xzArea, yzArea;
	int iMin, iMax, iNum, jMin, jMax, jNum, kMin, kMax, kNum;
	int iMIN, iMAX, iNUM, jMIN, jMAX, jNUM, kMIN, kMAX, kNUM;
	int xMinSur, xMaxSur, yMinSur, yMaxSur, zMinSur, zMaxSur;
	int subN, complexField, sourceType;
	unsigned long int n, N, T;
	int (*xtoi)(world, float), (*ytoj)(world, float), (*ztok)(world, float);
	float (*itox)(world, float), (*jtoy)(world, float), (*ktoz)(world, float);
	char ID[128];
	vfield *E, *H;
	coeffs *CE, *CH;
	source SE, SH;
	tfsf S;
};

world createWorld(dom, res, sur, char*, ...);
world cloneWorld(world);
void rebootWorld(world W);
void deleteWorld(world W);
void deleteCloneWorld(world W);
float ***makeField(int, int, int, int, int, int);
float complex ***makeComplexField(int, int, int, int, int, int);
#define removeField(F,i,j,k) free(&F[i][j][k]),free(&F[i][j]),free(&F[i]),F=0
void *h5read(int*, char*, ...);
void h5write(int*, float*, char*, ...);


/* object.c and shapes.c */

typedef int shape(object, float, float, float);

struct matter {
	struct {
		float x, y, z;
	} e;
	struct {
		float omega, gamma, freal, fimag;
	} P[MAXPOLES];
	float d[3][6];
};

struct object {
	shape *S;
	float f[3][3];
	void *p;
};

void putObjectArray(world, matter*, object*);
void putObjects(world, ...);
#define putObjects(...) ((putObjects)(__VA_ARGS__,(object){0}))
void removeObjects(world);
matter LT(matter, float);
matter nkFile(char*, float, float);
shape Box, Ball, Diamond;
shape RodX, RodY, RodZ, DRodX, DRodY, DRodZ, BiconeX, BiconeY, BiconeZ;
shape HRodXY, HRodXZ, HRodYZ, HRodYX, HRodZX, HRodZY;
shape Translation, Rotation;
shape Combination, Intersection, Difference;
shape Lattice, DLattice, CLatticeX, CLatticeY, CLatticeZ;
shape ImportContour;


/* update.c */

void copyBoundary(float***, float***, float, float, int, int, int, int, int, int, int, int, int);
#define copyBoundary(W,F,x,...) ((copyBoundary)(W->F[0].x,W->complexField?W->F[1].x:0,__VA_ARGS__))
void updateE(world);
void updateH(world);
void nonLinearSHG(world, world, float);
void nonLinearTWM(world, world, world, float);
int timer(long int, long int, char*, ...);
#define timer(...) ((timer)(__VA_ARGS__,""))


/* fields.c, slices.c and output.c */

typedef void (*method)(world, slice, float, char*, char*, unsigned long int*);
typedef unsigned long int colormap[];
struct slice {
	int iMin, iMax, jMin, jMax, kMin, kMax;
	int iMIN, iMAX, iNUM, jMIN, jMAX, jNUM, kMIN, kMAX, kNUM;
	float ***C, ***F;
	float ***(*interpolate)(world, slice, float***);
};

void exec(char*, ...);
void writeTxt(world, char*, ...);
void writeRow(world, char*, ...);
void writeSpectrum(world, int, float, float, char*, ...);
#define writeTxt(...) ((writeTxt)(__VA_ARGS__,""))
#define writeRow(...) ((writeRow)(__VA_ARGS__,INT2NAN(1)))
#define writeSpectrum(...) ((writeSpectrum)(__VA_ARGS__,INT2NAN(1)))
slice createSlice(world, float, float, float, float, float, float);
#define createSliceX(W,yCut,zCut) createSlice(W, W->xMin, W->xMax, yCut, yCut, zCut, zCut)
#define createSliceY(W,xCut,zCut) createSlice(W, xCut, xCut, W->yMin, W->yMax, zCut, zCut)
#define createSliceZ(W,xCut,yCut) createSlice(W, xCut, xCut, yCut, yCut, W->zMin, W->zMax)
#define createSliceXY(W,zCut) createSlice(W, W->xMin, W->xMax, W->yMin, W->yMax, zCut, zCut)
#define createSliceXZ(W,yCut) createSlice(W, W->xMin, W->xMax, yCut, yCut, W->zMin, W->zMax)
#define createSliceYZ(W,xCut) createSlice(W, xCut, xCut, W->yMin, W->yMax, W->zMin, W->zMax)
void deleteSlice(slice);
void (h5)(world, slice, float, char*, char*, colormap);
void (txt)(world, slice, float, char*, char*, colormap);
void (png)(world, slice, float, char*, char*, colormap);
void sliceSnap(world, field, slice, char*, char*, ...);
#define sliceSnap(W,F,S,...) ((sliceSnap)(W, F, S, #F, #S, __VA_ARGS__))
void sliceTimeAvg(world, field, slice, char*, char*, ...);
#define sliceTimeAvg(W,F,S,...) ((sliceTimeAvg)(W, F, S, #F, #S, __VA_ARGS__))
void sliceFreqDom(world, field, slice, char*, char*, int, float, ...);
#define sliceFreqDom(W,F,S,...) ((sliceFreqDom)(W, F, S, #F, #S, __VA_ARGS__))
float sliceMax(world, field, slice);
float sliceSum(world, field, slice);
float worldMax(world, field);
float worldSum(world, field);
float get(world, field, float, float, float);
float objectAbsorption(world, object);
float poyntingX(world W, float, float, float, float, float, ...);
float poyntingY(world W, float, float, float, float, float, ...);
float poyntingZ(world W, float, float, float, float, float, ...);
float poyntingOut(world, float, float, float, float, float, float);
#define poyntingX(...) ((poyntingX)(__VA_ARGS__,-INF,INF,-INF,INF))
#define poyntingY(...) ((poyntingY)(__VA_ARGS__,-INF,INF,-INF,INF))
#define poyntingZ(...) ((poyntingZ)(__VA_ARGS__,-INF,INF,-INF,INF))
#define poyntingIn(...) (-poyntingOut(__VA_ARGS__))
field RI, LogRI, ContourRI, Contour;
field raweEffx, raweEffy, raweEffz, eEff;
field rawEx, rawEy, rawEz, rawHx, rawHy, rawHz, rawJx, rawJy, rawJz;
field Ex, Ey, Ez, Hx, Hy, Hz, iEx, iEy, iEz, iHx, iHy, iHz, Jx, Jy, Jz;
field ExHy, ExHz, EyHz, EyHx, EzHx, EzHy, Sx, Sy, Sz;
field JE, EE, UE, HH, UH, U, LogJE, LogEE, LogUE, LogHH, LogUH, LogU;
field ScEx, ScEy, ScEz, ScEE, SoEx, SoEy, SoEz, SoEE;
field DivE, JxHy, JxHz, JyHz, JyHx, JzHx, JzHy, Fx, Fy, Fz;
field Txx, Txy, Txz, Tyx, Tyy, Tyz, Tzx, Tzy, Tzz;


/* phaser.c and farfld.c */

typedef void projmap(int, int, int, int, float[2]);

typedef struct {
	float A;
	float complex *Ex, *Ey, *Ez, *Hx, *Hy, *Hz;
} tfield;

typedef struct {
	float complex ***Ex, ***Ey, ***Ez, ***Hx, ***Hy, ***Hz;
} cfield;

struct phaser {
	world W;
	int iMin, iMax, jMin, jMax, kMin, kMax, kMIN, xNUM, yNUM;
	unsigned long int N;
	float f, A, *dE, *dH;
	float complex expikx, expiky, *ex, *ez;
	cfield xMin, xMax, yMin, yMax, zMin, zMax;
};

phaser createPhaser(world, float, float, float, float, float, float, float, ...);
#define createPhaser(...) ((createPhaser)(__VA_ARGS__,-INF,INF,-INF,INF,-INF,INF))
phaser createPhasers(world, float, float, float, float, float, float, float, float, float, ...);
#define createPhasers(...) ((createPhasers)(__VA_ARGS__,-INF,INF,-INF,INF,-INF,INF))
phaser choosePhaser(phaser, float);
void updatePhaser(world, phaser);
#define updatePhasers(...) updatePhaser(__VA_ARGS__)
void deletePhaser(phaser);
#define deletePhasers(...) deletePhaser(__VA_ARGS__)
tfield Ppol(phaser, float, float);
tfield Spol(phaser, float, float);
tfield Unpol(phaser, float, float);
float phaserPoynting(phaser);
float phaserPoyntingIn(phaser);
float phaserPoyntingOut(phaser);
float complex farField(phaser, tfield(phaser,float,float), float, float);
float farFieldI(phaser, tfield(phaser,float,float), float, float);
float farFieldFlux(phaser, tfield(phaser,float,float), float, float, float, float, float, ...);
#define farFieldFlux(...) ((farFieldFlux)(__VA_ARGS__,0,0))
void farFieldTheta(world, phaser, float, char*, ...);
void farFieldPhi(world, phaser, float, char*, ...);
void farFieldProfile(world, phaser, projmap, int, int, ...);
void writePoyntingSpectrum(world, char*, ...);
#define writePoyntingSpectrum(...) ((writePoyntingSpectrum)(__VA_ARGS__,0))
void writeFarFieldSpectrum(world, tfield(phaser,float,float), float, float, char*, ...);
#define writeFarFieldSpectrum(...) ((writeFarFieldSpectrum)(__VA_ARGS__,0))
projmap CylindericalMap, AzimuthalMap, NorthernAzimuthalMap, SouthernAzimuthalMap;


/* source.c */

typedef float waveform(world, source, float, float);
struct source {
	waveform *waveForm;
	float wavelength, supplement, resolution, amplitude, phase, frequency, deviation, peaktime;
	float ***Amplitude, ***Phase, ***F;
	int iMin, iMax, jMin, jMax, kMin, kMax, nMin, nMax;
	unsigned long int N;
	source S;
};
struct tfsf {
	int iMin, iMax, jMin, jMax, kMin, kMax;
	int iMIN, iMAX, jMIN, jMAX, kMIN, kMAX;
	world W;
	tfsf S;
};

waveform Sine, Pulse, Band, WhiteBand, BrownBand, PinkBand;
tfsf createTFSF(world, dom, sur);
void pointDipole(world, field, float, float, float, waveform, float, float, float, float, ...);
void guidedWaveX(world W, char *file, field F, float x, waveform waveForm, float wavelength, float supplement, float amplitude, float phase, ...);
void guidedWaveY(world W, char *file, field F, float x, waveform waveForm, float wavelength, float supplement, float amplitude, float phase, ...);
void guidedWaveZ(world W, char *file, field F, float x, waveform waveForm, float wavelength, float supplement, float amplitude, float phase, ...);
void planewave(world, dom, tfield(phaser,float,float), float, float, waveform, float, float, float, float, ...);
void focusedBeam(world, dom, field, float, float, float, float, waveform, float, float, float, float, float, ...);
void removeSources(world W);
#define pointDipole(...) ((pointDipole)(__VA_ARGS__,0,0,0,0,0))
#define guidedWaveX(...) ((guidedWaveX)(__VA_ARGS__,0,0,0,0,0))
#define guidedWaveY(...) ((guidedWaveY)(__VA_ARGS__,0,0,0,0,0))
#define guidedWaveZ(...) ((guidedWaveZ)(__VA_ARGS__,0,0,0,0,0))
#define planewave(...) ((planewave)(__VA_ARGS__,0,0,0,0,0))
#define focusedBeam(...) ((focusedBeam)(__VA_ARGS__,0,0,0,0,0,0))
