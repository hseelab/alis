/**********************************************************************
 * ALiS is Copyright (C) 2009-2019 by Ho-Seok Ee <hsee@kongju.ac.kr>. *
 * Redistribution and use with or without modification, are permitted *
 * under the terms of the Artistic License version 2.                 *
 **********************************************************************/

#include "alis_c.h"

#define	MINMAX xMin,xMax,yMin,yMax,zMin,zMax;\
	if (O.f[2][1]-O.f[2][0]) xMin=O.f[0][0],xMax=O.f[0][1],yMin=O.f[1][0],yMax=O.f[1][1],zMin=O.f[2][0],zMax=O.f[2][1];\
	else xMin=O.f[0][0]-O.f[1][0],yMin=O.f[0][1]-O.f[1][1]?O.f[1][1]:O.f[1][0],zMin=O.f[0][2]-O.f[1][2]?O.f[1][2]:O.f[1][0],\
		 xMax=O.f[0][0]+O.f[1][0],yMax=O.f[0][1]+O.f[1][1]?O.f[1][1]:O.f[1][0],zMax=O.f[0][2]+O.f[1][2]?O.f[1][2]:O.f[1][0]


/* Basic objects */

int Box(object O, float x, float y, float z) {
	float MINMAX;
	if (xMin <= x && x <= xMax && yMin <= y && y <= yMax && zMin <= z && z <= zMax) return 1;
	return 0;
}

int Ball(object O, float x, float y, float z) {
	float MINMAX;
	if (sq((2*x-xMin-xMax)/(xMax-xMin))+sq((2*y-yMin-yMax)/(yMax-yMin))+sq((2*z-zMin-zMax)/(zMax-zMin)) < 1) return 1;
	return 0;
}

int Diamond(object O, float x, float y, float z) {
	float MINMAX;
	if (fabsf((2*x-xMin-xMax)/(xMax-xMin))+fabsf((2*y-yMin-yMax)/(yMax-yMin))+fabsf((2*z-zMin-zMax)/(zMax-zMin)) <= 1) return 1;
	return 0;
}

int BiconeX(object O, float x, float y, float z) {
	float MINMAX;
	if (fabsf((2*x-xMin-xMax)/(xMax-xMin))+sqrt(sq((2*y-yMin-yMax)/(yMax-yMin))+sq((2*z-zMin-zMax)/(zMax-zMin))) < 1) return 1;
	return 0;
}

int BiconeY(object O, float x, float y, float z) {
	float MINMAX;
	if (fabsf((2*y-yMin-yMax)/(yMax-yMin))+sqrt(sq((2*z-zMin-zMax)/(zMax-zMin))+sq((2*x-xMin-xMax)/(xMax-xMin))) < 1) return 1;
	return 0;
}

int BiconeZ(object O, float x, float y, float z) {
	float MINMAX;
	if (fabsf((2*z-zMin-zMax)/(zMax-zMin))+sqrt(sq((2*x-xMin-xMax)/(xMax-xMin))+sq((2*y-yMin-yMax)/(yMax-yMin))) < 1) return 1;
	return 0;
}

int RodX(object O, float x, float y, float z) {
	float MINMAX;
	if (xMin <= x && x <= xMax && sq((2*y-yMin-yMax)/(yMax-yMin))+sq((2*z-zMin-zMax)/(zMax-zMin)) < 1) return 1;
	return 0;
}

int RodY(object O, float x, float y, float z) {
	float MINMAX;
	if (yMin <= y && y <= yMax && sq((2*x-xMin-xMax)/(xMax-xMin))+sq((2*z-zMin-zMax)/(zMax-zMin)) < 1) return 1;
	return 0;
}

int RodZ(object O, float x, float y, float z) {
	float MINMAX;
	if (zMin <= z && z <= zMax && sq((2*x-xMin-xMax)/(xMax-xMin))+sq((2*y-yMin-yMax)/(yMax-yMin)) < 1) return 1;
	return 0;
}

int DRodX(object O, float x, float y, float z) {
	float MINMAX;
	if (xMin <= x && x <= xMax && fabsf((2*y-yMin-yMax)/(yMax-yMin))+fabsf((2*z-zMin-zMax)/(zMax-zMin)) <= 1) return 1;
	return 0;
}

int DRodY(object O, float x, float y, float z) {
	float MINMAX;
	if (yMin <= y && y <= yMax && fabsf((2*z-zMin-zMax)/(zMax-zMin))+fabsf((2*x-xMin-xMax)/(xMax-xMin)) <= 1) return 1;
	return 0;
}

int DRodZ(object O, float x, float y, float z) {
	float MINMAX;
	if (zMin <= z && z <= zMax && fabsf((2*x-xMin-xMax)/(xMax-xMin))+fabsf((2*y-yMin-yMax)/(yMax-yMin)) <= 1) return 1;
	return 0;
}

int HRodXY(object O, float x, float y, float z) {
	float MINMAX;
	if (xMin <= x && x <= xMax && zMin < z && z < zMax && fabsf((2*z-zMin-zMax)/(zMax-zMin))+fabsf((2*y-yMin-yMax)/(yMax-yMin))*2 < 2) return 1;
	return 0;
}

int HRodYZ(object O, float x, float y, float z) {
	float MINMAX;
	if (yMin <= y && y <= yMax && xMin < x && x < xMax && fabsf((2*x-xMin-xMax)/(xMax-xMin))+fabsf((2*z-zMin-zMax)/(zMax-zMin))*2 < 2) return 1;
	return 0;
}

int HRodZX(object O, float x, float y, float z) {
	float MINMAX;
	if (zMin <= z && z <= zMax && yMin < y && y < yMax && fabsf((2*y-yMin-yMax)/(yMax-yMin))+fabsf((2*x-xMin-xMax)/(xMax-xMin))*2 < 2) return 1;
	return 0;
}

int HRodXZ(object O, float x, float y, float z) {
	float MINMAX;
	if (xMin <= x && x <= xMax && yMin < y && y < yMax && fabsf((2*y-yMin-yMax)/(yMax-yMin))+fabsf((2*z-zMin-zMax)/(zMax-zMin))*2 < 2) return 1;
	return 0;
}

int HRodYX(object O, float x, float y, float z) {
	float MINMAX;
	if (yMin <= y && y <= yMax && zMin < z && z < zMax && fabsf((2*z-zMin-zMax)/(zMax-zMin))+fabsf((2*x-xMin-xMax)/(xMax-xMin))*2 < 2) return 1;
	return 0;
}

int HRodZY(object O, float x, float y, float z) {
	float MINMAX;
	if (zMin <= z && z <= zMax && xMin < x && x < xMax && fabsf((2*x-xMin-xMax)/(xMax-xMin))+fabsf((2*y-yMin-yMax)/(yMax-yMin))*2 < 2) return 1;
	return 0;
}


/* Compound objects */

int Translation(object O, float x, float y, float z) {
	object *o = O.p;
	return o->S(*o, x-O.f[0][0], y-O.f[0][1], z-O.f[0][2]);
}

int Rotation(object O, float x, float y, float z) {
	object *o = O.p;
	x -= O.f[1][0], y -= O.f[1][1], z -= O.f[1][2];
	float cx = cosf(O.f[0][0]*PI/180.), sx = sinf(O.f[0][0]*PI/180.);
	float cy = cosf(O.f[0][1]*PI/180.), sy = sinf(O.f[0][1]*PI/180.);
	float cz = cosf(O.f[0][2]*PI/180.), sz = sinf(O.f[0][2]*PI/180.);
	return o->S(*o, cy*(cz*x+sz*y)-sy*z, cx*(cz*y-sz*x)+sx*(cy*z+sy*(cz*x+sz*y)), cx*(cy*z+sy*(cz*x+sz*y))-sx*(cz*y-sz*x));
}

int Combination(object O, float x, float y, float z) {
	object *o = O.p;
	for (int n=0, N=O.f[0][0]; n<N; n++) if (o[n].S(o[n], x, y, z)) return 1;
	return 0;
}

int Intersection(object O, float x, float y, float z) {
	object *o = O.p;
	for (int n=0, N=O.f[0][0]; n<N; n++) if (!o[n].S(o[n], x, y, z)) return 0;
	return 1;
}

int Difference(object O, float x, float y, float z) {
	object *o = O.p;
	if (!o[0].S(o[0], x, y, z)) return 0;
	for (int n=1, N=O.f[0][0]; n<N; n++) if (o[n].S(o[n], x, y, z)) return 0;
	return 1;
}

int Lattice(object O, float x, float y, float z) {
	object *o = O.p;
	float i=0, j=0, k=0;
	if (O.f[0][0]) {
		if (O.f[0][1]) x += 0.5 * (O.f[0][1]-1) * O.f[0][0];
		i = floor((x + 0.5 * (O.f[0][0] - o->f[0][0] - o->f[0][1])) / O.f[0][0]);
		x -= O.f[0][0] * (O.f[0][1] > 0 ? MAX(0, MIN(i, O.f[0][1]-1)) : i);
	}
	if (O.f[1][0]) {
		if (O.f[1][1]) y += 0.5 * (O.f[1][1]-1) * O.f[1][0];
		j = floor((y + 0.5 * (O.f[1][0] - o->f[1][0] - o->f[1][1])) / O.f[1][0]);
		y -= O.f[1][0] * (O.f[1][1] > 0 ? MAX(0, MIN(j, O.f[1][1]-1)) : j);
	}
	if (O.f[2][0]) {
		if (O.f[2][1]) z += 0.5 * (O.f[2][1]-1) * O.f[2][0];
		k = floor((z + 0.5 * (O.f[2][0] - o->f[2][0] - o->f[2][1])) / O.f[2][0]);
		z -= O.f[2][0] * (O.f[2][1] > 0 ? MAX(0, MIN(k, O.f[2][1]-1)) : k);
	}
	if (o->S(*o, x, y, z)) return 1;
	return 0;
}

int DLattice(object O, float x, float y, float z) {
	object *o = O.p;
	if (Lattice(O, x, y, z)) return 1;
	if (O.f[0][0]) O.f[0][1] -= O.f[0][1] ? 1 : 2;
	if (O.f[1][0]) O.f[1][1] -= O.f[1][1] ? 1 : 2;
	if (O.f[2][0]) O.f[2][1] -= O.f[2][1] ? 1 : 2;
	if (Lattice(O, x, y, z)) return 1;
	return 0;
}

int CLatticeX(object O, float x, float y, float z) {
	object *o = O.p;
	float y0 = O.f[1][0], z0 = O.f[1][1];
	for (int n=0, N=O.f[0][0]; n<N; n++) {
		float theta = 2 * PI * n / N;
		if (o->S(*o, x, (y-y0)*cosf(theta) + (z-z0)*sinf(theta) + y0, (z-z0)*cosf(theta) - (y-y0)*sinf(theta) + z0)) return 1;
	}
	return 0;
}

int CLatticeY(object O, float x, float y, float z) {
	object *o = O.p;
	float x0 = O.f[1][0], z0 = O.f[1][1];
	for (int n=0, N=O.f[0][0]; n<N; n++) {
		float theta = 2 * PI * n / N;
		if (o->S(*o, (x-x0)*cosf(theta) + (z-z0)*sinf(theta) + x0, y, (z-z0)*cosf(theta) - (x-x0)*sinf(theta) + z0)) return 1;
	}
	return 0;
}

int CLatticeZ(object O, float x, float y, float z) {
	object *o = O.p;
	float x0 = O.f[1][0], y0 = O.f[1][1];
	for (int n=0, N=O.f[0][0]; n<N; n++) {
		float theta = 2 * PI * n / N;
		if (o->S(*o, (x-x0)*cosf(theta) + (y-y0)*sinf(theta) + x0, (y-y0)*cosf(theta) - (x-x0)*sinf(theta) + y0, z)) return 1;
	}
	return 0;
}


/* Etc. */

int ImportContour(object O, float x, float y, float z)
{
	float **f = O.p;
	if (O.f[0][0] < x && x < O.f[0][1] && O.f[1][0] < y && y < O.f[1][1]) {
		int i = O.f[0][2] * (x - O.f[0][0]) / (O.f[0][1] - O.f[0][0]);
		int j = O.f[1][2] * (y - O.f[1][0]) / (O.f[1][1] - O.f[1][0]);
		if (O.f[2][0] <= f[i][j] && f[i][j] <= O.f[2][1]) return 1;
	}
	return 0;
}
