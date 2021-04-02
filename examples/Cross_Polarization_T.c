#include <alis.h>

float phaserPoyntingP(phaser P, cfield S, int k)
{
	float Sum=0;
	#pragma omp parallel for reduction(+:Sum)
	for (int i=P->iMin+1; i<=P->iMax; i++) for (int j=P->jMin+1; j<=P->jMax; j++)
		Sum -= (S.Ex[i][j][k]+S.Ex[i][j-1][k]+S.Ey[i][j][k]+S.Ey[i-1][j][k]) * conj(S.Hx[i][j][k]+S.Hx[i-1][j][k]-S.Hy[i][j][k]-S.Hy[i][j-1][k]);
	return 0.125 * Sum * P->W->dx * P->W->dy * sq(P->A);
}

float phaserPoyntingM(phaser P, cfield S, int k)
{
	float Sum=0;
	#pragma omp parallel for reduction(+:Sum)
	for (int i=P->iMin+1; i<=P->iMax; i++) for (int j=P->jMin+1; j<=P->jMax; j++)
		Sum += (S.Ex[i][j][k]+S.Ex[i][j-1][k]-S.Ey[i][j][k]-S.Ey[i-1][j][k]) * conj(S.Hx[i][j][k]+S.Hx[i-1][j][k]+S.Hy[i][j][k]+S.Hy[i][j-1][k]);
	return 0.125 * Sum * P->W->dx * P->W->dy * sq(P->A);
}


int main(int argc, char **argv)
{
	float a = atof(argv[1]);
	float x = atof(argv[2]);
	float y = atof(argv[3]);
	float z = atof(argv[4]);

	dom Dom = {{a/2}, {0, a*1.75/2}, {z/2+200}};
	res Res = {2.5, {5}, {1240, 16}};
	sur Sur = {{PBC}, {HBC}, {PML, PML}};
	world W = createWorld(Dom, Res, Sur, "%s_a%03.0f_x%03.0f_y%03.0f_z%03.0f", argv[0], a, x, y, z);

	object Nanowire  = {Translation, {{(x-floor(x*0.1)*10)/2, (y-floor(y*0.1)*10)/2}}, objects {{Combination, {3}, objects {
							{DLattice, {{a}, {a*1.75}}, objects {
								{Box,  {{x>y?(y-x)/2:-x/2, x>y?(x-y)/2:x/2}, {x>y?-y/2:(x-y)/2, x>y?y/2:(y-x)/2}, {-z/2, z/2}}}}},
							{DLattice, {{a}, {a*1.75}}, objects {{RodZ, {{-x/2, x>y?y-x/2:x/2}, {-y/2, x>y?y/2:x-y/2}, {-z/2, z/2}}}}},
							{DLattice, {{a}, {a*1.75}}, objects {{RodZ, {{x>y?x/2-y:-x/2, x/2}, {x>y?-y/2:y/2-x, y/2}, {-z/2, z/2}}}}}
						}}}};
	object Substrate = {Box, {{-INF, INF}, {-INF, INF}, {-INF, z/2}}};

	putObjects(W, Si, Nanowire, n(1.4), Substrate, Air);

	planewave(W, Dom, Ppol, 180,  0, Band, 400, 800, 1/sqrtf(2));
	planewave(W, Dom, Ppol, 180, 90, Band, 400, 800, 1/sqrtf(2));

	slice XY = createSliceXY(W, 0);
//	sliceSnap(W, LogRI, XY, png(gray,0), "");
//	exit(0);

	phaser T = createPhasers(W, 400, 800, 10, W->xMin, W->xMax, W->yMin, W->yMax, z/2+190, z/2+190);

	for (int n=1; timer(n, T->N); n++) {
		updateE(W);
		updateH(W);
		updatePhaser(W, T);
	}

	for (phaser p=T; p->f; p++) {
		writeRow(W, "", 1/p->f, phaserPoyntingP(p, p->zMax, p->kMax)/W->xyArea, phaserPoyntingM(p, p->zMax, p->kMax)/W->xyArea);
	}
}
