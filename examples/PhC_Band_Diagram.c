#include "alis.h"


int main(int argc, char **argv)
{
	float kx, ky, k=1.5;

	if (k<=0 || k>=3) kx = 0, ky = 0;
	else if (k < 1) kx =   0, ky =   k; // G - M
	else if (k < 2) kx = k-1, ky =   1; // M - K
	else if (k < 3) kx = 3-k, ky = 3-k; // K - G

	dom Dom = {{250}, {0, 250*sqrtf(3)}, {500}};
	res Res = {5, {10, 5*sqrtf(3), 10}, {1240}};
	sur Sur = {{BBC, 500*3/kx}, {BBC, 500*sqrtf(3)/ky}, {SYM, PML}};
	world W = createWorld(Dom, Res, Sur, "%s-%03.1f", argv[0], k);
	object Slab = {Difference, {2}, objects {
						{Box, {{-INF, INF}, {-INF, INF}, {-100, 100}}},
						{DLattice, {{500}, {500*sqrtf(3)}}, objects {
							{RodZ, {{-200, 200}, {-200, 200}, {-INF, INF}}},
						}},
					}};
	putObjects(W, n(3.5), Slab, Air);
	pointDipole(W, iHz, 200, 100*sqrtf(3), 0, Band, 750, 2000);

	slice XY = createSliceXY(W, 0);
	sliceSnap(W, ContourRI, XY, png(gray,0), "/%%");
	
	for (int n=1, N=750*2000/5/W->dt; timer(n, W->N+N); n++) {
		updateE(W);
		updateH(W);
		copyBoundary(W, H, x, cosf(PI*(ky+kx/3)), sinf(PI*(ky+kx/3)), W->iMIN, 0, W->jMIN,   W->jMIN,   W->kMIN, W->kMAX, W->iMAX, W->jNUM-1, 0);
		copyBoundary(W, H, z, cosf(PI*(ky+kx/3)), sinf(PI*(ky+kx/3)), W->iMIN, 0, W->jMIN,   W->jMIN,   W->kMIN, W->kMAX, W->iMAX, W->jNUM-1, 0);
		copyBoundary(W, H, x, cosf(PI*(ky-kx/3)),-sinf(PI*(ky-kx/3)), W->iMIN, 0, W->jMAX+1, W->jMAX+1, W->kMIN, W->kMAX, W->iMAX, 1-W->jNUM, 0);
		copyBoundary(W, H, z, cosf(PI*(ky-kx/3)),-sinf(PI*(ky-kx/3)), W->iMIN, 0, W->jMAX+1, W->jMAX+1, W->kMIN, W->kMAX, W->iMAX, 1-W->jNUM, 0);
		copyBoundary(W, H, x, cosf(PI*(ky-kx/3)), sinf(PI*(ky-kx/3)), 1, W->iMAX+1, W->jMIN,   W->jMIN,   W->kMIN, W->kMAX, W->iMIN, W->jNUM-1, 0);
		copyBoundary(W, H, z, cosf(PI*(ky-kx/3)), sinf(PI*(ky-kx/3)), 1, W->iMAX+1, W->jMIN,   W->jMIN,   W->kMIN, W->kMAX, W->iMIN, W->jNUM-1, 0);
		copyBoundary(W, H, x, cosf(PI*(ky+kx/3)),-sinf(PI*(ky+kx/3)), 1, W->iMAX+1, W->jMAX+1, W->jMAX+1, W->kMIN, W->kMAX, W->iMIN, 1-W->jNUM, 0);
		copyBoundary(W, H, z, cosf(PI*(ky+kx/3)),-sinf(PI*(ky+kx/3)), 1, W->iMAX+1, W->jMAX+1, W->jMAX+1, W->kMIN, W->kMAX, W->iMIN, 1-W->jNUM, 0);

		if (n > W->N) {
			sliceFreqDom(W, Hz, XY, N, 1596, png(dkbr,0), "/%%/");
			writeSpectrum(W, N, 750, 2000, "", get(W, Hz, 200, 0, 0));
		}
	}
}
