#include <alis.h>


int main(int argc, char **argv)
{
	dom DTF = {{500}, {500}, {300}}; // Domain of total-field area
	dom DSF = {{500}, {500}, {500}}; // Domain of scattered-field area
	res Res = {10, {20}, {1240}};
	sur Sur = {{BBC}, {SYM}, {PML}};
	world W = createWorld(DSF, Res, Sur, "%s", argv[0]);

	object Nanoblock = {Lattice, {{500}, {500}}, objects {{Box, {{0, 200}, {-100, 100}, {-100, 100}}}}};
	object Substrate = {Box, {{-INF, INF}, {-INF, INF}, {-INF,-100}}};
	putObjects(W, n(4), Nanoblock, n(1.5), Substrate, Air);
	planewave(W, DTF, Ppol, 55, 0, Sine, 530, 5);

	float abs=0, sca=0;
	slice XY = createSliceXY(W, 0);
	slice XZ = createSliceXZ(W, 0);
	slice YZ = createSliceYZ(W, 0);
	phaser P = createPhaser(W, 500,-2000, 2000,-2000, 2000);

	for (int n=1, N=10*W->T; timer(n, N); n++) {
		updateE(W);
		updateH(W);

		if (N-n < W->T) {
			updatePhaser(W, P);
			sliceSnap(W, EE, XY, 10, png(hot,2), "/%%/");
			sliceSnap(W, EE, XZ, 10, png(hot,2), "/%%/");
			sliceSnap(W, EE, YZ, 10, png(hot,2), "/%%/");
		}
	}

	farFieldProfile(W, P, AzimuthalMap, 360, 180, png(jet,0), "/Log-");
	// AzimuthalMap, SouthernAzimuthalMap, NorthernAzimuthalMap or CylindericalMap
}
