#include <alis.h>


int main(int argc, char **argv)
{
	dom Dom = {{500}, {500}, {-500, 2000}};
	res Res = {10, {20}, {1240, 1}};
	sur Sur = {{SYM, PML}, {PML, PML}, {PML, PML}};
	world W = createWorld(Dom, Res, Sur, "%s", argv[0]);

	object Layer1 = {Box, {{-INF, INF}, {-INF, INF}, {-80, 0}}};
	object Layer2 = {Box, {{-INF, INF}, {-INF, INF}, {-INF,-80}}};
	putObjects(W, n(2.4), Layer1, n(1.5), Layer2, Air);

	pointDipole(W, Ey, 0, 0, 1500, Sine, 500, 5, 0.5*sqrtf(3));
	pointDipole(W, Ez, 0, 0, 1500, Sine, 500, 5, 0.5);
	phaser P = createPhaser(W, 500);

	for (int n=1, N=50*W->T; timer(n, N); n++) {
		updateE(W);
		updateH(W);
		if (N-n < W->T) updatePhaser(W, P);
	}

	farFieldTheta(W, P, 90, "");
	farFieldProfile(W, P, AzimuthalMap, 360, 180, png(hot,0), "/");
	// AzimuthalMap, SouthernAzimuthalMap, NorthernAzimuthalMap or CylindericalMap
}
