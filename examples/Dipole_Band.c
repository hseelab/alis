#include <alis.h>


int main(int argc, char **argv)
{
	dom Dom = {{500}, {500}, {500}};
	res Res = {10, {20}, {1240}, {{2, 100}, {2, 100}, {2, 100}}};
	sur Sur = {{SYM, PML}, {SYM, PML}, {SYM, PML}};
	world W = createWorld(Dom, Res, Sur, "%s", argv[0]);

	pointDipole(W, Ex, 0, 0, 0, Band, 400, 800);
	phaser P = createPhasers(W, 400, 800, 5);

	for (int n=1; timer(n, P->N); n++) {
		updateE(W);
		updateH(W);
		updatePhaser(W, P);
		writeSpectrum(W, P->N, 400, 800, "/Near-field", get(W,Ex,0,0,250)*250*sqrtf(8*PI/3), get(W,Ex,0,0,500)*500*sqrtf(8*PI/3));
	}

	writePoyntingSpectrum(W, "/Power", P);
	writeFarFieldSpectrum(W, Unpol, 0, 0, "/Far-field", P);
	farFieldTheta(W, choosePhaser(P,500), 0, "/Far-field-theta");
	farFieldProfile(W, choosePhaser(P,500), AzimuthalMap, 360, 180, png(hot,0), "/");

	printf("Total radiated power: %f\n", phaserPoyntingOut(choosePhaser(P,500)));
	printf("Intensity to (0,0) / (3/(8*PI)): %f\n", farFieldI(choosePhaser(P,500),Unpol,0,0) * 8*PI/3);
	printf("Radiated power to hemisphere * 2: %f\n", farFieldFlux(choosePhaser(P,500),Unpol,0,0,90) * 2);
}
