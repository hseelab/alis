#include <alis.h>


int main(int argc, char **argv)
{
	dom DTF = {{180}, {0}, {180}};
	dom DSF = {{200}, {0}, {200}};
	res Res = {2.5, {5}, {1240}};
	sur Sur = {{SYM, PML}, {SYM}, {PML, PML}};
	world W = createWorld(DSF, Res, Sur, "%s", argv[0]);

	object Nanowire  = {RodY, {{ -40,  40}, {-INF, INF}, { -40,  40}}};
	object Substrate = {Box,  {{-INF, INF}, {-INF, INF}, {-INF, -40}}};
	putObjects(W, Si, Nanowire, n(1.45), Substrate, Air);

	planewave(W, DTF, Spol, 0, 0, Band, 400, 1000);
	slice XZ = createSliceXZ(W, 0);
	phaser A = createPhasers(W, 400, 1000, 10, -170, 170, -INF, INF, -170, 170);
	phaser S = createPhasers(W, 400, 1000, 10, -190, 190, -INF, INF, -190, 190);

	for (int n=1; timer(n, S->N); n++) {
		updateE(W);
		updateH(W);
		updatePhaser(W, A);
		updatePhaser(W, S);
		sliceFreqDom(W, ScEy, XZ, S->N, 473, png(dkbr,0), "/%%-473/");
		sliceFreqDom(W, ScEy, XZ, S->N, 800, png(dkbr,0), "/%%-800/");
	}

	writePoyntingSpectrum(W, "", A, S);
}
