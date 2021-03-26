#include <alis.h>


int main(int argc, char **argv)
{
	dom Dom = {{500}, {500}, {-300, 500}};
	res Res = {5, {10}, {1240}};
	sur Sur = {{SYM, PML}, {PML, PML}, {PML, PML}};
	world W = createWorld(Dom, Res, Sur, "%s", argv[0]);

	object Waveguide = {Box, {{ -50,  50}, {-INF, INF}, { -50,  50}}};
	object Substrate = {Box, {{-INF, INF}, {-INF, INF}, {-INF, -50}}};
	putObjects(W, n(4), Waveguide, n(1.5), Substrate, Air);

	guidedWaveY(W, "Waveguide_Dispersion-3.0/rawEx-XZ/07", Ex, -500, Sine, 625, 10, 1, 0);
	guidedWaveY(W, "Waveguide_Dispersion-3.0/rawEz-XZ/07", Ez, -500, Sine, 625, 10, 1, 0);

	slice XY = createSliceXY(W, 0);
	slice YZ = createSliceYZ(W, 0);

	for (int n=1, N=30*W->T; timer(n, N); n++) {
		updateE(W);
		updateH(W);

		if (N-n < 2*W->T) {
			sliceSnap(W, Ez, XY, 25, png(dkbr,-1), "/%%/");
			sliceSnap(W, Ez, YZ, 25, png(dkbr,-1), "/%%/");
		}
	}
}
