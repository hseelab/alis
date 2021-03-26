#include <alis.h>


int main(int argc, char **argv)
{
	float k = atof(argv[1]);

	dom Dom = {{500}, {0}, {-300, 500}};
	res Res = {5, {10}, {1240}};
	sur Sur = {{SYM, PML}, {BBC, 1000/k}, {PML, PML}};
	world W = createWorld(Dom, Res, Sur, "%s-%03.1f", argv[0], k);

	object Waveguide = {Box, {{ -50,  50}, {-INF, INF}, { -50,  50}}};
	object Substrate = {Box, {{-INF, INF}, {-INF, INF}, {-INF, -50}}};
	putObjects(W, n(4), Waveguide, n(1.5), Substrate, Air);

	pointDipole(W, Ez, 0, 0, 0, Band, 400, 1000);

	slice XY = createSliceXY(W, 0);
	slice XZ = createSliceXZ(W, 0);
	slice YZ = createSliceYZ(W, 0);

	for (int n=1, N=400*1000/2/W->dt; timer(n, W->N+N); n++) {
		updateE(W);
		updateH(W);

		if (n > W->N) {
			sliceFreqDom(W, rawEx, XZ, N, 625, h5, "/%%/");
			sliceFreqDom(W, rawEz, XZ, N, 625, h5, "/%%/");
			sliceFreqDom(W, Ex, XZ, N, 625, png(dkbr,0), "/%%/");
			sliceFreqDom(W, Ez, XZ, N, 625, png(dkbr,0), "/%%/");
			writeSpectrum(W, N, 400, 1000, "", get(W, Ez, 0, 0, 0));
		}
	}
}
