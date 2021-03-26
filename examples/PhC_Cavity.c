#include <alis.h>


int main(int argc, char **argv)
{
	dom Dom = {{2250}, {2250}, {1000}};
	res Res = {10, {20}};
	sur Sur = {{SYM, PML}, {SYM, PML}, {SYM, PML}};
	world W = createWorld(Dom, Res, Sur, "%s", argv[0]);

	object Slab = {Difference, {2}, objects {
						{Box, {{-INF, INF}, {-INF, INF}, {-100, 100}}},
						{Difference, {2}, objects {
							{DLattice, {{500, 9}, {500*sqrtf(3), 5}}, objects {
								{RodZ, {{-175, 175}, {-175, 175}, {-INF, INF}}},
							}},
							{RodZ, {{-175, 175}, {-175, 175}, {-INF, INF}}},
						}},
					}};
	putObjects(W, n(3.4), Slab, Air);
	pointDipole(W, Hz, 200, 0, 0, Band, 750, 2000);

	slice XY = createSliceXY(W, 0);
	slice XZ = createSliceXZ(W, 0);
	slice YZ = createSliceYZ(W, 0);
	sliceSnap(W, LogRI, XY, h5, "/%%");

	for (int n=1, N=750*2000/10/W->dt; timer(n, W->N+N); n++) {
		updateE(W);
		updateH(W);

		if (n > W->N) {
			writeRow(W, "/Time", W->t, get(W, Hz, 200, 200, 0));
			writeSpectrum(W, N, 750, 2000, "/Spectrum", get(W, Hz, 200, 200, 0));
		}
		if (n+2*W->T > W->N+N) {
			sliceSnap(W, EE, XY, 15, png(hot,-1), "/%%/");
			sliceTimeAvg(W, HH, XY, 2*W->T, png(jet,0), "/%%/");
			sliceFreqDom(W, Hz, XY, 2*W->T, 1500, png(dkbr,0), "/%%/");
		}
	}
}
