#include <alis.h>


int main(int argc, char **argv)
{
	float g = 20;
	float r = 20;
	float d = 100;

	dom Dom = {{250}, {200}, {-150, 150}};
	res Res = {1.25, {2.5}, {1240,-1}, {{4, g/2+d}, {4, r+(d-2*r)/sqrtf(3)}, {4, 15}}};
	sur Sur = {{SYM, PML}, {SYM, PML}, {PML, PML}};
	world W = createWorld(Dom, Res, Sur, "%s", argv[0]);

	object Bowtie = {Combination, {3}, objects {
						{RodZ, {{g/2, g/2+r*2}, {-r, r}, {-15, 15}}},
						{RodZ, {{g/2-r*2+d, g/2+d}, {-r+(d-2*r)/sqrtf(3), r+(d-2*r)/sqrtf(3)}, {-15, 15}}},
						{Intersection, {2}, objects {
							{Box,   {{g/2+r/2, g/2+d}, {-(d-r/2)/sqrtf(3), (d-r/2)/sqrtf(3)}, {-15, 15}}},
							{DRodZ, {{g/2-r, g/2-r*2+d*2}, {-(d-r/2)/sqrtf(3), (d-r/2)/sqrtf(3)}, {-15, 15}}},
						}}
					}};
	object Substrate = {Box,  {{-INF, INF}, {-INF, INF}, {-INF, -15}}};
	putObjects(W, Au, Bowtie, n(1.5), Substrate, Air);

	pointDipole(W, Ex, 0, 0, 0, Band, 400, 1000);

	slice XY = createSliceXY(W, 0);
	slice XZ = createSliceXZ(W, 0);

	for (int n=1, N=400*1000/10/W->dt; timer(n, W->N+N); n++) {
		updateE(W);
		updateH(W);

		if (n > W->N) {
			writeRow(W, "/Time", W->t, get(W, Ex, 0, 0, 0));
			writeSpectrum(W, N, 400, 1000, "/Spectrum", get(W, Ex, 0, 0, 0));
		}
		if (n > W->N) {
			sliceFreqDom(W, Ex, XY, N, 700, png(dkbr,0), "/%%/");
			sliceFreqDom(W, Ex, XZ, N, 700, png(dkbr,0), "/%%/");
			sliceFreqDom(W, Ez, XY, N, 700, png(dkbr,0), "/%%/");
			sliceFreqDom(W, Ez, XZ, N, 700, png(dkbr,0), "/%%/");
		}
	}
}
