#include <alis.h>


int main(int argc, char **argv)
{
	int d[2];
	float **f = h5read(d, "Contour_filename_without_[.h5]");

	dom Dom = {{1500}, {1500}, {500}};
	res Res = {10, {20}, {1240}};
	sur Sur = {{PML, PML}, {PML, PML}, {SYM, PML}};
	world W = createWorld(Dom, Res, Sur, "%s", argv[0]);

	object Slab = {Difference, {2}, objects {
						{Box, {{-INF, INF}, {-INF, INF}, {-100, 100}}},
						{Contour, {{-1400, 1400, d[0]}, {-1400, 1400, d[1]}, {0, 80}}, f}
					}};
	putObjects(W, n(3.4), Slab, Air);
	pointDipole(W, Hz, 0, 0, 0, Pulse, 1500, 500);

	slice XY = createSliceXY(W, 0);
	slice XZ = createSliceXZ(W, 0);
	slice YZ = createSliceYZ(W, 0);
	sliceSnap(W, LogRI, XY, png(gray,0), "/%%");
	sliceSnap(W, LogRI, XZ, png(gray,0), "/%%");
	sliceSnap(W, LogRI, YZ, png(gray,0), "/%%");

	for (int n=1, N=2048; timer(n, W->N+N); n++) {
		updateE(W);
		updateH(W);
		if (W->N+N-n < 2*W->T) {
			sliceSnap(W, Hz, XY, 15, png(dkbr,-1), "/%%-");
		}
	}
}
