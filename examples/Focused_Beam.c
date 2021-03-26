#include <alis.h>


int main(int argc, char **argv)
{
	dom Dom = {{3000}, {3000}, {3000}};
	res Res = {20, {40}, {1240}};
	sur Sur = {{SYM, PML}, {SYM, PML}, {PML}};
	world W = createWorld(Dom, Res, Sur, "%s", argv[0]);

	focusedBeam(W, Dom, Ex, 0, INF, 0, 0, Sine, 500, 3, 2);

	slice XY = createSliceXY(W, 0);
	slice XZ = createSliceXZ(W, 0);
	slice YZ = createSliceYZ(W, 0);

	float inten = 0, power = 0;

	for (int n=1, N=10000/W->dt; timer(n, N); n++) {
		updateE(W);
		inten -= get(W,Sz,0,0,0);
		power -= poyntingZ(W, 0);
		updateH(W);
		inten -= get(W,Sz,0,0,0);
		power -= poyntingZ(W, 0);

		if (n%W->T == 0) writeRow(W, "/Power", W->t, inten/W->T, power/(W->T*sq(500*2))), inten=power=0;
		if (N-n < W->T) {
			sliceSnap(W, Ex, XY, 10, png(dkbr,1), "/%%/");
			sliceSnap(W, Ex, XZ, 10, png(dkbr,1), "/%%/");
			sliceSnap(W, Ex, YZ, 10, png(dkbr,1), "/%%/");
		}
	}
}
