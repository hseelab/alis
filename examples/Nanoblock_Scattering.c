#include <alis.h>


int main(int argc, char **argv)
{
	dom DTF = {{300}, {300}, {300}}; // Domain of total-field area
	dom DSF = {{400}, {400}, {400}}; // Domain of scattered-field area
	res Res = {5, {10}, {1240}};
	sur Sur = {{PML}, {SYM, PML}, {PML}};
	world W = createWorld(DSF, Res, Sur, "%s", argv[0]);

	object Nanoblock = {Box, {{-200, 200}, {-200, 200}, {-200, 200}}};
	object Substrate = {Box, {{-INF, INF}, {-INF, INF}, {-INF,-200}}};
	putObjects(W, Au, Nanoblock, n(1.5), Substrate, Air);
	planewave(W, DTF, Ppol, 55, 0, Sine, 500, 20);

	float abs=0, sca=0;
	slice XY = createSliceXY(W, 0);
	slice XZ = createSliceXZ(W, 0);
	slice YZ = createSliceYZ(W, 0);

	for (int n=1, N=40*W->T; timer(n, N); n++) {
		updateH(W);
		abs += poyntingIn (W,-250, 250,-250, 250,-250, 250);
		sca += poyntingOut(W,-350, 350,-350, 350,-350, 350);
		updateE(W);
		abs += poyntingIn (W,-250, 250,-250, 250,-250, 250);
		sca += poyntingOut(W,-350, 350,-350, 350,-350, 350);

		if (n%W->T == 0) {
			writeRow(W, "/Cross-section", W->t, abs/W->T, sca/W->T, (abs+sca)/W->T);
			abs = sca = 0;
		}
		if (N-n < W->T) {
			sliceSnap(W, ScEE, XY, 10, png(hot,2), "/%%/");
			sliceSnap(W, ScEE, XZ, 10, png(hot,2), "/%%/");
			sliceSnap(W, ScEE, YZ, 10, png(hot,2), "/%%/");
		}
	}
}
