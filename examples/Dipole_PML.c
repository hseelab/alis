#include <alis.h>


int main(int argc, char **argv)
{
	dom Dom = {{1000}, {1000}, {1000}};
	res Res = {10, {20}, {1240}};
	sur Sur = {{SYM, PML}, {SYM, PML}, {SYM, PML}, {24}};
	world W = createWorld(Dom, Res, Sur, "%s", argv[0]);

//	putObjects(W, Air);

	pointDipole(W, Ez, 0, 0, 0, Pulse, 400, 800, 10000000);

	slice XY = createSlice(W, W->xMIN, W->xMAX, W->yMIN, W->yMAX, 0, 0);
	slice XZ = createSlice(W, W->xMIN, W->xMAX, 0, 0, W->zMIN, W->zMAX);
	slice YZ = createSlice(W, 0, 0, W->yMIN, W->yMAX, W->zMIN, W->zMAX);

	for (int n=1, N=3000/W->dt; timer(n, N); n++) {
		updateE(W);
		updateH(W);
		sliceSnap(W, Ez, XY, 10, png(dkbr,1), "/%%/");
		sliceSnap(W, Ez, XZ, 10, png(dkbr,1), "/%%/");
		sliceSnap(W, Ez, YZ, 10, png(dkbr,1), "/%%/");
	}
}
