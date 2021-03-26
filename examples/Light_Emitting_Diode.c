#include <alis.h>


int main(int argc, char **argv)
{
	dom Dom = {{500}, {500}, {-200, 1500}};
	res Res = {10, {20}, {1240}};
	sur Sur = {{SYM, PBC}, {SYM, PBC}, {PEC, PML}};
	world W = createWorld(Dom, Res, Sur, "%s", argv[0]);

	matter InGaN  = nk(2.46, 0.045, 1240/450);
	object MQWU   = {Box, {{-INF, INF}, {-INF, INF}, {550, 600}}};
	object MQWL   = {Box, {{-INF, INF}, {-INF, INF}, {400, 450}}};
	object Slab   = {Box, {{-INF, INF}, {-INF, INF}, {0, 1000}}};
	object Mirror = {Box, {{-INF, INF}, {-INF, INF}, {-INF, 0}}};
	putObjects(W, InGaN, MQWU, InGaN, MQWL, n(2.46), Slab, Ag, Mirror, n(1.4));

	pointDipole(W, Ex, 0, 0, 500, Pulse, 450, 25);

	float inc=0, ext=0, abs=0;
	slice XY = createSliceXY(W, 0);
	slice XZ = createSliceXZ(W, 0);
	slice YZ = createSliceYZ(W, 0);

	for (int n=1, N=500000/W->dt; timer(n, N); n++) {
		updateE(W);
		inc -= 2 * W->dx * W->dy * W->dz * get(W, JE, 0, 0, 500);
		ext += poyntingZ(W, W->zMax);
		abs += objectAbsorption(W, Mirror);
		abs += objectAbsorption(W, MQWU);
		abs += objectAbsorption(W, MQWL);

		updateH(W);
		inc -= 2 * W->dx * W->dy * W->dz * get(W, JE, 0, 0, 500);
		ext += poyntingZ(W, W->zMax);
		abs += objectAbsorption(W, Mirror);
		abs += objectAbsorption(W, MQWU);
		abs += objectAbsorption(W, MQWL);

		if (n%W->T==0) writeRow(W, "/Absolute", W->dt*n, ext, abs, inc);
		if (n%W->T==0) writeRow(W, "/Relative", W->dt*n, ext/inc, abs/inc, (ext+abs)/inc);
	}
}
