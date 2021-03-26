#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#define MALLOC(p,n) (p)=calloc((n),sizeof(*(p)))

typedef struct { float *dx, *dy, *dz, ***x, ***y, ***z; } vfield;

float ***makefield(int N)
{
	float ***f;
	MALLOC(f, N+1);
	MALLOC(f[0], (N+1)*(N+1));
	MALLOC(f[0][0], (N+1)*(N+1)*(N+1));
	#pragma omp parallel for
	for (int i=0; i<=N; i++) { f[i]=f[0]+(N+1)*(i-0); for (int j=0; j<=N; j++) f[i][j]=f[0][0]+(N+1)*((N+1)*(i-0)+(j-0)); }
	return f;
}

static void updateF(vfield F, vfield G, int m, int n, int N)
{
	#pragma omp parallel for
	for (int i=n; i<n+N; i++) for (int j=n; j<n+N; j++) {
		#pragma GCC ivdep
		for (int k=n; k<n+N; k++) {
			F.x[i][j][k] += F.dz[k] * (G.y[i][j][k] - G.y[i][j][k+m]) + F.dy[j] * (G.z[i][j+m][k] - G.z[i][j][k]);
			F.y[i][j][k] += F.dx[i] * (G.z[i][j][k] - G.z[i+m][j][k]) + F.dz[k] * (G.x[i][j][k+m] - G.x[i][j][k]);
			F.z[i][j][k] += F.dy[j] * (G.x[i][j][k] - G.x[i][j+m][k]) + F.dx[i] * (G.y[i+m][j][k] - G.y[i][j][k]);
		}
	}
}

int main(int argc, char **argv)
{
	vfield E, H;
	double T = (argc>1) ? atof(argv[1]) : 1;
	long n = 0, N = (argc>2) ? atoi(argv[2]) : 200;
	struct timeval start, stop;
	MALLOC(E.dx,N+1); MALLOC(E.dy,N+1); MALLOC(E.dz,N+1); MALLOC(H.dx,N+1); MALLOC(H.dy,N+1); MALLOC(H.dz,N+1);
	E.x=makefield(N); E.y=makefield(N); E.z=makefield(N); H.x=makefield(N); H.y=makefield(N); H.z=makefield(N);

	gettimeofday(&start, 0);
	do {
		n++;
		updateF(E, H, 1, 0, N);
		updateF(H, E,-1, 1, N);
		gettimeofday(&stop, 0);
	} while ((stop.tv_sec - start.tv_sec) + (stop.tv_usec - start.tv_usec) * 1.0E-6 < T);
	printf("%ld\n", n);
}
