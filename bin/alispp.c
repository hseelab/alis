/**********************************************************************
 * ALiS is Copyright (C) 2009-2019 by Ho-Seok Ee <hsee@kongju.ac.kr>. *
 * Redistribution and use with or without modification, are permitted *
 * under the terms of the Artistic License version 2.                 *
 **********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>


float getNumber(char *str)
{
	char no[16];
	for (int i = 0; *str; str++) {
		if ('0' <= *str && *str <='9' || *str=='.' && '0' <= *(str+1) && *(str+1) <='9') no[i++] = *str;
		else if (i) no[i] = '\0', i = 0;
	}
	return atof(no);
}


int isPeak(float *f, int n, int min, int max)
{
	while (min--) if (f[n-min] < f[n-min-1]) return 0;
	while (max--) if (f[n+max] < f[n+max+1]) return 0;
	return 1;
}


int main(int argc, char **argv)
{
	FILE *file;
	int n=5, w=1;
	float f=0, m=0;

	for (int i=1; i<argc; i++) {
		if (argv[i][0]=='-') {
			if (argv[i][1]=='f' && !(f=atof(argv[i]+2))) { if (f=atof(argv[i+1])) i++; else f=1; }
			if (argv[i][1]=='m' && !(m=atof(argv[i]+2))) { if (m=atof(argv[i+1])) i++; else m=0; }
			if (argv[i][1]=='n' && !(n=atoi(argv[i]+2))) { if (n=atoi(argv[i+1])) i++; else n=1; }
			if (argv[i][1]=='w' && !(w=atoi(argv[i]+2))) { if (w=atoi(argv[i+1])) i++; else w=1; }
		}

		else if (file=fopen(argv[i], "r")) {
			char str[1024];
			int N=0;
			float max=0, *l=0, *r=0, *L=0, *R=0;

			for (int j=0; fgets(str, 1024, file); N++) {
				if (j <= N) {
					j += 1024;
					l = realloc(l, j*sizeof(float));
					r = realloc(r, j*sizeof(float));
				}
				char *s = str;
				r[N] = 0;
				l[N] = atof(s);
				do r[N] += atof(s=1+strchr(s,'\t'));
				while (strchr(s,'\t'));
				if (max < r[N]) max = r[N];
			}

			l = realloc(l, N*sizeof(float));
			r = realloc(r, N*sizeof(float));
			L = calloc(n, sizeof(float));
			R = calloc(n, sizeof(float));
			fclose(file);

			for (int j=0, k=0; j<N; j++) {
				if (r[j] >= m*max && isPeak(r, j, j<w?j:w, N-j-1<w?N-j-1:w)) {
					float min = max;
					for (int i=0; i<n; i++) if (min >= R[i]) k=i, min=R[i];
					if (R[k] < r[j]) {
						for (int i=k; i<n-1; i++) L[i]=L[i+1], R[i]=R[i+1];
						L[n-1]=l[j], R[n-1]=r[j];
					}
				}
			}

			printf("%f", getNumber(argv[i]));
			if (f) for (int j=n-1; j>=0; j--) printf("\t%f", L[j] ? f/L[j] : 0);
			else for (int j=n-1; j>=0; j--) printf("\t%f", L[j]);
			printf("\n");
		}

		else {
			printf("File not found.\n");
			exit(0);
		}
	}
}
