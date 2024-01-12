#include "fd.h"


float *transform3(float ***a, int nrl, int nrh, int ncl, int nch, int ndl, int ndh) {
	int i, j, d, idx, nrow = nrh - nrl + 1, ncol = nch - ncl + 1, ndep = ndh - ndl + 1;
	float *p = (float *)malloc((size_t)((nrow * ncol * ndep) * sizeof(float)));
	/*
	for (i = nrl; i <= nrh; i++)
	    for (j = ncl; j <= nch; j++)
	        for (d = ndl; d <= ndh; d++) printf("i:%d, j:%d, k:%d, value:%f\n", i, j, d, a[i][j][d]);
	*/
	//p += 1;
	p -= (nrl * ndep * ncol +  ncl * ndep + ndl);
	//idx = 0;
	for (i = nrl; i <= nrh; i++)
		    for (j = ncl; j <= nch; j++)
		        for (d = ndl; d <= ndh; d++) {
		        	idx = i * ndep * ncol + j * ndep + d;
		        	p[idx] = a[i][j][d];
		        	//idx++;
		        }

	return p;
}


float *transform_3(float ***a, float ***b, float ***c, int nrl, int nrh, int ncl, int nch, int ndl, int ndh) {
	int i, j, d, idx, nrow = nrh - nrl + 1, ncol = nch - ncl + 1, ndep = ndh - ndl + 1;
	float *p = (float *)malloc((size_t)((3 * nrow * ncol * ndep) * sizeof(float)));
	/*
	for (i = nrl; i <= nrh; i++)
	    for (j = ncl; j <= nch; j++)
	        for (d = ndl; d <= ndh; d++) printf("i:%d, j:%d, k:%d, value:%f\n", i, j, d, a[i][j][d]);
	*/
	//p += 1;
	p -= 3 * (nrl * ndep * ncol +  ncl * ndep + ndl);
	//idx = 0;
	for (i = nrl; i <= nrh; i++)
		    for (j = ncl; j <= nch; j++)
		        for (d = ndl; d <= ndh; d++) {
		        	idx = i * ndep * ncol + j * ndep + d;
		        	p[idx * 3 + 0] = a[i][j][d];
		        	p[idx * 3 + 1] = b[i][j][d];
		        	p[idx * 3 + 2] = c[i][j][d];
		        }

	return p;
}


float *transform_6(float ***a, float ***b, float ***c, float ***e, float ***f, float ***g, int nrl, int nrh, int ncl, int nch, int ndl, int ndh) {
	int i, j, d, idx, nrow = nrh - nrl + 1, ncol = nch - ncl + 1, ndep = ndh - ndl + 1;
	float *p = (float *)malloc((size_t)((6 * nrow * ncol * ndep) * sizeof(float)));
	/*
	for (i = nrl; i <= nrh; i++)
	    for (j = ncl; j <= nch; j++)
	        for (d = ndl; d <= ndh; d++) printf("i:%d, j:%d, k:%d, value:%f\n", i, j, d, a[i][j][d]);
	*/
	//p += 1;
	p -= 6 * (nrl * ndep * ncol +  ncl * ndep + ndl);
	//idx = 0;
	for (i = nrl; i <= nrh; i++)
		    for (j = ncl; j <= nch; j++)
		        for (d = ndl; d <= ndh; d++) {
		        	idx = i * ndep * ncol + j * ndep + d;
		        	p[idx * 6 + 0] = a[i][j][d];
		        	p[idx * 6 + 1] = b[i][j][d];
		        	p[idx * 6 + 2] = c[i][j][d];
		        	p[idx * 6 + 3] = e[i][j][d];
		        	p[idx * 6 + 4] = f[i][j][d];
		        	p[idx * 6 + 5] = g[i][j][d];
		        }

	return p;
}

//fuse the 6 array which has 4 dimensions into a array only has one dimension
float *transform4_6(float ****a, float ****b, float ****c, float ****e, float ****f, float ****g, int nvl, int nvh, int nrl, int nrh, int ncl, int nch, int ndl, int ndh) {
	int i, j, d, v, idx, nval=nvh-nvl+1, nrow=nrh-nrl+1,ncol=nch-ncl+1,ndep=ndh-ndl+1;
	float *p = (float *)malloc((size_t)((6 * nval * nrow * ncol * ndep) * sizeof(float)));
	/*
	for (i = nrl; i <= nrh; i++)
	    for (j = ncl; j <= nch; j++)
	        for (d = ndl; d <= ndh; d++) printf("i:%d, j:%d, k:%d, value:%f\n", i, j, d, a[i][j][d]);
	*/
	//p += 1;
	p -= 6 * (nvl * nrow * ncol * ndep + nrl * ndep * ncol +  ncl * ndep + ndl);
	//idx = 0;
	for (v = nvl; v <= nvh; v++)
		for (i = nrl; i <= nrh; i++)
		    for (j = ncl; j <= nch; j++)
		        for (d = ndl; d <= ndh; d++) {
		        	idx = v * nrow * ncol * ndep + i * ncol * ndep + j * ndep + d;
		        	p[idx * 6 + 0] = a[v][i][j][d];
		        	p[idx * 6 + 1] = b[v][i][j][d];
		        	p[idx * 6 + 2] = c[v][i][j][d];
		        	p[idx * 6 + 3] = e[v][i][j][d];
		        	p[idx * 6 + 4] = f[v][i][j][d];
		        	p[idx * 6 + 5] = g[v][i][j][d];
		        }

	return p;
}


void inversion(float *p, float ***a, int nrl, int nrh, int ncl, int nch, int ndl, int ndh) {
	int i, j, d, idx;
	idx = 0;
	for (i = nrl; i <= nrh; i++)
			    for (j = ncl; j <= nch; j++)
			        for (d = ndl; d <= ndh; d++) {
			        	a[i][j][d] = p[idx];
			        	idx++;
			        }
}


void inverse_3(float *Fv, float ***vx, float ***vy, float ***vz, int nrl, int nrh, int ncl, int nch, int ndl, int ndh) {
	int i, j, d, idx, nrow = nrh - nrl + 1, ncol = nch - ncl + 1, ndep = ndh - ndl + 1;

	/*printf("row:%d, col:%d, dep:%d, row:%d, col:%d, dep:%d\n",  nrl, nrh,  ncl, nch, ndl, ndh);
	for (i = nrl; i <= nrh; i++)
		for (j = ncl; j <= nch; j++)
		    for (d = ndl; d <= ndh; d++) {
				idx = i * ndep * ncol + j * ndep + d;
				if (Fv[idx * 3 + 0] != 0)
				     printf("--> i:%d, j:%d, d:%d, p1:%f, p2:%f, p3:%f\n", i, j, d, Fv[3 * idx + 0], Fv[3 * idx + 1], Fv[3 * idx + 2]);

			}*/

	for (i = nrl; i <= nrh; i++)
		for (j = ncl; j <= nch; j++)
		    for (d = ndl; d <= ndh; d++) {
				idx = i * ndep * ncol + j * ndep + d;
				vx[i][j][d] = Fv[idx * 3 + 0];
				vy[i][j][d] = Fv[idx * 3 + 1];
				vz[i][j][d] = Fv[idx * 3 + 2];
				//if (Fv[idx * 3 + 0] != 0)
				     //printf("--> i:%d, j:%d, d:%d, p1:%f, p2:%f, p3:%f\n", i, j, d, Fv[3 * idx + 0], Fv[3 * idx + 1], Fv[3 * idx + 2]);

			}
}


void inverse_6(float *p, float ***a, float ***b, float ***c, float ***e, float ***f, float ***g, int nrl, int nrh, int ncl, int nch, int ndl, int ndh) {
	int i, j, d, idx, ncol = nch - ncl + 1, ndep = ndh - ndl + 1;
	for (i = nrl; i <= nrh; i++)
		for (j = ncl; j <= nch; j++)
		    for (d = ndl; d <= ndh; d++) {
				idx = i * ndep * ncol + j * ndep + d;
				a[i][j][d] = p[idx * 6 + 0];
				b[i][j][d] = p[idx * 6 + 1];
				c[i][j][d] = p[idx * 6 + 2];
				e[i][j][d] = p[idx * 6 + 3];
				f[i][j][d] = p[idx * 6 + 4];
				g[i][j][d] = p[idx * 6 + 5];
				//printf("--> i:%d, j:%d, d:%d, p1:%f, p2:%f, p3:%f\n", i, j, d, p[3 * idx + 0], p[3 * idx + 1], p[3 * idx + 2]);

			}
}



void inverse_6_4(float *p, float ****a, float ****b, float ****c, float ****e, float ****f, float ****g, int nvl, int nvh, int nrl, int nrh, int ncl, int nch, int ndl, int ndh) {
	int v,i, j, d, idx, nrow = nrh - nrl + 1, ncol = nch - ncl + 1, ndep = ndh - ndl + 1;
	for (v = nvl; v <= nvh; v++)
		for (i = nrl; i <= nrh; i++)
		    for (j = ncl; j <= nch; j++)
		        for (d = ndl; d <= ndh; d++) {
		        	idx = v * nrow * ncol * ndep + i * ncol * ndep + j * ndep + d;
		        	a[v][i][j][d] = p[idx * 6 + 0];
		        	b[v][i][j][d] = p[idx * 6 + 1];
		        	c[v][i][j][d] = p[idx * 6 + 2];
		        	e[v][i][j][d] = p[idx * 6 + 3];
		        	f[v][i][j][d] = p[idx * 6 + 4];
		        	g[v][i][j][d] = p[idx * 6 + 5];
		        }
}




void free_trans(float *a, int nrl, int nrh, int ncl, int nch, int ndl, int ndh) {
	int ncol = nch - ncl + 1, ndep = ndh - ndl + 1;
	free(a + (nrl * ndep * ncol +  ncl * ndep + ndl));
}



void free_trans_3(float *a, int nrl, int nrh, int ncl, int nch, int ndl, int ndh) {
	int ncol = nch - ncl + 1, ndep = ndh - ndl + 1;
	free(a + 3 * (nrl * ndep * ncol +  ncl * ndep + ndl));
}


void free_trans_6(float *a, int nrl, int nrh, int ncl, int nch, int ndl, int ndh) {
	int ncol = nch - ncl + 1, ndep = ndh - ndl + 1;
	free(a + 6 * (nrl * ndep * ncol +  ncl * ndep + ndl));
}


void free_trans_6_4(float *a, int nvl, int nvh, int nrl, int nrh, int ncl, int nch, int ndl, int ndh) {
	int nval = nvh - nvl + 1, nrow = nrh - nrl + 1, ncol = nch - ncl + 1, ndep = ndh - ndl + 1;
	free(a + 6 * (nvl * nrow * ncol * ndep + nrl * ndep * ncol +  ncl * ndep + ndl));
}


float *fuse_3(float *a, float *b, float *c, int size) {
    float *p = (float *)malloc((size_t)(size * 3 * sizeof(float)));
    int i;
    for (i = 0; i < size; i++) {
        p[i * 3 + 0] = a[i];
        p[i * 3 + 1] = b[i];
        p[i * 3 + 2] = c[i];
    }

    return p;
}


float *fuse_6(float *a0, float *a1, float *a2, float *a3, float *a4, float *a5, int size) {
    float *p = (float *)malloc((size_t)(size * 6 * sizeof(float)));
    int i;
    for (i = 0; i < size; i++) {
        p[i * 6 + 0] = a0[i];
        p[i * 6 + 1] = a1[i];
        p[i * 6 + 2] = a2[i];
        p[i * 6 + 3] = a3[i];
        p[i * 6 + 4] = a4[i];
        p[i * 6 + 5] = a5[i];
    }

    return p;
}


void unfuse_3(float *p, float *a, float *b, float *c, int size) {
	int i;
    for (i = 0; i < size; i++) {
        a[i] = p[i * 3 + 0];
        b[i] = p[i * 3 + 1];
        c[i] = p[i * 3 + 2];
    }
}


void unfuse_6(float *p, float *a, float *b, float *c, float *d, float *e, float *f, int size) {
	int i;
    for (i = 0; i < size; i++) {
        a[i] = p[i * 6 + 0];
        b[i] = p[i * 6 + 1];
        c[i] = p[i * 6 + 2];
        d[i] = p[i * 6 + 3];
        e[i] = p[i * 6 + 4];
        f[i] = p[i * 6 + 5];

    }
}
