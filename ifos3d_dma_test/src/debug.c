/*
 * debug.c
 *
 *  Created on: Oct 17, 2017
 *      Author: bingo
 */


#include "fd.h"
#include <stdbool.h>

void debug_v_3(float ***a, float ***b, float ***c, int nx1, int nx2, int ny1, int ny2, int nz1, int nz2, bool tag) {
	int i, j, k;
	printf("nx1:%dn nx2:%d, ny1:%d, ny2:%d, nz1:%d, nz2:%d\n", nx1, nx2, ny1, ny2, nz1, nz2);
	for (j =  ny1; j <= ny2; j++) {
		for (i = nx1; i <= nx2; i++) {
			for (k = nz1; k <= nz2; k++) {
				printf("--> %s: %s: line:%d , j:%d, i:%d, k:%d, vx:%f, vy:%f, vz:%f\n",
											__FILE__, __FUNCTION__, __LINE__, j, i, k, a[j][i][k], b[j][i][k], c[j][i][k]);
				/*
				if (tag == true){
					if (j == ny2 && i == nx2) {
						printf("--> %s: %s: line:%d , j:%d, i:%d, k:%d, vx:%f, vy:%f, vz:%f\n",
						       __FILE__, __FUNCTION__, __LINE__, j, i, k, a[j][i][k], b[j][i][k], c[j][i][k]);
					}
				} else {
					printf("--> %s: %s: line:%d , j:%d, i:%d, k:%d, vx:%f, vy:%f, vz:%f\n",
							__FILE__, __FUNCTION__, __LINE__, j, i, k, a[j][i][k], b[j][i][k], c[j][i][k]);
				}*/

			}
		}
	}

}


void debug_s_3(float ***a, float ***b, float ***c, float ***d, float ***e, float ***f,  int nx1, int nx2, int ny1, int ny2, int nz1, int nz2, bool tag) {
	int i, j, k;
	for (j =  ny1; j <= ny2; j++) {
		for (i = nx1; i <= nx2; i++) {
			for (k = nz1; k <= nz2; k++) {
				//printf("j:%d, i:%d, k:%d, vx:%f, vy:%f, vz:%f, Fvx:%f, fvy:%f, fvz:%f, --> Frp: %f\n",
						//j, i, k, vx[j][i][k], vy[j][i][k], vz[j][i][k], Fv[idx * 3 + 0], Fv[idx * 3 + 1], Fv[idx * 3 + 2] ,Frp[idx * 3 + 2]);
			    //printf("-->%s: %s: line:%d , j:%d, i:%d, k:%d, vx:%f, vy:%f, vz:%f, Fvx:%f, fvy:%f, fvz:%f, --> Frp: %f\n",
			    		//__FILE__, __FUNCTION__, __LINE__, j, i, k, vx[j][i][k], vy[j][i][k], vz[j][i][k], Fv[idx * 3 + 0], Fv[idx * 3 + 1], Fv[idx * 3 + 2] ,Frp[idx * 3 + 2]);

			    //printf("j:%d, i:%d, k:%d, sxx:%f, syy:%f, szz:%f, sxy:%f, syz:%f, sxz:%f\n",
				       //j, i, k, Fs[idx * 3 + 0], Fs[idx * 3 + 1], Fs[idx * 3 + 2] , Fs[idx * 3 + 3], Fs[idx * 3 + 4], Fs[idx * 3 + 5]);
				if (tag == true){
					if (j == ny2 && i == nx2) {
						printf("####>> %s: %s: line:%d , j:%d, i:%d, k:%d, sxx:%f, syy:%f, szz:%f, sxy:%f, syz:%f, sxz:%f\n",
								__FILE__, __FUNCTION__, __LINE__, j, i, k, a[j][i][k],  b[j][i][k], c[j][i][k], d[j][i][k], e[j][i][k], f[j][i][k]);
					}
				} else {
					printf("####>> %s: %s: line:%d , j:%d, i:%d, k:%d, sxx:%f, syy:%f, szz:%f, sxy:%f, syz:%f, sxz:%f\n",
							__FILE__, __FUNCTION__, __LINE__, j, i, k, a[j][i][k],  b[j][i][k], c[j][i][k], d[j][i][k], e[j][i][k], f[j][i][k]);
				}



			}
		}
	}

}
