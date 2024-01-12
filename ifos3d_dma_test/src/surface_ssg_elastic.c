/*------------------------------------------------------------------------
 * Copyright (C) 2015 For the list of authors, see file AUTHORS.
 *
 * This file is part of IFOS3D.
 * 
 * IFOS3D is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 2.0 of the License only.
 * 
 * IFOS3D is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with IFOS3D. See file COPYING and/or 
 * <http://www.gnu.org/licenses/gpl-2.0.html>.
--------------------------------------------------------------------------*/

/*------------------------------------------------------------------------
 *   stress free surface condition, elastic case
 *  ----------------------------------------------------------------------*/

#include "fd.h"

void surface_elastic(int ndepth, float *** u, float *** pi,	float *Fs, float *Fv) {

	int i, k ,j, fdoh,m,idx,idx_n,idx_p;
	float  vxx, vyy, vzz;
	float f, g, h;


	extern int NX, NY, NZ, FDORDER, FDCOEFF, ABS_TYPE;
	register float b1, b2, b3, b4, b5, b6;
	extern float DT, DX, DY, DZ;

	j=ndepth;     /* The free surface is located exactly in y=ndepth*dh meter!! */

	/*dthalbe=DT/2.0;*/
	fdoh=FDORDER/2;
	int ll = 1;
	if(ABS_TYPE==1 && FDORDER==2) {ll=2; printf("ABS:%d, FDOR:%d\n", ABS_TYPE, FDORDER);}
	int nrow = NY + 2 * (ll*FDORDER/2) + 1, ncol = NX + 2 * (ll*FDORDER/2), ndep = NZ + 2 * (ll*FDORDER/2);
    int size = nrow * ncol * ndep;

	float *Fpi = transform3(pi, 1, NY + 1, 1, NX + 1, 1, NZ + 1);
    float *Fu = transform3(u, 1, NY + 1, 1, NX + 1, 1, NZ + 1);
	int strip = ndep;
	int slice = ncol * ndep;
	int strip_n = NZ;
	int slice_n = NX * NZ;
	int strip_p = NZ + 1;
	int slice_p = (NX + 1) * (NZ + 1);

	switch (FDORDER){

	case 4 :

		b1=9.0/8.0; b2=-1.0/24.0; /* Taylor coefficients */

		if(FDCOEFF==2) {b1=1.1382; b2=-0.046414;} /* Holberg coefficients E=0.1 %*/



		for (i=1; i<=NX; i++) {
			for (k=1; k<=NZ; k++) {
				//Mirroring the components of the stress tensor to make
					  //a stress free surface (method of imaging, Levander, 1988)
				for (m=1; m<=fdoh; m++) {
					idx = j * slice + i * strip + k;
					idx_n = j * slice_n + i * strip_n + k;
					idx_p = j * slice_p + i * strip_p + k;
					Fs[idx * 6 + 1] = 0.0;
					Fs[(idx - m * slice) * 6 + 1] = - Fs[(idx + m * slice) * 6 + 1];
					Fs[(idx - m * slice) * 6 + 3] = - Fs[(idx + (m - 1) * slice) * 6 + 3];
					Fs[(idx - m * slice) * 6 + 4] = - Fs[(idx + (m - 1) * slice) * 6 + 4];
				}

				vxx = (b1 * (Fv[idx * 3 + 0] - Fv[(idx - strip) * 3 + 0]) + b2 * (Fv[(idx + strip) * 3 + 0] - Fv[(idx - 2 * strip) * 3 + 0])) / DX;
				vyy = (b1 * (Fv[idx * 3 + 1] - Fv[(idx - slice) * 3 + 1]) + b2 * (Fv[(idx + slice) * 3 + 1] - Fv[(idx - 2 * slice) * 3 + 1])) / DY;
				vzz = (b1 * (Fv[idx * 3 + 2] - Fv[(idx - 1) * 3 + 2]) + b2 * (Fv[(idx + 1) * 3 + 2] - Fv[(idx - 2) * 3 + 2])) / DZ;

				f = 2.0*Fu[idx_p];
				g = Fpi[idx_p];
				h=-(DT*(g-f)*(g-f)*(vxx+vzz)/g)-(DT*(g-f)*vyy);
				Fs[idx * 6 + 0]+=h;
				Fs[idx * 6 + 2]+=h;
			}

		}

		break;

	}


    free_trans(Fpi, 1, NY + 1, 1, NX + 1, 1, NZ + 1);
    free_trans(Fu, 1, NY + 1, 1, NX + 1, 1, NZ + 1);
}
