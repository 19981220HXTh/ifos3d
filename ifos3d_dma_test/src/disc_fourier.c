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

/*-------------------------------------------------------------------------
 * calculation of a discrete Fourier transfomation on the fly:
 * Fouriercomponents of forward or backpropagated wavefields are summed up for each frequency
 * S. Butezr 2013
 --------------------------------------------------------------------------*/
 

#include "fd.h"
void discfourier(int nx1, int nx2, int ny1, int ny2, int nz1, int nz2, int nt, float *Fv, float *F_Fv, float *finv, int nf, int ntast, int pshot1, int back) {
  

	extern float DT,TIME;
	extern int NX, NY, NZ, BOUNDARY, FDORDER, NFMAX, ABS_TYPE, MYID;
	int i, j, k, l,m,idx,idx_f;
	double trig1,trig2;
	float t=0.0;
	//printf("--#------------------>   FUNCTION:%s, LINE:%d, MYID:%d\n", __FUNCTION__, __LINE__, MYID);
	if(back==0) t=nt*DT;
	if(back==1) t=TIME-nt*DT;
	
	int ll = 1;
	if(ABS_TYPE==1 && FDORDER==2){ll=2;}

	int nval = NFMAX, nrow = NY + 2 * (ll*FDORDER/2) + 1, ncol = NX + 2 * (ll*FDORDER/2), ndep = NZ + 2 * (ll*FDORDER/2);



    int strip = ndep;
    int slice = ncol * ndep;
    int cube = nrow * ncol * ndep;



	for(l=1;l<=nf;l++){
		m=(pshot1)*nf+l;
		trig1=0.0;
		trig1=cos(2.0*t*finv[l-1]*M_PI)*DT*ntast;
		trig2=0.0;
		trig2=sin(2.0*t*finv[l-1]*M_PI)*DT*ntast;

		for (j=ny1;j<=ny2;j++){
			for (i=nx1;i<=nx2;i++){
				for (k=nz1;k<=nz2;k++){
					idx = j * slice + i * strip + k;
					idx_f = m * cube + idx;
                    F_Fv[idx_f * 6 + 0] += Fv[idx * 3 + 0] * trig1;
                    F_Fv[idx_f * 6 + 1] += Fv[idx * 3 + 1] * trig1;
                    F_Fv[idx_f * 6 + 2] += Fv[idx * 3 + 2] * trig1;


                    F_Fv[idx_f * 6 + 3] += Fv[idx * 3 + 0] * trig2;
                    F_Fv[idx_f * 6 + 4] += Fv[idx * 3 + 1] * trig2;
                    F_Fv[idx_f * 6 + 5] += Fv[idx * 3 + 2] * trig2;


				}
			}
		}
	}



}
