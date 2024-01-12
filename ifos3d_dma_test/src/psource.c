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
 *   generate explosive source at source nodes
 *  ----------------------------------------------------------------------*/

#include "fd.h"

void psource(int nt, float *Fs, float **  srcpos_loc, float ** signals, int nsrc) {


    extern float DX, DY, DZ, ABS_TYPE;
	extern int  FDORDER,  FDCOEFF;//,LOG,
    extern int NX, NY, NZ;
    extern int NT;
	int i, j, k, l, idx;
	float amp=0.0;

	int ll = 1;
	if(ABS_TYPE==1 && FDORDER==2){ll=2;}

	int nrow = NY + 2 * (ll*FDORDER/2) + 1, ncol = NX + 2 * (ll*FDORDER/2), ndep = NZ + 2 * (ll*FDORDER/2);
	int size = nrow * ncol * ndep;
	int strip = ndep;
	int slice = ncol * ndep;

	/* adding source wavelet to stress components 
	   (explosive source) at source points */


	/*
	for (l=1;l<=nsrc;l++) {
		if((int)srcpos_loc[7][l]==1){
		i=(int)srcpos_loc[1][l];
		j=(int)srcpos_loc[2][l];
		k=(int)srcpos_loc[3][l];
			if(nt==1){amp=signals[l][nt+1]/(2.0*DX*DY*DZ);}
			if((nt>1)&&(nt<NT)){amp=(signals[l][nt+1]-signals[l][nt-1])/(2.0*DX*DY*DZ);}
			if(nt==NT){amp=-signals[l][nt-1]/(2.0*DX*DY*DZ);}
		sxx[j][i][k]+=amp;
		syy[j][i][k]+=amp;
		szz[j][i][k]+=amp;
	    }
	}
	*/

	for (l=1;l<=nsrc;l++) {
		if((int)srcpos_loc[7][l]==1) {
			i=(int)srcpos_loc[1][l];
			j=(int)srcpos_loc[2][l];
			k=(int)srcpos_loc[3][l];
			if(nt==1) {amp=signals[l][nt+1]/(2.0*DX*DY*DZ);}
			if((nt>1)&&(nt<NT)) {amp=(signals[l][nt+1]-signals[l][nt-1])/(2.0*DX*DY*DZ);}
			if(nt==NT) {amp=-signals[l][nt-1]/(2.0*DX*DY*DZ);}
			idx = j * slice + i * strip + k;
			Fs[6 * idx + 0] += amp;
			Fs[6 * idx + 1] += amp;
			Fs[6 * idx + 2] += amp;
	    }
	}


}
