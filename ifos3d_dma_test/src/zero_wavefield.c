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

/* -----------------------------------------------------------------------
 * Initialise wavefield with zero
 -------------------------------------------------------------------------*/

#include "fd.h"
#include <unistd.h>
void zero_wavefield( int NX, int NY, int NZ, float *Fv, float *Fs, float *** rxx, float *** ryy, float *** rzz, float *** rxy, float *** ryz, float *** rxz, float *** psi_sxx_x, float *** psi_sxy_x, float *** psi_sxz_x, float *** psi_sxy_y, float *** psi_syy_y, float *** psi_syz_y, float *** psi_sxz_z, float *** psi_syz_z, float *** psi_szz_z, float *** psi_vxx, float *** psi_vyx, float *** psi_vzx, float *** psi_vxy, float *** psi_vyy, float *** psi_vzy, float *** psi_vxz, float *** psi_vyz, float *** psi_vzz) {


    extern int FDORDER, ABS_TYPE, FW, POS[4],L, MYID;
    int nx1, ny1, nz1, nx2, ny2, nz2, a,b,l, i, j, k;


	int ll = 1;
	if(ABS_TYPE==1 && FDORDER==2) {ll=2; printf("ABS:%d, FDOR:%d\n", ABS_TYPE, FDORDER);}
	int nrow = NY + 2 * (ll*FDORDER/2) + 1, ncol = NX + 2 * (ll*FDORDER/2), ndep = NZ + 2 * (ll*FDORDER/2);
    int size = nrow * ncol * ndep;



	int strip = ndep;
	int slice = ncol * ndep;
	int strip_r = NZ;
	int slice_r = NX * NZ;
	int strip_z = 2 * FW;
	int slice_z = NX * 2 * FW;
	int strip_x = NZ;
	int slice_x = 2 * FW * NZ;
	int strip_y = NZ;
	int slice_y = NX * NZ;

    l=1;
    if(ABS_TYPE==1 && FDORDER==2){l=2;}
	
	if(POS[2]==0){
	  a=0;
	  b=1;}
	else{
	  a=1;
	  b=l;}
	

    ny1=0-l*FDORDER/2;
    ny2=NY+l*FDORDER/2;
    nx1=1-l*FDORDER/2;
    nx2=NX+l*FDORDER/2;
    nz1=1-l*FDORDER/2;
    nz2=NZ+l*FDORDER/2;


    for (j = ny1; j <= ny2; j++) {
    	for (i = nx1; i <= nx2; i++) {
    		for (k = nz1; k <= nz2; k++) {
    			int idx = j * slice + i * strip + k;
    			Fv[idx * 3 + 0] = 0.0;
    			Fv[idx * 3 + 1] = 0.0;
    			Fv[idx * 3 + 2] = 0.0;
    			Fs[idx * 6 + 0] = 0.0;
    			Fs[idx * 6 + 1] = 0.0;
    			Fs[idx * 6 + 2] = 0.0;
    			Fs[idx * 6 + 3] = 0.0;
    			Fs[idx * 6 + 4] = 0.0;
    			Fs[idx * 6 + 5] = 0.0;

    		}
    	}
    }
	

    if (L) {
    	float *Fr = transform_6(rxx, ryy, rzz, rxy, ryz, rxz, 1, NY, 1, NX, 1, NZ);
    	for (j=1;j<=NY;j++) {
    	    for (i=1;i<=NX;i++) {
    	    	for (k=1;k<=NZ;k++) {
    	    		int idx_r = j * slice_r + i * strip_r + k;
    	    		Fr[idx_r * 6 + 0] = 0.0;
    	    		Fr[idx_r * 6 + 1] = 0.0;
    	    		Fr[idx_r * 6 + 2] = 0.0;
    	    		Fr[idx_r * 6 + 3] = 0.0;
    	    		Fr[idx_r * 6 + 4] = 0.0;
    	    		Fr[idx_r * 6 + 5] = 0.0;
    	    	}
    	    }
    	}

        inverse_6(Fr, rxx, ryy, rzz, rxy, ryz, rxz, 1, NY, 1, NX, 1, NZ);
        free_trans_6(Fr, 1, NY, 1, NX, 1, NZ);
    }
    

    if(ABS_TYPE==1){

    	float *F_psi_z = transform_6(psi_sxz_z, psi_syz_z, psi_szz_z, psi_vxz, psi_vyz, psi_vzz, 1, NY, 1, NX, 1, 2 * FW);
    	float *F_psi_x = transform_6(psi_sxx_x, psi_sxy_x, psi_sxz_x, psi_vxx, psi_vyx, psi_vzx, 1, NY, 1, 2 * FW, 1, NZ);
    	float *F_psi_y = transform_6(psi_sxy_y, psi_syy_y, psi_syz_y, psi_vxy, psi_vyy, psi_vyz, 1, 2 * FW, 1, NX, 1, NZ);

    	for (j=1;j<=NY;j++){
    		for (i=1;i<=NX;i++){
    			for (k=1;k<=2*FW;k++){
    				int idx_z = j * slice_z + i * strip_z + k;
    	    		F_psi_z[idx_z * 6 + 0] = 0.0;
    	    		F_psi_z[idx_z * 6 + 1] = 0.0;
    	    		F_psi_z[idx_z * 6 + 2] = 0.0;
    	    		F_psi_z[idx_z * 6 + 3] = 0.0;
    	    		F_psi_z[idx_z * 6 + 4] = 0.0;
    	    		F_psi_z[idx_z * 6 + 5] = 0.0;

    			}
    		}
    	}

    	for (j=1;j<=NY;j++) {
    		for (i=1;i<=2*FW;i++) {
    			for (k=1;k<=NZ;k++){
    				int idx_x = j * slice_x + i * strip_x + k;
    	    		F_psi_x[idx_x * 6 + 0] = 0.0;
    	    		F_psi_x[idx_x * 6 + 1] = 0.0;
    	    		F_psi_x[idx_x * 6 + 2] = 0.0;
    	    		F_psi_x[idx_x * 6 + 3] = 0.0;
    	    		F_psi_x[idx_x * 6 + 4] = 0.0;
    	    		F_psi_x[idx_x * 6 + 5] = 0.0;

    			}
    		}
    	}
    	for (j=1;j<=2*FW;j++){
    		for (i=1;i<=NX;i++){
    			for (k=1;k<=NZ;k++){
    				int idx_y = j * slice_y + i * strip_y + k;
    	    		F_psi_y[idx_y * 6 + 0] = 0.0;
    	    		F_psi_y[idx_y * 6 + 1] = 0.0;
    	    		F_psi_y[idx_y * 6 + 2] = 0.0;
    	    		F_psi_y[idx_y * 6 + 3] = 0.0;
    	    		F_psi_y[idx_y * 6 + 4] = 0.0;
    	    		F_psi_y[idx_y * 6 + 5] = 0.0;

    			}
    		}
    	}

        inverse_6(F_psi_x, psi_sxx_x, psi_sxy_x, psi_sxz_x, psi_vxx, psi_vyx, psi_vzx, 1, NY, 1, 2 * FW, 1, NZ);
        inverse_6(F_psi_y, psi_sxy_y, psi_syy_y, psi_syz_y, psi_vxy, psi_vyy, psi_vzy, 1, 2 * FW, 1, NX, 1, NZ);
        inverse_6(F_psi_z, psi_sxz_z, psi_syz_z, psi_szz_z, psi_vxz, psi_vyz, psi_vzz, 1, NY, 1, NX, 1, 2 * FW);
        free_trans_6(F_psi_x, 1, NY, 1, 2 * FW, 1, NZ);
        free_trans_6(F_psi_y, 1, 2 * FW, 1, NX, 1, NZ);
        free_trans_6(F_psi_z, 1, NY, 1, NX, 1, 2 * FW);
    }




}
