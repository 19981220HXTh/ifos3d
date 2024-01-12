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
 *   Write 3D snapshot for current timestep  to disk                                   
 *  ----------------------------------------------------------------------*/

#include "fd.h"


void snap(FILE *fp, int nt, int nsnap, int format, int type, float *Fv, float *Fs, float ***u, float ***pi, int idx, int idy, int idz, int nx1, int ny1, int nz1, int nx2, int ny2, int nz2) {

	/* 
	different data formats of output available:
	format=1  :  SU (IEEE)
	format=2  :  ASCII
	format=3  :  BINARY (IEEE)
	
	different types:
	type=1 : values in vx, vy, and vz
	type=2 : -(sxx+syy+szz) (pressure field)
	type=3 : divergence of vx, vy and vz (energy of compressional waves)
	         and curl of vx, vy and vz (energy of shear waves)
	type=4 : both particle velocities (type=1) and energy (type=3)
	*/


	
	char xfile[STRING_SIZE], yfile[STRING_SIZE], zfile[STRING_SIZE];
	char rotfile[STRING_SIZE], ext[8], wm[1];
	char  divfile[STRING_SIZE], pfile[STRING_SIZE];
	FILE *fpx1, *fpy1, *fpz1, *fpx2, *fpy2;
	int i,j,k, idx_n;
	float a=0.0, amp, dh24x, dh24y, dh24z, vyx, vxy, vxx, vyy, vzx, vyz, vxz, vzy, vzz;
	/*extern FILE *FP;*/

	extern float DX, DY, DZ, DT;
	extern char SNAP_FILE[STRING_SIZE];
	extern int MYID, POS[4], SNAP_PLANE, LOG, ABS_TYPE;
	extern int  FDORDER,  FDCOEFF;//,LOG,
    extern int NX, NY, NZ;


	int ll = 1;
	if(ABS_TYPE==1 && FDORDER==2){ll=2; printf("ABS:%d, FDOR:%d\n", ABS_TYPE, FDORDER);}
	int nrow = NY + 2 * (ll*FDORDER/2) + 1, ncol = NX + 2 * (ll*FDORDER/2), ndep = NZ + 2 * (ll*FDORDER/2);
	int size = nrow * ncol * ndep;
	int strip = ndep;
	int slice = ncol * ndep;


	switch(format){
	case 1: 
		sprintf(ext,".su");
		break;
	case 2: 
		sprintf(ext,".asc");
		break;
	case 3: 
		sprintf(ext,".bin");
		break;
	}


	sprintf(xfile,"%s%s.x.%i%i%i",SNAP_FILE,ext,POS[1],POS[2],POS[3]);
	sprintf(yfile,"%s%s.y.%i%i%i",SNAP_FILE,ext,POS[1],POS[2],POS[3]);
	sprintf(zfile,"%s%s.z.%i%i%i",SNAP_FILE,ext,POS[1],POS[2],POS[3]);
	sprintf(divfile,"%s%s.div.%i%i%i",SNAP_FILE,ext,POS[1],POS[2],POS[3]);
	sprintf(rotfile,"%s%s.rot.%i%i%i",SNAP_FILE,ext,POS[1],POS[2],POS[3]);
	sprintf(pfile,"%s%s.p.%i%i%i",SNAP_FILE,ext,POS[1],POS[2],POS[3]);

    if (LOG) {fprintf(fp,"\n\n PE %d is writing snapshot-data at T=%fs to \n",MYID,nt*DT);}

      

	if (nsnap==1) 
		sprintf(wm,"w");
	else 
		sprintf(wm,"a");

	switch(type){
	case 1 :
		fprintf(fp,"\t%s\n", xfile);
		fprintf(fp,"\t%s\n", yfile);
		fprintf(fp,"\t%s\n\n", zfile);
		fpx1=fopen(xfile,wm);
		fpy1=fopen(yfile,wm);
		fpz1=fopen(zfile,wm);
		
		/*for (k=ny1;k<=ny2;k+=idy)
			for (i=nx1;i<=nx2;i+=idx)
				for (j=nz1;j<=nz2;j+=idz){*/
		/*
		for (j=ny1;j<=ny2;j+=idy){
		for (i=nx1;i<=nx2;i+=idx){
			for (k=nz1;k<=nz2;k+=idz){
					writedsk(fpx1,vx[j][i][k],format);
					writedsk(fpy1,vy[j][i][k],format);
					writedsk(fpz1,vz[j][i][k],format);
					}}
				}
		*/
		for (j=ny1;j<=ny2;j+=idy) {
			for (i=nx1;i<=nx2;i+=idx ){
				for (k=nz1;k<=nz2;k+=idz) {
					int idx_n = j * slice + i * strip + k;
					writedsk(fpx1,Fv[idx_n * 3 + 0],format);
					writedsk(fpy1,Fv[idx_n * 3 + 1],format);
					writedsk(fpz1,Fv[idx_n * 3 + 2],format);
				}
			}
		}

		fclose(fpx1);
		fclose(fpy1);
		fclose(fpz1);
		break;

	}

}


