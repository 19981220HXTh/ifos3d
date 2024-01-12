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
 *   updating velocity values by a staggered grid finite difference scheme of
 *   nth order accuracy in space and second order accuracy in time
 *  ----------------------------------------------------------------------*/

#include "fd.h"
#include<sys/time.h>

#define MASK




extern SLAVE_FUN(update_v_kernel_fusion_slave)(Param_vel *);

double update_v(int nx1, int nx2, int ny1, int ny2, int nz1, int nz2, int nt, float *Fv, float *Fs, float  ***  rho,  float  *** rjp, float  *** rkp, float  *** rip, float **  srcpos_loc, float ** signals, float ** signaly, float ** signalz, int nsrc, float *** absorb_coeff, int back, float *** buffertop_to_bot, MPI_Request * req_send, MPI_Request * req_rec) {


	/*extern FILE *FP;*/
	extern float DT, DX, DY, DZ, ALPHA, BETA;
	extern int NX, NY, NZ;
	double time=0.0; /*, time1=0.0;*/
	/*double time2=0.0;*/
	extern int  FDORDER,  ABS_TYPE, FDCOEFF; extern int MYID;  // ,LOG,*/



	int i, j, k, l;
	float  amp, alpha_rad, beta_rad;
	float b1, b2, b3, b4, b5, b6, dx, dy, dz;
	float sxx_x, sxy_y, sxz_z, syy_y, sxy_x, syz_z;
	float szz_z, sxz_x, syz_y;

	int ll = 1;
	if(ABS_TYPE==1 && FDORDER==2){ll=2; printf("ABS:%d, FDOR:%d\n", ABS_TYPE, FDORDER);}
	 //printf("debugging MYID:%d %s: ------------------------------------------------------------->\n", MYID,  __FUNCTION__);
    /************************************************* part of array fusion  ****************************************/

	//int nrow = NY + 2 * (ll*FDORDER/2) + 1;
	int ncol = NX + 2 * (ll*FDORDER/2), ndep = NZ + 2 * (ll*FDORDER/2);
    //int size = nrow * ncol * ndep;
	float *Frp = transform_3(rip, rjp, rkp, 1, NY, 1, NX, 1, NZ);
	int strip = ndep;
	int slice = ncol * ndep;
	int strip_rp = NZ;
	int slice_rp = NZ * NX;




	int dim_x = nx2 - nx1 + 1;
	int dim_y = ny2 - ny1 + 1;
	int dim_z = nz2 - nz1 + 1;
	//if(MYID == 0)
		//printf("---###-----------------------------------> nz1: %d, nz2:%d dim_z: %d\n", nz1, nz2, dim_z);


	Param_vel param;


	/*unsigned int * test_float;*/



	/*if (LOG)
	if (MYID==0) time1=MPI_Wtime();*/


    switch (FDORDER){

	case 4 :

	    dx=DT/DX;
		dy=DT/DY;
		dz=DT/DZ;
		b1=9.0/8.0; b2=-1.0/24.0; /* Taylor coefficients*/
		if(FDCOEFF==2){
		b1=1.1382; b2=-0.046414;} /* Holberg coefficients E=0.1 %*/

		// float *vel = Fv + 3 * (ny1 * slice + nx1 * strip + nz1);
		// float *stress = Fs + 6 * (ny1 * slice + nx1 * strip + nz1);
		// float *rp = Frp + 3 * (ny1 * slice_rp + nx1 * strip_rp + nz1) ;

#ifdef  MASK
        int ny;
		// if(NY % 2 == 0){
		// 	ny = NY/2;
		// }else{
		// 	ny = NY/2 + 1;
		// }
		ny = ((NY % 2 == 0) ? (NY/2):(NY/2+1));
		param.dim_x = dim_x;
		param.dim_y = ny;
		param.dim_z = dim_z;
		param.slice = slice;
		param.strip = strip;
		param.slice_rp = slice_rp;
		param.strip_rp = strip_rp;
		param.dx = dx;
		param.dy = dy;
		param.dz = dz;
		param.FDCOEFF = FDCOEFF;
		param.Fv = Fv;
		param.Fs = Fs;
		param.Frp = Frp;
		// gettimeofday(&start1, NULL);
		// athread_spawn(update_v_kernel_fusion_new_slave, &param);
	  	athread_spawn(update_v_kernel_fusion_slave, &param);
		athread_join();
		// gettimeofday(&end1, NULL);
		param.dim_y = dim_y - ny;
		param.Fv = Fv + 3 * ny * slice;
		param.Fs = Fs + 6 * ny * slice;
		param.Frp = Frp + 3 * ny * slice_rp;
		// gettimeofday(&start2, NULL);
		// athread_spawn(update_v_kernel_fusion_new_slave, &param);
	  	athread_spawn(update_v_kernel_fusion_slave, &param);
		// gettimeofday(&end2, NULL);
		exchange_v2_1(Fv, buffertop_to_bot, req_send, req_rec);
		// gettimeofday(&start3, NULL);
		athread_join();
	// 	gettimeofday(&end3, NULL);
	// 	time_update_v_1= (end1.tv_sec - start1.tv_sec)*1000000 + (end1.tv_usec-start1.tv_usec);
	// 	time_update_v_2= (end3.tv_sec - start2.tv_sec)*1000000 + (end3.tv_usec-start2.tv_usec);
	// 	time_exchange_v_1= (start3.tv_sec - end2.tv_sec)*1000000 + (start3.tv_usec-end2.tv_usec);
	// 	if (MYID == 0){
	// 	printf("time test--  update_v_1: %lf us, update_v_2: %lf us, exchange_v_1: %lf us\n",time_update_v_1, time_update_v_2, time_exchange_v_1);
	//   }
#else

		param.dim_x = dim_x;
		param.dim_y = dim_y;
		param.dim_z = dim_z;
		param.slice = slice;
		param.strip = strip;
		param.slice_rp = slice_rp;
		param.strip_rp = strip_rp;
		param.dx = dx;
		param.dy = dy;
		param.dz = dz;
		param.FDCOEFF = FDCOEFF;
		// param.vel = vel;
		// param.stress = stress;
		// param.rp = rp;
		param.Fv = Fv;
		param.Fs = Fs;
		param.Frp = Frp;

		//printf("---###------------------------------------------------------> %s,  %d,  MYID: %d\n", __FUNCTION__, __LINE__, MYID);

	  	athread_spawn(update_v_kernel_fusion_slave, &param);
	  	athread_join();
	  	//printf("---###------------------------------------------------------> %s,  %d,  MYID: %d\n", __FUNCTION__, __LINE__, MYID);
#endif	 
	  	break;


     }
	/* Adding body force components to corresponding particle velocities */



	if(back==0){
	for (l=1;l<=nsrc;l++) {
		i=(int)srcpos_loc[1][l];
		j=(int)srcpos_loc[2][l];
		k=(int)srcpos_loc[3][l];
		int idx = j * slice + i * strip + k;
		int idx_rp = j * slice_rp + i * strip_rp + k;
		amp=DT*signals[l][nt]/(DX*DY*DZ);


		switch ((int)srcpos_loc[7][l]){
		case 2 :
			Fv[idx * 3 + 0]+=amp*Frp[idx_rp * 3 + 0];

			/*/(rip[j][i][k]);  *//* single force in x , DX^3 because of force density*/

			/*test_float = (unsigned int*) (void*) &vx[j][i][k];
			if(((*test_float & 0x7f800000) == 0x0) && ((*test_float & 0x007fffff) != 0x0))fprintf(FP,"Achtung:ILLEGALER FLOAT");*/

			break;
		case 3 :
			Fv[idx * 3 + 2]+=amp*Frp[idx_rp * 3 + 1];  /* single force in z  */
			break;
		case 4 :
			Fv[idx * 3 + 1]+=amp*Frp[idx_rp * 3 + 2];  /* single force in y, vertical direction*/
			break;
		case 5 :
			alpha_rad=ALPHA*PI/180; /* custom force */
			beta_rad=BETA*PI/180;
			Fv[idx * 3 + 0]+=cos(beta_rad)*cos(alpha_rad)*amp;
			Fv[idx * 3 + 1]+=sin(alpha_rad)*amp;
			Fv[idx * 3 + 2]+=sin(beta_rad)*cos(alpha_rad)*amp;
			break;
		}
	}



	}

	if (back==1){

		for (l=1;l<=nsrc;l++) {
			i=(int)srcpos_loc[1][l];
			j=(int)srcpos_loc[2][l];
			k=(int)srcpos_loc[3][l];
			int idx = j * slice + i * strip + k;
			Fv[idx * 3 + 0]+=signals[l][nt];
			Fv[idx * 3 + 2]+=signalz[l][nt];
			Fv[idx * 3 + 1]+=signaly[l][nt];
		}
	}




	/* absorbing boundary condition (exponential damping) */

	if (ABS_TYPE==2){
		for (j=ny1;j<=ny2;j++){
			for (i=nx1;i<=nx2;i++){
				for (k=nz1;k<=nz2;k++){
					int idx = j * slice + i * strip + k;
					Fv[idx * 3 + 0]*=absorb_coeff[j][i][k];
					Fv[idx * 3 + 1]*=absorb_coeff[j][i][k];
					Fv[idx * 3 + 2]*=absorb_coeff[j][i][k];

					Fs[idx * 6 + 0]*=absorb_coeff[j][i][k];
					Fs[idx * 6 + 1]*=absorb_coeff[j][i][k];
					Fs[idx * 6 + 2]*=absorb_coeff[j][i][k];
					Fs[idx * 6 + 3]*=absorb_coeff[j][i][k];
					Fs[idx * 6 + 4]*=absorb_coeff[j][i][k];
					Fs[idx * 6 + 5]*=absorb_coeff[j][i][k];


			        }
		        }
	        }
        }


        /*if (LOG)
	if (MYID==0){
		time2=MPI_Wtime();
		time=time2-time1;
		fprintf(FP," Real time for particle velocity update: \t %4.2f s.\n",time);
	}*/


	/*
	for (j=1;j<=NY;j++){
		for (i=1;i<=NX;i++){
			for (k=1;k<=NZ;k++){
	            int idx = j * slice + i * strip + k;
	            int idx_rp = j * slice_rp + i * strip_rp + k;
	            if (vz[j][i][k] != 0)
	            	if(MYID == 0)
			    	printf("########>line: %d, j:%d, i:%d, k:%d, vz:%f, Fv:%f, sxx:%f, Fs:%f\n", __LINE__,  j, i, k, vz[j][i][k], Fv[idx * 3 + 2], sxx[j][i][k], Fs[idx * 6 + 0]);



		    //if (MYID == 0 && (szz_z + sxz_x + syz_y) != 0)
		    	//printf("s:%f\n",(szz_z + sxz_x + syz_y));
            	//printf("j:%d, i:%d, k:%d, vx:%f, vy:%f, vz: %f, rip:%f\n", j, i, k, vx[j][i][k], vy[j][i][k], vz[j][i][k], rip[j][i][k]);

		    }
	    }
    }
    */

	//debug_v_3(vx, vy, vz, 0-ll*FDORDER/2,NY+ll*FDORDER/2,1-ll*FDORDER/2,NX+ll*FDORDER/2,1-ll*FDORDER/2,NZ+ll*FDORDER/2, true);

	//printf("---###------------------------------------------------------> %s,  %d,  MYID: %d\n", __FUNCTION__, __LINE__, MYID);

    free_trans_3(Frp, 1, NY, 1, NX, 1, NZ );




	return time;

}
