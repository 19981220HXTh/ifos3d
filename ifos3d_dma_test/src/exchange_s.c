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
 * exchange of stress values at grid boundaries between processors
 * when using the standard staggered grid
 * ----------------------------------------------------------------------*/

#include "fd.h"

double exchange_s(float *Fs, float *** bufferlef_to_rig, float *** bufferrig_to_lef, float *** buffertop_to_bot, float *** bufferbot_to_top, float *** bufferfro_to_bac, float *** bufferbac_to_fro, MPI_Request * req_send, MPI_Request * req_rec) {

	extern int NX, NY, NZ, POS[4], NPROCX, NPROCY, NPROCZ, BOUNDARY, FDORDER,  INDEX[7], ABS_TYPE, MYID;    /*MYID,LOG,*/
	extern const int TAG1,TAG2,TAG3,TAG4,TAG5,TAG6;

	MPI_Status status;	
	int i, j, k, l, n, nf1, nf2, const_n;
	double time=0.0; /*time1=0.0;*/

	
    nf1=(3*FDORDER/2)-1;
	nf2=nf1-1;
	int idx_s;
	int ll = 1;
	if(ABS_TYPE==1 && FDORDER==2){ll=2;}
	int nrow = NY + 2 * (ll*FDORDER/2) + 1, ncol = NX + 2 * (ll*FDORDER/2), ndep = NZ + 2 * (ll*FDORDER/2);
	int size = nrow * ncol * ndep;
	float *Ftb = transform3(buffertop_to_bot, 1, NX, 1, NZ,  1, nf2);
	float *Fbt = transform3(bufferbot_to_top, 1, NX, 1, NZ,  1, nf1);
	float *Flr = transform3(bufferlef_to_rig, 1, NY, 1, NZ,  1, nf2);
	float *Frl = transform3(bufferrig_to_lef, 1, NY, 1, NZ,  1, nf1);
	float *Ffb = transform3(bufferfro_to_bac, 1, NY, 1, NX,  1, nf2);
	float *Fbf = transform3(bufferbac_to_fro, 1, NY, 1, NX,  1, nf1);
    int strip = ndep;
    int slice = ncol * ndep;
    int strip_tb = nf2;
    int slice_tb = NZ * nf2;
    int strip_bt = nf1;
    int slice_bt = NZ * nf1;
    int strip_lr = nf2;
    int slice_lr = NZ * nf2;
    int strip_rl = nf1;
    int slice_rl = NZ * nf1;
    int strip_fb = nf2;
    int slice_fb = NX * nf2;
    int strip_bf = nf1;
    int slice_bf = NX * nf1;


	/*if (LOG)
	if (MYID==0) time1=MPI_Wtime();*/

	/* top-bottom -----------------------------------------------------------*/	
	

	if (POS[2]!=0) {	// no boundary exchange at top of global grid
		int n = 1;
        for (l=1;l<=FDORDER/2-1;l++) {
		    for (i=1;i<=NX;i++) {
			    // storage of top of local volume into buffer
			    for (k=1;k<=NZ;k++) {
			    	idx_s =  l * slice + i * strip + k;
			    	Ftb[i * slice_tb + k * strip_tb + n + 0] = Fs[6 * idx_s + 1];
			    	Ftb[i * slice_tb + k * strip_tb + n + 1] = Fs[6 * idx_s + 3];
			    	Ftb[i * slice_tb + k * strip_tb + n + 2] = Fs[6 * idx_s + 4];
			    }
		    }
		    n = n + 3;
		}
		for (l = FDORDER/2; l<=(FDORDER/2);l++) {
		    for (i=1;i<=NX;i++) {
			    // storage of top of local volume into buffer
			    for (k=1;k<=NZ;k++) {
			    	idx_s =  l * slice + i * strip + k;
			    	Ftb[i * slice_tb + k * strip_tb + n + 0] = Fs[6 * idx_s + 1];
			    }
		    }
		    n = n + 1;
		}
	}



    if (POS[2] != NPROCY - 1) {	// no boundary exchange at top of global grid
		int n = 1;
		for (l=FDORDER/2-1;l>=1;l--) {
		    for (i=1;i<=NX;i++) {
				// storage of top of local volume into buffer
				for (k=1;k<=NZ;k++) {
			    	idx_s =  l * slice + i * strip + k;
			    	Fbt[i * slice_bt + k * strip_bt + n + 0] = Fs[6 * (idx_s + (NY -2 * l  + 1) * slice) + 1];
			    	Fbt[i * slice_bt + k * strip_bt + n + 1] = Fs[6 * (idx_s + (NY -2 * l  + 1) * slice) + 3];
			    	Fbt[i * slice_bt + k * strip_bt + n + 2] = Fs[6 * (idx_s + (NY -2 * l  + 1) * slice) + 4];
				}
			}
		    n = n + 3;
		}
		for (l = FDORDER/2; l>=(FDORDER/2);l--) {
			for (i=1;i<=NX;i++) {
			    // storage of top of local volume into buffer
				for (k=1;k<=NZ;k++) {
			    	idx_s =  l * slice + i * strip + k;
			    	Fbt[i * slice_bt + k * strip_bt + n + 0] = Fs[6 * (idx_s + (NY -2 * l  + 1) * slice) + 3];
			    	Fbt[i * slice_bt + k * strip_bt + n + 1] = Fs[6 * (idx_s + (NY -2 * l  + 1) * slice) + 4];
				}
			}
			n = n + 2;
		}
	}


	/* persistent communication see comm_ini.c*/
	/*for (i=4;i<=5;i++){*/
		/* send and reveive values at edges of the local grid */
		/*MPI_Start(&req_send[i]);
		MPI_Wait(&req_send[i],&status);
		MPI_Start(&req_rec[i]);
		MPI_Wait(&req_rec[i],&status);
	}*/			
	
	
	MPI_Bsend(&Ftb[NZ * nf2 + nf2 + 1],NX*NZ*nf2,MPI_FLOAT,INDEX[3],TAG5,MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Recv(&Ftb[NZ * nf2 + nf2 + 1], NX*NZ*nf2,MPI_FLOAT,INDEX[4],TAG5,MPI_COMM_WORLD,&status);
	MPI_Bsend(&Fbt[NZ * nf1 + nf1 + 1],NX*NZ*nf1,MPI_FLOAT,INDEX[4],TAG6,MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Recv(&Fbt[NZ * nf1 + nf1 + 1], NX*NZ*nf1,MPI_FLOAT,INDEX[3],TAG6,MPI_COMM_WORLD,&status);
	



	if (POS[2]!=NPROCY-1) {	// no boundary exchange at top of global grid
		int n = 1;
        for (l=1;l<=FDORDER/2-1;l++) {
		    for (i=1;i<=NX;i++) {
			    // storage of top of local volume into buffer
			    for (k=1;k<=NZ;k++) {
			    	idx_s =  l * slice + i * strip + k;
			    	Fs[6 * (idx_s + NY * slice) + 1] = Ftb[i * slice_tb + k * strip_tb + n + 0];
			    	Fs[6 * (idx_s + NY * slice) + 3] = Ftb[i * slice_tb + k * strip_tb + n + 1];
			    	Fs[6 * (idx_s + NY * slice) + 4] = Ftb[i * slice_tb + k * strip_tb + n + 2];
			    }
		    }
		    n = n + 3;
		}
		for (l = FDORDER/2; l<=(FDORDER/2);l++) {
		    for (i=1;i<=NX;i++) {
			    // storage of top of local volume into buffer
			    for (k=1;k<=NZ;k++) {
			    	idx_s =  l * slice + i * strip + k;
			    	Fs[6 * (idx_s + NY * slice) + 1] = Ftb[i * slice_tb + k * strip_tb + n + 0];
			    }
		    }
		    n = n + 1;
		}

	}




	if (POS[2]!=0) {	// no boundary exchange at top of global grid
		int n = 1;
        for (l=FDORDER/2-1;l>=1;l--) {
		    for (i=1;i<=NX;i++) {
			    // storage of top of local volume into buffer
			    for (k=1;k<=NZ;k++) {
			    	idx_s =  l * slice + i * strip + k;
			    	Fs[6 * (idx_s + (1 -2 * l) * slice) + 1] = Fbt[i * slice_bt + k * strip_bt + n + 0];
			    	Fs[6 * (idx_s + (1 -2 * l) * slice) + 3] = Fbt[i * slice_bt + k * strip_bt + n + 1];
			    	Fs[6 * (idx_s + (1 -2 * l) * slice) + 4] = Fbt[i * slice_bt + k * strip_bt + n + 2];
			    }
		    }
		    n = n + 3;
		}
		for (l = FDORDER/2; l>=(FDORDER/2);l--) {
		    for (i=1;i<=NX;i++) {
			    // storage of top of local volume into buffer
			    for (k=1;k<=NZ;k++) {
			    	idx_s =  l * slice + i * strip + k;
			    	Fs[6 * (idx_s + (1 -2 * l) * slice) + 3] = Fbt[i * slice_bt + k * strip_bt + n + 0];
			    	Fs[6 * (idx_s + (1 -2 * l) * slice) + 4] = Fbt[i * slice_bt + k * strip_bt + n + 1];
			    }
		    }
		    n = n + 2;
		}

	}



	

	
	/* left-right -----------------------------------------------------------*/	
	
	


	if ((BOUNDARY) || (POS[1]!=0)) {	// no boundary exchange at left edge of global grid
		int n = 1;
	    for (j=1;j<=NY;j++){
	    	n = 1;
		    for (l=1;l<=(FDORDER/2-1);l++){
			    // storage of left edge of local volume into buffer
			    for (k=1;k<=NZ;k++) {
			    	idx_s =  j * slice + l * strip + k;
			    	Flr[j * slice_lr + k * strip_lr + n + 0] = Fs[6 * idx_s + 0];
			    	Flr[j * slice_lr + k * strip_lr + n + 1] = Fs[6 * idx_s + 3];
			    	Flr[j * slice_lr + k * strip_lr + n + 2] = Fs[6 * idx_s + 5];
			    }

			    n += 3;
		    }
	    }
	    const_n = n;
	    for (j=1;j<=NY;j++){
	    	n = const_n;
		    for (l=(FDORDER/2);l<=FDORDER/2;l++){
			    // storage of left edge of local volume into buffer
			    for (k=1;k<=NZ;k++) {
			    	idx_s =  j * slice + l * strip + k;
			    	Flr[j * slice_lr + k * strip_lr + n + 0] = Fs[6 * idx_s + 0];
			    }

			    n += 1;
		    }
	    }
	}


	if ((BOUNDARY) || (POS[1]!=NPROCX-1)) {	// no boundary exchange at left edge of global grid
		int n = 1;
	    for (j=1;j<=NY;j++){
	    	n = 1;
		    for (l=(FDORDER/2-1);l>=1;l--) {
			    // storage of left edge of local volume into buffer
			    for (k=1;k<=NZ;k++) {
			    	idx_s =  j * slice + l * strip + k;
			    	Frl[j * slice_rl + k * strip_rl + n + 0] = Fs[6 * (idx_s + (NX -2 * l  + 1) * strip) + 0];
			    	Frl[j * slice_rl + k * strip_rl + n + 1] = Fs[6 * (idx_s + (NX -2 * l  + 1) * strip) + 3];
			    	Frl[j * slice_rl + k * strip_rl + n + 2] = Fs[6 * (idx_s + (NX -2 * l  + 1) * strip) + 5];
			    }

			    n += 3;
		    }
	    }
	    const_n = n;
	    for (j=1;j<=NY;j++){
	    	n = const_n;
		    for (l=FDORDER/2;l>=(FDORDER/2);l--){
			    // storage of left edge of local volume into buffer
			    for (k=1;k<=NZ;k++) {
			    	idx_s =  j * slice + l * strip + k;
			    	Frl[j * slice_rl + k * strip_rl + n + 0] = Fs[6 * (idx_s + (NX -2 * l  + 1) * strip) + 3];
			    	Frl[j * slice_rl + k * strip_rl + n + 1] = Fs[6 * (idx_s + (NX -2 * l  + 1) * strip) + 5];
			    }

			    n += 1;
		    }
	    }
	}


	/* persistent communication see comm_ini.c*/
	/*for (i=0;i<=1;i++){*/
		/* send and reveive values at edges of the local grid */
		/*MPI_Start(&req_send[i]);
		MPI_Wait(&req_send[i],&status);
		MPI_Start(&req_rec[i]);
		MPI_Wait(&req_rec[i],&status);
	}*/			
	
 	
	
	MPI_Bsend(&Flr[NZ * nf2 + nf2 + 1],NY*NZ*nf2,MPI_FLOAT,INDEX[1],TAG1,MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Recv(&Flr[NZ * nf2 + nf2 + 1], NY*NZ*nf2,MPI_FLOAT,INDEX[2],TAG1,MPI_COMM_WORLD,&status);
	MPI_Bsend(&Frl[NZ * nf1 + nf1 + 1],NY*NZ*nf1,MPI_FLOAT,INDEX[2],TAG2,MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Recv(&Frl[NZ * nf1 + nf1 + 1], NY*NZ*nf1,MPI_FLOAT,INDEX[1],TAG2,MPI_COMM_WORLD,&status);
	




	if ((BOUNDARY) || (POS[1]!=NPROCX-1)) {	// no boundary exchange at left edge of global grid
		int n = 1;
	    for (j=1;j<=NY;j++){
	    	n = 1;
		    for (l=1;l<=(FDORDER/2-1);l++){
			    // storage of left edge of local volume into buffer
			    for (k=1;k<=NZ;k++) {
			    	idx_s =  j * slice + l * strip + k;
			    	Fs[6 * (idx_s + NX * strip) + 0] = Flr[j * slice_lr + k * strip_lr + n + 0];
			    	Fs[6 * (idx_s + NX * strip) + 3] = Flr[j * slice_lr + k * strip_lr + n + 1];
			    	Fs[6 * (idx_s + NX * strip) + 5] = Flr[j * slice_lr + k * strip_lr + n + 2];
			    }

			    n += 3;
		    }
	    }
	    const_n = n;
	    for (j=1;j<=NY;j++){
	    	n = const_n;
		    for (l=(FDORDER/2);l<=FDORDER/2;l++){
			    // storage of left edge of local volume into buffer
			    for (k=1;k<=NZ;k++) {
			    	idx_s =  j * slice + l * strip + k;
			    	Fs[6 * (idx_s + NX * strip) + 0] = Flr[j * slice_lr + k * strip_lr + n + 0];
			    }

			    n += 2;
		    }
	    }

	}

	
	if ((BOUNDARY) || (POS[1]!=0)) {	// no boundary exchange at left edge of global grid
		int n = 1;
	    for (j=1;j<=NY;j++){
	    	n = 1;
		    for (l=(FDORDER/2-1);l>=1;l--) {
			    // storage of left edge of local volume into buffer
			    for (k=1;k<=NZ;k++) {
			    	idx_s =  j * slice + l * strip + k;
			    	Fs[6 * (idx_s + (1 -2 * l) * strip) + 0] = Frl[j * slice_rl + k * strip_rl + n + 0];
			    	Fs[6 * (idx_s + (1 -2 * l) * strip) + 3] = Frl[j * slice_rl + k * strip_rl + n + 1];
			    	Fs[6 * (idx_s + (1 -2 * l) * strip) + 5] = Frl[j * slice_rl + k * strip_rl + n + 2];
			    }

			    n += 3;
		    }
	    }
	    const_n = n;
	    for (j=1;j<=NY;j++){
	    	n = const_n;
		    for (l=FDORDER/2;l>=(FDORDER/2);l--){
			    // storage of left edge of local volume into buffer
			    for (k=1;k<=NZ;k++) {
			    	idx_s =  j * slice + l * strip + k;
			    	Fs[6 * (idx_s + (1 -2 * l) * strip) + 3] = Frl[j * slice_rl + k * strip_rl + n + 0];
			    	Fs[6 * (idx_s + (1 -2 * l) * strip) + 5] = Frl[j * slice_rl + k * strip_rl + n + 1];
			    }

			    n += 1;
		    }
	    }
	}
	
	

	
	
	/* front-back -----------------------------------------------------------*/




	if ((BOUNDARY) || (POS[3]!=0))	//* no boundary exchange at front side of global grid
		for (j=1;j<=NY;j++) {
			for (i=1;i<=NX;i++) {   //* storage of front side of local volume into buffer
				n=1;
			    for (l=1;l<=FDORDER/2;l++){
			    	idx_s =  j * slice + i * strip + l;
			    	Ffb[j * slice_fb + i * strip_fb + (n++)] = Fs[6 * idx_s + 2];
			    }
			    for (l=1;l<=(FDORDER/2-1);l++) {
			    	idx_s =  j * slice + i * strip + l;
			    	Ffb[j * slice_fb + i * strip_fb + (n++)] = Fs[6 * idx_s + 4];
			    	Ffb[j * slice_fb + i * strip_fb + (n++)] = Fs[6 * idx_s + 5];
			    }
			}
		}


	if ((BOUNDARY) || (POS[3]!=NPROCZ-1))	//* no boundary exchange at back side of global grid
		for (j=1;j<=NY;j++) {
			for (i=1;i<=NX;i++) {
					//* storage of back side of local volume into buffer
				n=1;
				for (l=FDORDER/2;l>=1;l--) {
					idx_s =  j * slice + i * strip + l;
					Fbf[j * slice_bf + i * strip_bf + (n++)] = Fs[6 * (idx_s + (NZ -2 * l  + 1)) + 4];
					Fbf[j * slice_bf + i * strip_bf + (n++)] = Fs[6 * (idx_s + (NZ -2 * l  + 1)) + 5];
				}

				for (l=(FDORDER/2-1);l>=1;l--) {
					idx_s =  j * slice + i * strip + l;
                    Fbf[j * slice_bf + i * strip_bf + (n++)] = Fs[6 * (idx_s + (NZ -2 * l  + 1)) + 2];
				}

			}
		}






	/* persistent communication see comm_ini.c*/
	/*for (i=2;i<=3;i++){*/
		/* send and reveive values at edges of the local grid */
		/*MPI_Start(&req_send[i]);
		MPI_Wait(&req_send[i],&status);
		MPI_Start(&req_rec[i]);
		MPI_Wait(&req_rec[i],&status);
	}*/			

	
	MPI_Bsend(&Ffb[NX * nf2 + nf2 + 1],NX*NY*nf2,MPI_FLOAT,INDEX[5],TAG3,MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Recv(&Ffb[NX * nf2 + nf2 + 1], NX*NY*nf2,MPI_FLOAT,INDEX[6],TAG3,MPI_COMM_WORLD,&status);
	MPI_Bsend(&Fbf[NX * nf1 + nf1 + 1],NX*NY*nf1,MPI_FLOAT,INDEX[6],TAG4,MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Recv(&Fbf[NX * nf1 + nf1 + 1], NX*NY*nf1,MPI_FLOAT,INDEX[5],TAG4,MPI_COMM_WORLD,&status);




	if ((BOUNDARY) || (POS[3]!=NPROCZ-1))  //* no boundary exchange at back side of global grid
		for (j=1;j<=NY;j++) {
			for (i=1;i<=NX;i++) {
				n=1;
				for (l=1;l<=FDORDER/2;l++) {
					idx_s =  j * slice + i * strip + l;
					Fs[6 * (idx_s + NZ) + 2] = Ffb[j * slice_fb + i * strip_fb + (n++)];

				}

				for (l=1;l<=(FDORDER/2-1);l++) {
					idx_s =  j * slice + i * strip + l;
					Fs[6 * (idx_s + NZ) + 4] = Ffb[j * slice_fb + i * strip_fb + (n++)];
					Fs[6 * (idx_s + NZ) + 5] = Ffb[j * slice_fb + i * strip_fb + (n++)];
				}
			}
		}



	/* no exchange if periodic boundary condition is applied */
	if ((BOUNDARY) || (POS[3]!=0))	/* no boundary exchange at front side of global grid */
		for (j=1;j<=NY;j++) {
			for (i=1;i<=NX;i++) {
				n=1;
				for (l=FDORDER/2;l>=1;l--) {
					idx_s =  j * slice + i * strip + l;
					Fs[6 * (idx_s + (1 -2 * l)) + 4] = Fbf[j * slice_bf + i * strip_bf + (n++)];
					Fs[6 * (idx_s + (1 -2 * l)) + 5] = Fbf[j * slice_bf + i * strip_bf + (n++)];
				}

				for (l=(FDORDER/2-1);l>=1;l--) {
					idx_s =  j * slice + i * strip + l;
					Fs[6 * (idx_s + (1 -2 * l)) + 2] = Fbf[j * slice_bf + i * strip_bf + (n++)];
				}
			}
		}



	free_trans(Ftb, 1, NX, 1, NZ, 1, nf2);
	free_trans(Fbt, 1, NX, 1, NZ, 1, nf1);
	free_trans(Flr, 1, NY, 1, NZ, 1, nf2);
	free_trans(Frl, 1, NY, 1, NZ, 1, nf1);
	free_trans(Ffb, 1, NY, 1, NX, 1, nf2);
	free_trans(Fbf, 1, NY, 1, NX, 1, nf1);

	/*if (LOG)
	if (MYID==0){
		time2=MPI_Wtime();
		time=time2-time1;
		fprintf(FP," Real time for stress tensor exchange: \t\t %4.2f s.\n",time);
	}*/
	return time;

}



double exchange_s2_1(float *Fs, float *** buffertop_to_bot, MPI_Request * req_send, MPI_Request * req_rec) {
    
	extern int NX, NY, NZ, POS[4], NPROCX, NPROCY, NPROCZ, BOUNDARY, FDORDER,  INDEX[7], ABS_TYPE, MYID;    /*MYID,LOG,*/
	extern const int TAG1,TAG2,TAG3,TAG4,TAG5,TAG6;

	MPI_Status status;	
	int i, j, k, l, n, nf1, nf2, const_n;
	double time=0.0; /*time1=0.0;*/

	
    nf1=(3*FDORDER/2)-1;
	nf2=nf1-1;
	int idx_s;
	int ll = 1;
	if(ABS_TYPE==1 && FDORDER==2){ll=2;}
	int nrow = NY + 2 * (ll*FDORDER/2) + 1, ncol = NX + 2 * (ll*FDORDER/2), ndep = NZ + 2 * (ll*FDORDER/2);
	int size = nrow * ncol * ndep;
	float *Ftb = transform3(buffertop_to_bot, 1, NX, 1, NZ,  1, nf2);
	int ny = ((NY % 2 == 0) ? (NY/2):(NY/2+1));
	float Flr[(ny+1)*(NZ+1)*(nf2+1)],Frl[(ny+1)*(NZ+1)*(nf1+1)],Ffb[(ny+1)*(NX+1)*(nf2+1)],Fbf[(ny+1)*(NX+1)*(nf1+1)];
    int strip = ndep;
    int slice = ncol * ndep;
    int strip_tb = nf2;
    int slice_tb = NZ * nf2;
    int strip_bt = nf1;
    int slice_bt = NZ * nf1;
    int strip_lr = nf2;
    int slice_lr = NZ * nf2;
    int strip_rl = nf1;
    int slice_rl = NZ * nf1;
    int strip_fb = nf2;
    int slice_fb = NX * nf2;
    int strip_bf = nf1;
    int slice_bf = NX * nf1;


	/*if (LOG)
	if (MYID==0) time1=MPI_Wtime();*/

	/* top-bottom -----------------------------------------------------------*/	
	

	if (POS[2]!=0) {	// no boundary exchange at top of global grid
		int n = 1;
        for (l=1;l<=FDORDER/2-1;l++) {
		    for (i=1;i<=NX;i++) {
			    // storage of top of local volume into buffer
			    for (k=1;k<=NZ;k++) {
			    	idx_s =  l * slice + i * strip + k;
			    	Ftb[i * slice_tb + k * strip_tb + n + 0] = Fs[6 * idx_s + 1];
			    	Ftb[i * slice_tb + k * strip_tb + n + 1] = Fs[6 * idx_s + 3];
			    	Ftb[i * slice_tb + k * strip_tb + n + 2] = Fs[6 * idx_s + 4];
			    }
		    }
		    n = n + 3;
		}
		for (l = FDORDER/2; l<=(FDORDER/2);l++) {
		    for (i=1;i<=NX;i++) {
			    // storage of top of local volume into buffer
			    for (k=1;k<=NZ;k++) {
			    	idx_s =  l * slice + i * strip + k;
			    	Ftb[i * slice_tb + k * strip_tb + n + 0] = Fs[6 * idx_s + 1];
			    }
		    }
		    n = n + 1;
		}
	}



    // if (POS[2] != NPROCY - 1) {	// no boundary exchange at top of global grid
	// 	int n = 1;
	// 	for (l=FDORDER/2-1;l>=1;l--) {
	// 	    for (i=1;i<=NX;i++) {
	// 			// storage of top of local volume into buffer
	// 			for (k=1;k<=NZ;k++) {
	// 		    	idx_s =  l * slice + i * strip + k;
	// 		    	Fbt[i * slice_bt + k * strip_bt + n + 0] = Fs[6 * (idx_s + (NY -2 * l  + 1) * slice) + 1];
	// 		    	Fbt[i * slice_bt + k * strip_bt + n + 1] = Fs[6 * (idx_s + (NY -2 * l  + 1) * slice) + 3];
	// 		    	Fbt[i * slice_bt + k * strip_bt + n + 2] = Fs[6 * (idx_s + (NY -2 * l  + 1) * slice) + 4];
	// 			}
	// 		}
	// 	    n = n + 3;
	// 	}
	// 	for (l = FDORDER/2; l>=(FDORDER/2);l--) {
	// 		for (i=1;i<=NX;i++) {
	// 		    // storage of top of local volume into buffer
	// 			for (k=1;k<=NZ;k++) {
	// 		    	idx_s =  l * slice + i * strip + k;
	// 		    	Fbt[i * slice_bt + k * strip_bt + n + 0] = Fs[6 * (idx_s + (NY -2 * l  + 1) * slice) + 3];
	// 		    	Fbt[i * slice_bt + k * strip_bt + n + 1] = Fs[6 * (idx_s + (NY -2 * l  + 1) * slice) + 4];
	// 			}
	// 		}
	// 		n = n + 2;
	// 	}
	// }


	/* persistent communication see comm_ini.c*/
	/*for (i=4;i<=5;i++){*/
		/* send and reveive values at edges of the local grid */
		/*MPI_Start(&req_send[i]);
		MPI_Wait(&req_send[i],&status);
		MPI_Start(&req_rec[i]);
		MPI_Wait(&req_rec[i],&status);
	}*/			
	
	
	MPI_Bsend(&Ftb[NZ * nf2 + nf2 + 1],NX*NZ*nf2,MPI_FLOAT,INDEX[3],TAG5,MPI_COMM_WORLD);
	


	// if (POS[2]!=0) {	// no boundary exchange at top of global grid
	// 	int n = 1;
    //     for (l=FDORDER/2-1;l>=1;l--) {
	// 	    for (i=1;i<=NX;i++) {
	// 		    // storage of top of local volume into buffer
	// 		    for (k=1;k<=NZ;k++) {
	// 		    	idx_s =  l * slice + i * strip + k;
	// 		    	Fs[6 * (idx_s + (1 -2 * l) * slice) + 1] = Fbt[i * slice_bt + k * strip_bt + n + 0];
	// 		    	Fs[6 * (idx_s + (1 -2 * l) * slice) + 3] = Fbt[i * slice_bt + k * strip_bt + n + 1];
	// 		    	Fs[6 * (idx_s + (1 -2 * l) * slice) + 4] = Fbt[i * slice_bt + k * strip_bt + n + 2];
	// 		    }
	// 	    }
	// 	    n = n + 3;
	// 	}
	// 	for (l = FDORDER/2; l>=(FDORDER/2);l--) {
	// 	    for (i=1;i<=NX;i++) {
	// 		    // storage of top of local volume into buffer
	// 		    for (k=1;k<=NZ;k++) {
	// 		    	idx_s =  l * slice + i * strip + k;
	// 		    	Fs[6 * (idx_s + (1 -2 * l) * slice) + 3] = Fbt[i * slice_bt + k * strip_bt + n + 0];
	// 		    	Fs[6 * (idx_s + (1 -2 * l) * slice) + 4] = Fbt[i * slice_bt + k * strip_bt + n + 1];
	// 		    }
	// 	    }
	// 	    n = n + 2;
	// 	}

	// }



	

	
	/* left-right -----------------------------------------------------------*/	
	
	


	if ((BOUNDARY) || (POS[1]!=0)) {	// no boundary exchange at left edge of global grid
		int n = 1;
	    for (j=1;j<=ny;j++){
	    	n = 1;
		    for (l=1;l<=(FDORDER/2-1);l++){
			    // storage of left edge of local volume into buffer
			    for (k=1;k<=NZ;k++) {
			    	idx_s =  j * slice + l * strip + k;
			    	Flr[j * slice_lr + k * strip_lr + n + 0] = Fs[6 * idx_s + 0];
			    	Flr[j * slice_lr + k * strip_lr + n + 1] = Fs[6 * idx_s + 3];
			    	Flr[j * slice_lr + k * strip_lr + n + 2] = Fs[6 * idx_s + 5];
			    }

			    n += 3;
		    }
	    }
	    const_n = n;
	    for (j=1;j<=ny;j++){
	    	n = const_n;
		    for (l=(FDORDER/2);l<=FDORDER/2;l++){
			    // storage of left edge of local volume into buffer
			    for (k=1;k<=NZ;k++) {
			    	idx_s =  j * slice + l * strip + k;
			    	Flr[j * slice_lr + k * strip_lr + n + 0] = Fs[6 * idx_s + 0];
			    }

			    n += 1;
		    }
	    }
	}


	if ((BOUNDARY) || (POS[1]!=NPROCX-1)) {	// no boundary exchange at left edge of global grid
		int n = 1;
	    for (j=1;j<=ny;j++){
	    	n = 1;
		    for (l=(FDORDER/2-1);l>=1;l--) {
			    // storage of left edge of local volume into buffer
			    for (k=1;k<=NZ;k++) {
			    	idx_s =  j * slice + l * strip + k;
			    	Frl[j * slice_rl + k * strip_rl + n + 0] = Fs[6 * (idx_s + (NX -2 * l  + 1) * strip) + 0];
			    	Frl[j * slice_rl + k * strip_rl + n + 1] = Fs[6 * (idx_s + (NX -2 * l  + 1) * strip) + 3];
			    	Frl[j * slice_rl + k * strip_rl + n + 2] = Fs[6 * (idx_s + (NX -2 * l  + 1) * strip) + 5];
			    }

			    n += 3;
		    }
	    }
	    const_n = n;
	    for (j=1;j<=ny;j++){
	    	n = const_n;
		    for (l=FDORDER/2;l>=(FDORDER/2);l--){
			    // storage of left edge of local volume into buffer
			    for (k=1;k<=NZ;k++) {
			    	idx_s =  j * slice + l * strip + k;
			    	Frl[j * slice_rl + k * strip_rl + n + 0] = Fs[6 * (idx_s + (NX -2 * l  + 1) * strip) + 3];
			    	Frl[j * slice_rl + k * strip_rl + n + 1] = Fs[6 * (idx_s + (NX -2 * l  + 1) * strip) + 5];
			    }

			    n += 1;
		    }
	    }
	}


	/* persistent communication see comm_ini.c*/
	/*for (i=0;i<=1;i++){*/
		/* send and reveive values at edges of the local grid */
		/*MPI_Start(&req_send[i]);
		MPI_Wait(&req_send[i],&status);
		MPI_Start(&req_rec[i]);
		MPI_Wait(&req_rec[i],&status);
	}*/			
	
 	
	
	MPI_Bsend(&Flr[NZ * nf2 + nf2 + 1],ny*NZ*nf2,MPI_FLOAT,INDEX[1],TAG1,MPI_COMM_WORLD);
	MPI_Bsend(&Frl[NZ * nf1 + nf1 + 1],ny*NZ*nf1,MPI_FLOAT,INDEX[2],TAG2,MPI_COMM_WORLD);
	
	
	

	
	
	/* front-back -----------------------------------------------------------*/




	if ((BOUNDARY) || (POS[3]!=0))	//* no boundary exchange at front side of global grid
		for (j=1;j<=ny;j++) {
			for (i=1;i<=NX;i++) {   //* storage of front side of local volume into buffer
				n=1;
			    for (l=1;l<=FDORDER/2;l++){
			    	idx_s =  j * slice + i * strip + l;
			    	Ffb[j * slice_fb + i * strip_fb + (n++)] = Fs[6 * idx_s + 2];
			    }
			    for (l=1;l<=(FDORDER/2-1);l++) {
			    	idx_s =  j * slice + i * strip + l;
			    	Ffb[j * slice_fb + i * strip_fb + (n++)] = Fs[6 * idx_s + 4];
			    	Ffb[j * slice_fb + i * strip_fb + (n++)] = Fs[6 * idx_s + 5];
			    }
			}
		}


	if ((BOUNDARY) || (POS[3]!=NPROCZ-1))	//* no boundary exchange at back side of global grid
		for (j=1;j<=ny;j++) {
			for (i=1;i<=NX;i++) {
					//* storage of back side of local volume into buffer
				n=1;
				for (l=FDORDER/2;l>=1;l--) {
					idx_s =  j * slice + i * strip + l;
					Fbf[j * slice_bf + i * strip_bf + (n++)] = Fs[6 * (idx_s + (NZ -2 * l  + 1)) + 4];
					Fbf[j * slice_bf + i * strip_bf + (n++)] = Fs[6 * (idx_s + (NZ -2 * l  + 1)) + 5];
				}

				for (l=(FDORDER/2-1);l>=1;l--) {
					idx_s =  j * slice + i * strip + l;
                    Fbf[j * slice_bf + i * strip_bf + (n++)] = Fs[6 * (idx_s + (NZ -2 * l  + 1)) + 2];
				}

			}
		}






	/* persistent communication see comm_ini.c*/
	/*for (i=2;i<=3;i++){*/
		/* send and reveive values at edges of the local grid */
		/*MPI_Start(&req_send[i]);
		MPI_Wait(&req_send[i],&status);
		MPI_Start(&req_rec[i]);
		MPI_Wait(&req_rec[i],&status);
	}*/			

	
	MPI_Bsend(&Ffb[NX * nf2 + nf2 + 1],NX*ny*nf2,MPI_FLOAT,INDEX[5],TAG3,MPI_COMM_WORLD);
	MPI_Bsend(&Fbf[NX * nf1 + nf1 + 1],NX*ny*nf1,MPI_FLOAT,INDEX[6],TAG4,MPI_COMM_WORLD);
    

	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Recv(&Ftb[NZ * nf2 + nf2 + 1], NX*NZ*nf2,MPI_FLOAT,INDEX[4],TAG5,MPI_COMM_WORLD,&status);

	if (POS[2]!=NPROCY-1) {	// no boundary exchange at top of global grid
		int n = 1;
        for (l=1;l<=FDORDER/2-1;l++) {
		    for (i=1;i<=NX;i++) {
			    // storage of top of local volume into buffer
			    for (k=1;k<=NZ;k++) {
			    	idx_s =  l * slice + i * strip + k;
			    	Fs[6 * (idx_s + NY * slice) + 1] = Ftb[i * slice_tb + k * strip_tb + n + 0];
			    	Fs[6 * (idx_s + NY * slice) + 3] = Ftb[i * slice_tb + k * strip_tb + n + 1];
			    	Fs[6 * (idx_s + NY * slice) + 4] = Ftb[i * slice_tb + k * strip_tb + n + 2];
			    }
		    }
		    n = n + 3;
		}
		for (l = FDORDER/2; l<=(FDORDER/2);l++) {
		    for (i=1;i<=NX;i++) {
			    // storage of top of local volume into buffer
			    for (k=1;k<=NZ;k++) {
			    	idx_s =  l * slice + i * strip + k;
			    	Fs[6 * (idx_s + NY * slice) + 1] = Ftb[i * slice_tb + k * strip_tb + n + 0];
			    }
		    }
		    n = n + 1;
		}

	}

	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Recv(&Flr[NZ * nf2 + nf2 + 1], ny*NZ*nf2,MPI_FLOAT,INDEX[2],TAG1,MPI_COMM_WORLD,&status);
	MPI_Recv(&Frl[NZ * nf1 + nf1 + 1], ny*NZ*nf1,MPI_FLOAT,INDEX[1],TAG2,MPI_COMM_WORLD,&status);
	




	if ((BOUNDARY) || (POS[1]!=NPROCX-1)) {	// no boundary exchange at left edge of global grid
		int n = 1;
	    for (j=1;j<=ny;j++){
	    	n = 1;
		    for (l=1;l<=(FDORDER/2-1);l++){
			    // storage of left edge of local volume into buffer
			    for (k=1;k<=NZ;k++) {
			    	idx_s =  j * slice + l * strip + k;
			    	Fs[6 * (idx_s + NX * strip) + 0] = Flr[j * slice_lr + k * strip_lr + n + 0];
			    	Fs[6 * (idx_s + NX * strip) + 3] = Flr[j * slice_lr + k * strip_lr + n + 1];
			    	Fs[6 * (idx_s + NX * strip) + 5] = Flr[j * slice_lr + k * strip_lr + n + 2];
			    }

			    n += 3;
		    }
	    }
	    const_n = n;
	    for (j=1;j<=ny;j++){
	    	n = const_n;
		    for (l=(FDORDER/2);l<=FDORDER/2;l++){
			    // storage of left edge of local volume into buffer
			    for (k=1;k<=NZ;k++) {
			    	idx_s =  j * slice + l * strip + k;
			    	Fs[6 * (idx_s + NX * strip) + 0] = Flr[j * slice_lr + k * strip_lr + n + 0];
			    }

			    n += 2;
		    }
	    }

	}

	
	if ((BOUNDARY) || (POS[1]!=0)) {	// no boundary exchange at left edge of global grid
		int n = 1;
	    for (j=1;j<=ny;j++){
	    	n = 1;
		    for (l=(FDORDER/2-1);l>=1;l--) {
			    // storage of left edge of local volume into buffer
			    for (k=1;k<=NZ;k++) {
			    	idx_s =  j * slice + l * strip + k;
			    	Fs[6 * (idx_s + (1 -2 * l) * strip) + 0] = Frl[j * slice_rl + k * strip_rl + n + 0];
			    	Fs[6 * (idx_s + (1 -2 * l) * strip) + 3] = Frl[j * slice_rl + k * strip_rl + n + 1];
			    	Fs[6 * (idx_s + (1 -2 * l) * strip) + 5] = Frl[j * slice_rl + k * strip_rl + n + 2];
			    }

			    n += 3;
		    }
	    }
	    const_n = n;
	    for (j=1;j<=ny;j++){
	    	n = const_n;
		    for (l=FDORDER/2;l>=(FDORDER/2);l--){
			    // storage of left edge of local volume into buffer
			    for (k=1;k<=NZ;k++) {
			    	idx_s =  j * slice + l * strip + k;
			    	Fs[6 * (idx_s + (1 -2 * l) * strip) + 3] = Frl[j * slice_rl + k * strip_rl + n + 0];
			    	Fs[6 * (idx_s + (1 -2 * l) * strip) + 5] = Frl[j * slice_rl + k * strip_rl + n + 1];
			    }

			    n += 1;
		    }
	    }
	}



    MPI_Barrier(MPI_COMM_WORLD);
	MPI_Recv(&Ffb[NX * nf2 + nf2 + 1], NX*ny*nf2,MPI_FLOAT,INDEX[6],TAG3,MPI_COMM_WORLD,&status);
	MPI_Recv(&Fbf[NX * nf1 + nf1 + 1], NX*ny*nf1,MPI_FLOAT,INDEX[5],TAG4,MPI_COMM_WORLD,&status);

	if ((BOUNDARY) || (POS[3]!=NPROCZ-1))  //* no boundary exchange at back side of global grid
		for (j=1;j<=ny;j++) {
			for (i=1;i<=NX;i++) {
				n=1;
				for (l=1;l<=FDORDER/2;l++) {
					idx_s =  j * slice + i * strip + l;
					Fs[6 * (idx_s + NZ) + 2] = Ffb[j * slice_fb + i * strip_fb + (n++)];

				}

				for (l=1;l<=(FDORDER/2-1);l++) {
					idx_s =  j * slice + i * strip + l;
					Fs[6 * (idx_s + NZ) + 4] = Ffb[j * slice_fb + i * strip_fb + (n++)];
					Fs[6 * (idx_s + NZ) + 5] = Ffb[j * slice_fb + i * strip_fb + (n++)];
				}
			}
		}



	/* no exchange if periodic boundary condition is applied */
	if ((BOUNDARY) || (POS[3]!=0))	/* no boundary exchange at front side of global grid */
		for (j=1;j<=ny;j++) {
			for (i=1;i<=NX;i++) {
				n=1;
				for (l=FDORDER/2;l>=1;l--) {
					idx_s =  j * slice + i * strip + l;
					Fs[6 * (idx_s + (1 -2 * l)) + 4] = Fbf[j * slice_bf + i * strip_bf + (n++)];
					Fs[6 * (idx_s + (1 -2 * l)) + 5] = Fbf[j * slice_bf + i * strip_bf + (n++)];
				}

				for (l=(FDORDER/2-1);l>=1;l--) {
					idx_s =  j * slice + i * strip + l;
					Fs[6 * (idx_s + (1 -2 * l)) + 2] = Fbf[j * slice_bf + i * strip_bf + (n++)];
				}
			}
		}



	free_trans(Ftb, 1, NX, 1, NZ, 1, nf2);

	/*if (LOG)
	if (MYID==0){
		time2=MPI_Wtime();
		time=time2-time1;
		fprintf(FP," Real time for stress tensor exchange: \t\t %4.2f s.\n",time);
	}*/
	return time;
	
}

double exchange_s2_2(float *Fs, float *** bufferbot_to_top, MPI_Request * req_send, MPI_Request * req_rec) {
    
	extern int NX, NY, NZ, POS[4], NPROCX, NPROCY, NPROCZ, BOUNDARY, FDORDER,  INDEX[7], ABS_TYPE, MYID;    /*MYID,LOG,*/
	extern const int TAG1,TAG2,TAG3,TAG4,TAG5,TAG6;

	MPI_Status status;	
	int i, j, k, l, n, nf1, nf2, const_n;
	double time=0.0; /*time1=0.0;*/

	
    nf1=(3*FDORDER/2)-1;
	nf2=nf1-1;
	int idx_s;
	int ll = 1;
	if(ABS_TYPE==1 && FDORDER==2){ll=2;}
	int nrow = NY + 2 * (ll*FDORDER/2) + 1, ncol = NX + 2 * (ll*FDORDER/2), ndep = NZ + 2 * (ll*FDORDER/2);
	int size = nrow * ncol * ndep;
	float *Fbt = transform3(bufferbot_to_top, 1, NX, 1, NZ,  1, nf1);
	int ny = NY/2;
	int ny1 = ((NY % 2 == 0) ? (NY/2):(NY/2+1));
	float Flr[(ny+1)*(NZ+1)*(nf2+1)],Frl[(ny+1)*(NZ+1)*(nf1+1)],Ffb[(ny+1)*(NX+1)*(nf2+1)],Fbf[(ny+1)*(NX+1)*(nf1+1)];
    int strip = ndep;
    int slice = ncol * ndep;
    int strip_tb = nf2;
    int slice_tb = NZ * nf2;
    int strip_bt = nf1;
    int slice_bt = NZ * nf1;
    int strip_lr = nf2;
    int slice_lr = NZ * nf2;
    int strip_rl = nf1;
    int slice_rl = NZ * nf1;
    int strip_fb = nf2;
    int slice_fb = NX * nf2;
    int strip_bf = nf1;
    int slice_bf = NX * nf1;


	/*if (LOG)
	if (MYID==0) time1=MPI_Wtime();*/

	/* top-bottom -----------------------------------------------------------*/	
	

	// if (POS[2]!=0) {	// no boundary exchange at top of global grid
	// 	int n = 1;
    //     for (l=1;l<=FDORDER/2-1;l++) {
	// 	    for (i=1;i<=NX;i++) {
	// 		    // storage of top of local volume into buffer
	// 		    for (k=1;k<=NZ;k++) {
	// 		    	idx_s =  l * slice + i * strip + k;
	// 		    	Ftb[i * slice_tb + k * strip_tb + n + 0] = Fs[6 * idx_s + 1];
	// 		    	Ftb[i * slice_tb + k * strip_tb + n + 1] = Fs[6 * idx_s + 3];
	// 		    	Ftb[i * slice_tb + k * strip_tb + n + 2] = Fs[6 * idx_s + 4];
	// 		    }
	// 	    }
	// 	    n = n + 3;
	// 	}
	// 	for (l = FDORDER/2; l<=(FDORDER/2);l++) {
	// 	    for (i=1;i<=NX;i++) {
	// 		    // storage of top of local volume into buffer
	// 		    for (k=1;k<=NZ;k++) {
	// 		    	idx_s =  l * slice + i * strip + k;
	// 		    	Ftb[i * slice_tb + k * strip_tb + n + 0] = Fs[6 * idx_s + 1];
	// 		    }
	// 	    }
	// 	    n = n + 1;
	// 	}
	// }



    if (POS[2] != NPROCY - 1) {	// no boundary exchange at top of global grid
		int n = 1;
		for (l=FDORDER/2-1;l>=1;l--) {
		    for (i=1;i<=NX;i++) {
				// storage of top of local volume into buffer
				for (k=1;k<=NZ;k++) {
			    	idx_s =  l * slice + i * strip + k;
			    	Fbt[i * slice_bt + k * strip_bt + n + 0] = Fs[6 * (idx_s + (NY -2 * l  + 1) * slice) + 1];
			    	Fbt[i * slice_bt + k * strip_bt + n + 1] = Fs[6 * (idx_s + (NY -2 * l  + 1) * slice) + 3];
			    	Fbt[i * slice_bt + k * strip_bt + n + 2] = Fs[6 * (idx_s + (NY -2 * l  + 1) * slice) + 4];
				}
			}
		    n = n + 3;
		}
		for (l = FDORDER/2; l>=(FDORDER/2);l--) {
			for (i=1;i<=NX;i++) {
			    // storage of top of local volume into buffer
				for (k=1;k<=NZ;k++) {
			    	idx_s =  l * slice + i * strip + k;
			    	Fbt[i * slice_bt + k * strip_bt + n + 0] = Fs[6 * (idx_s + (NY -2 * l  + 1) * slice) + 3];
			    	Fbt[i * slice_bt + k * strip_bt + n + 1] = Fs[6 * (idx_s + (NY -2 * l  + 1) * slice) + 4];
				}
			}
			n = n + 2;
		}
	}


	/* persistent communication see comm_ini.c*/
	/*for (i=4;i<=5;i++){*/
		/* send and reveive values at edges of the local grid */
		/*MPI_Start(&req_send[i]);
		MPI_Wait(&req_send[i],&status);
		MPI_Start(&req_rec[i]);
		MPI_Wait(&req_rec[i],&status);
	}*/			
	
	
	// MPI_Bsend(&Ftb[NZ * nf2 + nf2 + 1],NX*NZ*nf2,MPI_FLOAT,INDEX[3],TAG5,MPI_COMM_WORLD);
	// MPI_Barrier(MPI_COMM_WORLD);
	// MPI_Recv(&Ftb[NZ * nf2 + nf2 + 1], NX*NZ*nf2,MPI_FLOAT,INDEX[4],TAG5,MPI_COMM_WORLD,&status);
	MPI_Bsend(&Fbt[NZ * nf1 + nf1 + 1],NX*NZ*nf1,MPI_FLOAT,INDEX[4],TAG6,MPI_COMM_WORLD);



	

	
	/* left-right -----------------------------------------------------------*/	
	
	


	if ((BOUNDARY) || (POS[1]!=0)) {	// no boundary exchange at left edge of global grid
		int n = 1;
	    for (j=ny1+1;j<=NY;j++){
	    	n = 1;
		    for (l=1;l<=(FDORDER/2-1);l++){
			    // storage of left edge of local volume into buffer
			    for (k=1;k<=NZ;k++) {
			    	idx_s =  j * slice + l * strip + k;
			    	Flr[(j - ny1) * slice_lr + k * strip_lr + n + 0] = Fs[6 * idx_s + 0];
			    	Flr[(j - ny1) * slice_lr + k * strip_lr + n + 1] = Fs[6 * idx_s + 3];
			    	Flr[(j - ny1) * slice_lr + k * strip_lr + n + 2] = Fs[6 * idx_s + 5];
			    }

			    n += 3;
		    }
	    }
	    const_n = n;
	    for (j=ny1+1;j<=NY;j++){
	    	n = const_n;
		    for (l=(FDORDER/2);l<=FDORDER/2;l++){
			    // storage of left edge of local volume into buffer
			    for (k=1;k<=NZ;k++) {
			    	idx_s =  j * slice + l * strip + k;
			    	Flr[(j - ny1) * slice_lr + k * strip_lr + n + 0] = Fs[6 * idx_s + 0];
			    }

			    n += 1;
		    }
	    }
	}


	if ((BOUNDARY) || (POS[1]!=NPROCX-1)) {	// no boundary exchange at left edge of global grid
		int n = 1;
	    for (j=ny1+1;j<=NY;j++){
	    	n = 1;
		    for (l=(FDORDER/2-1);l>=1;l--) {
			    // storage of left edge of local volume into buffer
			    for (k=1;k<=NZ;k++) {
			    	idx_s =  j * slice + l * strip + k;
			    	Frl[(j - ny1) * slice_rl + k * strip_rl + n + 0] = Fs[6 * (idx_s + (NX -2 * l  + 1) * strip) + 0];
			    	Frl[(j - ny1) * slice_rl + k * strip_rl + n + 1] = Fs[6 * (idx_s + (NX -2 * l  + 1) * strip) + 3];
			    	Frl[(j - ny1) * slice_rl + k * strip_rl + n + 2] = Fs[6 * (idx_s + (NX -2 * l  + 1) * strip) + 5];
			    }

			    n += 3;
		    }
	    }
	    const_n = n;
	    for (j=ny1+1;j<=NY;j++){
	    	n = const_n;
		    for (l=FDORDER/2;l>=(FDORDER/2);l--){
			    // storage of left edge of local volume into buffer
			    for (k=1;k<=NZ;k++) {
			    	idx_s =  j * slice + l * strip + k;
			    	Frl[(j - ny1) * slice_rl + k * strip_rl + n + 0] = Fs[6 * (idx_s + (NX -2 * l  + 1) * strip) + 3];
			    	Frl[(j - ny1) * slice_rl + k * strip_rl + n + 1] = Fs[6 * (idx_s + (NX -2 * l  + 1) * strip) + 5];
			    }

			    n += 1;
		    }
	    }
	}


	/* persistent communication see comm_ini.c*/
	/*for (i=0;i<=1;i++){*/
		/* send and reveive values at edges of the local grid */
		/*MPI_Start(&req_send[i]);
		MPI_Wait(&req_send[i],&status);
		MPI_Start(&req_rec[i]);
		MPI_Wait(&req_rec[i],&status);
	}*/			
	
 	
	
	MPI_Bsend(&Flr[NZ * nf2 + nf2 + 1],ny*NZ*nf2,MPI_FLOAT,INDEX[1],TAG1,MPI_COMM_WORLD);
	MPI_Bsend(&Frl[NZ * nf1 + nf1 + 1],ny*NZ*nf1,MPI_FLOAT,INDEX[2],TAG2,MPI_COMM_WORLD);
	
	
	/* front-back -----------------------------------------------------------*/




	if ((BOUNDARY) || (POS[3]!=0))	//* no boundary exchange at front side of global grid
		for (j=ny1+1;j<=NY;j++) {
			for (i=1;i<=NX;i++) {   //* storage of front side of local volume into buffer
				n=1;
			    for (l=1;l<=FDORDER/2;l++){
			    	idx_s =  j * slice + i * strip + l;
			    	Ffb[(j - ny1) * slice_fb + i * strip_fb + (n++)] = Fs[6 * idx_s + 2];
			    }
			    for (l=1;l<=(FDORDER/2-1);l++) {
			    	idx_s =  j * slice + i * strip + l;
			    	Ffb[(j - ny1) * slice_fb + i * strip_fb + (n++)] = Fs[6 * idx_s + 4];
			    	Ffb[(j - ny1) * slice_fb + i * strip_fb + (n++)] = Fs[6 * idx_s + 5];
			    }
			}
		}


	if ((BOUNDARY) || (POS[3]!=NPROCZ-1))	//* no boundary exchange at back side of global grid
		for (j=ny1+1;j<=NY;j++) {
			for (i=1;i<=NX;i++) {
					//* storage of back side of local volume into buffer
				n=1;
				for (l=FDORDER/2;l>=1;l--) {
					idx_s =  j * slice + i * strip + l;
					Fbf[(j - ny1) * slice_bf + i * strip_bf + (n++)] = Fs[6 * (idx_s + (NZ -2 * l  + 1)) + 4];
					Fbf[(j - ny1) * slice_bf + i * strip_bf + (n++)] = Fs[6 * (idx_s + (NZ -2 * l  + 1)) + 5];
				}

				for (l=(FDORDER/2-1);l>=1;l--) {
					idx_s =  j * slice + i * strip + l;
                    Fbf[(j - ny1) * slice_bf + i * strip_bf + (n++)] = Fs[6 * (idx_s + (NZ -2 * l  + 1)) + 2];
				}

			}
		}






	/* persistent communication see comm_ini.c*/
	/*for (i=2;i<=3;i++){*/
		/* send and reveive values at edges of the local grid */
		/*MPI_Start(&req_send[i]);
		MPI_Wait(&req_send[i],&status);
		MPI_Start(&req_rec[i]);
		MPI_Wait(&req_rec[i],&status);
	}*/			

	
	MPI_Bsend(&Ffb[NX * nf2 + nf2 + 1],NX*ny*nf2,MPI_FLOAT,INDEX[5],TAG3,MPI_COMM_WORLD);
	MPI_Bsend(&Fbf[NX * nf1 + nf1 + 1],NX*ny*nf1,MPI_FLOAT,INDEX[6],TAG4,MPI_COMM_WORLD);

	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Recv(&Fbt[NZ * nf1 + nf1 + 1], NX*NZ*nf1,MPI_FLOAT,INDEX[3],TAG6,MPI_COMM_WORLD,&status);
	

	if (POS[2]!=0) {	// no boundary exchange at top of global grid
		int n = 1;
        for (l=FDORDER/2-1;l>=1;l--) {
		    for (i=1;i<=NX;i++) {
			    // storage of top of local volume into buffer
			    for (k=1;k<=NZ;k++) {
			    	idx_s =  l * slice + i * strip + k;
			    	Fs[6 * (idx_s + (1 -2 * l) * slice) + 1] = Fbt[i * slice_bt + k * strip_bt + n + 0];
			    	Fs[6 * (idx_s + (1 -2 * l) * slice) + 3] = Fbt[i * slice_bt + k * strip_bt + n + 1];
			    	Fs[6 * (idx_s + (1 -2 * l) * slice) + 4] = Fbt[i * slice_bt + k * strip_bt + n + 2];
			    }
		    }
		    n = n + 3;
		}
		for (l = FDORDER/2; l>=(FDORDER/2);l--) {
		    for (i=1;i<=NX;i++) {
			    // storage of top of local volume into buffer
			    for (k=1;k<=NZ;k++) {
			    	idx_s =  l * slice + i * strip + k;
			    	Fs[6 * (idx_s + (1 -2 * l) * slice) + 3] = Fbt[i * slice_bt + k * strip_bt + n + 0];
			    	Fs[6 * (idx_s + (1 -2 * l) * slice) + 4] = Fbt[i * slice_bt + k * strip_bt + n + 1];
			    }
		    }
		    n = n + 2;
		}

	}

	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Recv(&Frl[NZ * nf1 + nf1 + 1], ny*NZ*nf1,MPI_FLOAT,INDEX[1],TAG2,MPI_COMM_WORLD,&status);
	MPI_Recv(&Flr[NZ * nf2 + nf2 + 1], ny*NZ*nf2,MPI_FLOAT,INDEX[2],TAG1,MPI_COMM_WORLD,&status);
	




	if ((BOUNDARY) || (POS[1]!=NPROCX-1)) {	// no boundary exchange at left edge of global grid
		int n = 1;
	    for (j=ny1+1;j<=NY;j++){
	    	n = 1;
		    for (l=1;l<=(FDORDER/2-1);l++){
			    // storage of left edge of local volume into buffer
			    for (k=1;k<=NZ;k++) {
			    	idx_s =  j * slice + l * strip + k;
			    	Fs[6 * (idx_s + NX * strip) + 0] = Flr[(j - ny1) * slice_lr + k * strip_lr + n + 0];
			    	Fs[6 * (idx_s + NX * strip) + 3] = Flr[(j - ny1) * slice_lr + k * strip_lr + n + 1];
			    	Fs[6 * (idx_s + NX * strip) + 5] = Flr[(j - ny1) * slice_lr + k * strip_lr + n + 2];
			    }

			    n += 3;
		    }
	    }
	    const_n = n;
	    for (j=ny1+1;j<=NY;j++){
	    	n = const_n;
		    for (l=(FDORDER/2);l<=FDORDER/2;l++){
			    // storage of left edge of local volume into buffer
			    for (k=1;k<=NZ;k++) {
			    	idx_s =  j * slice + l * strip + k;
			    	Fs[6 * (idx_s + NX * strip) + 0] = Flr[(j - ny1) * slice_lr + k * strip_lr + n + 0];
			    }

			    n += 2;
		    }
	    }

	}

	
	if ((BOUNDARY) || (POS[1]!=0)) {	// no boundary exchange at left edge of global grid
		int n = 1;
	    for (j=ny1+1;j<=NY;j++){
	    	n = 1;
		    for (l=(FDORDER/2-1);l>=1;l--) {
			    // storage of left edge of local volume into buffer
			    for (k=1;k<=NZ;k++) {
			    	idx_s =  j * slice + l * strip + k;
			    	Fs[6 * (idx_s + (1 -2 * l) * strip) + 0] = Frl[(j - ny1) * slice_rl + k * strip_rl + n + 0];
			    	Fs[6 * (idx_s + (1 -2 * l) * strip) + 3] = Frl[(j - ny1) * slice_rl + k * strip_rl + n + 1];
			    	Fs[6 * (idx_s + (1 -2 * l) * strip) + 5] = Frl[(j - ny1) * slice_rl + k * strip_rl + n + 2];
			    }

			    n += 3;
		    }
	    }
	    const_n = n;
	    for (j=ny1+1;j<=NY;j++){
	    	n = const_n;
		    for (l=FDORDER/2;l>=(FDORDER/2);l--){
			    // storage of left edge of local volume into buffer
			    for (k=1;k<=NZ;k++) {
			    	idx_s =  j * slice + l * strip + k;
			    	Fs[6 * (idx_s + (1 -2 * l) * strip) + 3] = Frl[(j - ny1) * slice_rl + k * strip_rl + n + 0];
			    	Fs[6 * (idx_s + (1 -2 * l) * strip) + 5] = Frl[(j - ny1) * slice_rl + k * strip_rl + n + 1];
			    }

			    n += 1;
		    }
	    }
	}
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Recv(&Ffb[NX * nf2 + nf2 + 1], NX*ny*nf2,MPI_FLOAT,INDEX[6],TAG3,MPI_COMM_WORLD,&status);
	MPI_Recv(&Fbf[NX * nf1 + nf1 + 1], NX*ny*nf1,MPI_FLOAT,INDEX[5],TAG4,MPI_COMM_WORLD,&status);
	





	if ((BOUNDARY) || (POS[3]!=NPROCZ-1))  //* no boundary exchange at back side of global grid
		for (j=ny1+1;j<=NY;j++) {
			for (i=1;i<=NX;i++) {
				n=1;
				for (l=1;l<=FDORDER/2;l++) {
					idx_s =  j * slice + i * strip + l;
					Fs[6 * (idx_s + NZ) + 2] = Ffb[(j - ny1) * slice_fb + i * strip_fb + (n++)];

				}

				for (l=1;l<=(FDORDER/2-1);l++) {
					idx_s =  j * slice + i * strip + l;
					Fs[6 * (idx_s + NZ) + 4] = Ffb[(j - ny1) * slice_fb + i * strip_fb + (n++)];
					Fs[6 * (idx_s + NZ) + 5] = Ffb[(j - ny1) * slice_fb + i * strip_fb + (n++)];
				}
			}
		}



	/* no exchange if periodic boundary condition is applied */
	if ((BOUNDARY) || (POS[3]!=0))	/* no boundary exchange at front side of global grid */
		for (j=ny1+1;j<=NY;j++) {
			for (i=1;i<=NX;i++) {
				n=1;
				for (l=FDORDER/2;l>=1;l--) {
					idx_s =  j * slice + i * strip + l;
					Fs[6 * (idx_s + (1 -2 * l)) + 4] = Fbf[(j - ny1) * slice_bf + i * strip_bf + (n++)];
					Fs[6 * (idx_s + (1 -2 * l)) + 5] = Fbf[(j - ny1) * slice_bf + i * strip_bf + (n++)];
				}

				for (l=(FDORDER/2-1);l>=1;l--) {
					idx_s =  j * slice + i * strip + l;
					Fs[6 * (idx_s + (1 -2 * l)) + 2] = Fbf[(j - ny1) * slice_bf + i * strip_bf + (n++)];
				}
			}
		}



	free_trans(Fbt, 1, NX, 1, NZ, 1, nf1);


	/*if (LOG)
	if (MYID==0){
		time2=MPI_Wtime();
		time=time2-time1;
		fprintf(FP," Real time for stress tensor exchange: \t\t %4.2f s.\n",time);
	}*/
	return time;
	
}