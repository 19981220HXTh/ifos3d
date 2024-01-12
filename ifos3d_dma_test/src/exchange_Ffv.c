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
 * exchange of wavefield velocities (frequency donain) at grid boundaries between processors
 * when using the standard staggered grid
 *  ----------------------------------------------------------------------*/

#include "fd.h"

double exchange_Ffv(float *F_Fv, int nf, float *** bufferlef_to_rig, float *** bufferrig_to_lef, float *** buffertop_to_bot, float *** bufferbot_to_top, float *** bufferfro_to_bac, float *** bufferbac_to_fro, MPI_Request * req_send, MPI_Request * req_rec, int ntr_hess) {


	extern int NX, NY, NZ, POS[4], NPROCX, NPROCY, NPROCZ, BOUNDARY, FDORDER, INDEX[7], NFMAX, ABS_TYPE;
	extern const int TAG1,TAG2,TAG3,TAG4,TAG5,TAG6;

	MPI_Status status;
	int i, j, k, l, n, nf1, nf2, m, idx, idx_f, const_n;
	double time=0.0; /*, time1=0.0;*/

    nf1=3*FDORDER/2-1;
	nf2=nf1-1;

	int ll = 1;
	if(ABS_TYPE==1 && FDORDER==2){ll=2;}
	int nval = NFMAX, nrow = NY + 2 * (ll*FDORDER/2) + 1, ncol = NX + 2 * (ll*FDORDER/2), ndep = NZ + 2 * (ll*FDORDER/2);

	float *Ftb = transform3(buffertop_to_bot, 1, NX, 1, NZ,  1, 2 * nf1);
	float *Fbt = transform3(bufferbot_to_top, 1, NX, 1, NZ,  1, 2 * nf2);
	float *Flr = transform3(bufferlef_to_rig, 1, NY, 1, NZ,  1, 2 * nf1);
	float *Frl = transform3(bufferrig_to_lef, 1, NY, 1, NZ,  1, 2 * nf2);
	float *Ffb = transform3(bufferfro_to_bac, 1, NY, 1, NX,  1, 2 * nf1);
	float *Fbf = transform3(bufferbac_to_fro, 1, NY, 1, NX,  1, 2 * nf2);


    int strip = ndep;
    int slice = ncol * ndep;
    int cube = nrow * ncol * ndep;
    int strip_tb = 2 * nf1;
    int slice_tb = NZ * 2 * nf1;
    int strip_bt = 2 * nf2;
    int slice_bt = NZ * 2 * nf2;
    int strip_lr = 2 * nf1;
    int slice_lr = NZ * 2 * nf1;
    int strip_rl = 2 * nf2;
    int slice_rl = NZ * 2 * nf2;
    int strip_fb = 2 * nf1;
    int slice_fb = NX * 2 * nf1;
    int strip_bf = 2 * nf2;
    int slice_bf = NX * 2 * nf2;



	/*if (LOG){
	if (MYID==0) time1=MPI_Wtime();}*/

	for(m=1;m<=nf*(ntr_hess+1);m++){
		/* top-bottom -----------------------------------------------------------*/


		if (POS[2]!=0) {	// no boundary exchange at top of global grid
			int n = 1;
	        for (l=1;l<=FDORDER/2-1;l++) {
			    for (i=1;i<=NX;i++) {
				    // storage of top of local volume into buffer
				    for (k=1;k<=NZ;k++) {
				    	idx_f = m * cube + l * slice + i * strip + k;
				    	Ftb[i * slice_tb + k * strip_tb + n + 0] = F_Fv[idx_f * 6 + 0];
				    	Ftb[i * slice_tb + k * strip_tb + n + 1] = F_Fv[idx_f * 6 + 1];
				    	Ftb[i * slice_tb + k * strip_tb + n + 2] = F_Fv[idx_f * 6 + 2];
				    	Ftb[i * slice_tb + k * strip_tb + n + 3] = F_Fv[idx_f * 6 + 3];
				    	Ftb[i * slice_tb + k * strip_tb + n + 4] = F_Fv[idx_f * 6 + 4];
				    	Ftb[i * slice_tb + k * strip_tb + n + 5] = F_Fv[idx_f * 6 + 5];
				    }
			    }
			    n = n + 6;
			}
			for (l = FDORDER/2; l<=(FDORDER/2);l++) {
			    for (i=1;i<=NX;i++) {
				    // storage of top of local volume into buffer
				    for (k=1;k<=NZ;k++) {
				    	idx_f = m * cube + l * slice + i * strip + k;
				    	Ftb[i * slice_tb + k * strip_tb + n + 0] = F_Fv[idx_f * 6 + 0];
				    	Ftb[i * slice_tb + k * strip_tb + n + 1] = F_Fv[idx_f * 6 + 2];
				    	Ftb[i * slice_tb + k * strip_tb + n + 2] = F_Fv[idx_f * 6 + 3];
				    	Ftb[i * slice_tb + k * strip_tb + n + 3] = F_Fv[idx_f * 6 + 5];
				    }
			    }
			    n = n + 4;
			}
		}


	    if (POS[2] != NPROCY - 1) {	// no boundary exchange at top of global grid
			int n = 1;
			for (l=FDORDER/2-1;l>=1;l--) {
			    for (i=1;i<=NX;i++) {
					// storage of top of local volume into buffer
					for (k=1;k<=NZ;k++) {
				    	idx_f = m * cube + l * slice + i * strip + k;
				    	Fbt[i * slice_bt + k * strip_bt + n + 0] = F_Fv[6 * (idx_f + (NY -2 * l  + 1) * slice) + 0];
				    	Fbt[i * slice_bt + k * strip_bt + n + 1] = F_Fv[6 * (idx_f + (NY -2 * l  + 1) * slice) + 1];
				    	Fbt[i * slice_bt + k * strip_bt + n + 2] = F_Fv[6 * (idx_f + (NY -2 * l  + 1) * slice) + 2];
				    	Fbt[i * slice_bt + k * strip_bt + n + 3] = F_Fv[6 * (idx_f + (NY -2 * l  + 1) * slice) + 3];
				    	Fbt[i * slice_bt + k * strip_bt + n + 4] = F_Fv[6 * (idx_f + (NY -2 * l  + 1) * slice) + 4];
				    	Fbt[i * slice_bt + k * strip_bt + n + 5] = F_Fv[6 * (idx_f + (NY -2 * l  + 1) * slice) + 5];
					}
				}
			    n = n + 6;
			}
			for (l = FDORDER/2; l>=(FDORDER/2);l--) {
				for (i=1;i<=NX;i++) {
				    // storage of top of local volume into buffer
					for (k=1;k<=NZ;k++) {
						idx_f = m * cube + l * slice + i * strip + k;
						Fbt[i * slice_bt + k * strip_bt + n + 0] = F_Fv[6 * (idx_f + (NY -2 * l  + 1) * slice) + 1];
						Fbt[i * slice_bt + k * strip_bt + n + 1] = F_Fv[6 * (idx_f + (NY -2 * l  + 1) * slice) + 4];
					}
				}
				n = n + 2;
			}
		}



		MPI_Sendrecv_replace(&Ftb[NZ * 2 * nf1 + 2 * nf1 + 1],2 * nf1 * NX * NZ, MPI_FLOAT, INDEX[3], TAG5, INDEX[4], TAG5, MPI_COMM_WORLD, &status);
		MPI_Sendrecv_replace(&Fbt[NZ * 2 * nf2 + 2 * nf2 + 1],2 * nf2 * NX * NZ, MPI_FLOAT, INDEX[4], TAG6, INDEX[3], TAG6, MPI_COMM_WORLD, &status);


		if (POS[2]!=NPROCY-1) {	// no boundary exchange at top of global grid
			int n = 1;
	        for (l=1;l<=FDORDER/2-1;l++) {
			    for (i=1;i<=NX;i++) {
				    // storage of top of local volume into buffer
				    for (k=1;k<=NZ;k++) {
				    	idx_f = m * cube + l * slice + i * strip + k;
				    	F_Fv[6 * (idx_f + NY * slice) + 0] = Ftb[i * slice_tb + k * strip_tb + n + 0];
				    	F_Fv[6 * (idx_f + NY * slice) + 1] = Ftb[i * slice_tb + k * strip_tb + n + 1];
				    	F_Fv[6 * (idx_f + NY * slice) + 2] = Ftb[i * slice_tb + k * strip_tb + n + 2];
				    	F_Fv[6 * (idx_f + NY * slice) + 3] = Ftb[i * slice_tb + k * strip_tb + n + 3];
				    	F_Fv[6 * (idx_f + NY * slice) + 4] = Ftb[i * slice_tb + k * strip_tb + n + 4];
				    	F_Fv[6 * (idx_f + NY * slice) + 5] = Ftb[i * slice_tb + k * strip_tb + n + 5];
				    }
			    }
			    n = n + 6;
			}
			for (l = FDORDER/2; l<=(FDORDER/2);l++) {
			    for (i=1;i<=NX;i++) {
				    // storage of top of local volume into buffer
				    for (k=1;k<=NZ;k++) {
				    	idx_f = m * cube + l * slice + i * strip + k;
				    	F_Fv[6 * (idx_f + NY * slice) + 0] = Ftb[i * slice_tb + k * strip_tb + n + 0];
				    	F_Fv[6 * (idx_f + NY * slice) + 2] = Ftb[i * slice_tb + k * strip_tb + n + 1];
				    	F_Fv[6 * (idx_f + NY * slice) + 3] = Ftb[i * slice_tb + k * strip_tb + n + 2];
				    	F_Fv[6 * (idx_f + NY * slice) + 5] = Ftb[i * slice_tb + k * strip_tb + n + 3];
				    	//if(MYID == 0)
				    	//printf("--##-------------------->funcname:%s, lineNum:%d, l:%d, i:%d, k, n:%d:%d\n", __FUNCTION__, __LINE__, l, i, k,n);
				    }
			    }
			    n = n + 4;
			}

		}


		if (POS[2]!=0) {	// no boundary exchange at top of global grid
			int n = 1;
	        for (l=FDORDER/2-1;l>=1;l--) {
			    for (i=1;i<=NX;i++) {
				    // storage of top of local volume into buffer
				    for (k=1;k<=NZ;k++) {
				    	idx_f = m * cube + l * slice + i * strip + k;
				    	F_Fv[6 * (idx_f + (1 -2 * l) * slice) + 0] = Fbt[i * slice_bt + k * strip_bt + n + 0];
				    	F_Fv[6 * (idx_f + (1 -2 * l) * slice) + 1] = Fbt[i * slice_bt + k * strip_bt + n + 1];
				    	F_Fv[6 * (idx_f + (1 -2 * l) * slice) + 2] = Fbt[i * slice_bt + k * strip_bt + n + 2];
				    	F_Fv[6 * (idx_f + (1 -2 * l) * slice) + 3] = Fbt[i * slice_bt + k * strip_bt + n + 3];
				    	F_Fv[6 * (idx_f + (1 -2 * l) * slice) + 4] = Fbt[i * slice_bt + k * strip_bt + n + 4];
				    	F_Fv[6 * (idx_f + (1 -2 * l) * slice) + 5] = Fbt[i * slice_bt + k * strip_bt + n + 5];
				    }
			    }
			    n = n + 6;
			}
			for (l = FDORDER/2; l>=(FDORDER/2);l--) {
			    for (i=1;i<=NX;i++) {
				    // storage of top of local volume into buffer
				    for (k=1;k<=NZ;k++) {
				    	idx_f = m * cube + l * slice + i * strip + k;
				    	F_Fv[6 * (idx_f + (1 -2 * l) * slice) + 1] = Fbt[i * slice_bt + k * strip_bt + n + 0];
				    	F_Fv[6 * (idx_f + (1 -2 * l) * slice) + 4] = Fbt[i * slice_bt + k * strip_bt + n + 1];
				    }
			    }
			    n = n + 2;
			}

		}







        /***************************************** left to right     ******************************************/
	    if ((BOUNDARY) || (POS[1]!=0)) {	// no boundary exchange at left edge of global grid
			int n = 1;
		    for (j=1;j<=NY;j++){
		    	n = 1;
			    for (l=1;l<=(FDORDER/2-1);l++){
				    // storage of left edge of local volume into buffer
				    for (k=1;k<=NZ;k++) {
				    	idx_f = m * cube + j * slice + l * strip + k;
				    	Flr[j * slice_lr + k * strip_lr + n + 0] = F_Fv[6 * idx_f + 0];
				    	Flr[j * slice_lr + k * strip_lr + n + 1] = F_Fv[6 * idx_f + 1];
				    	Flr[j * slice_lr + k * strip_lr + n + 2] = F_Fv[6 * idx_f + 2];
				    	Flr[j * slice_lr + k * strip_lr + n + 3] = F_Fv[6 * idx_f + 3];
				    	Flr[j * slice_lr + k * strip_lr + n + 4] = F_Fv[6 * idx_f + 4];
				    	Flr[j * slice_lr + k * strip_lr + n + 5] = F_Fv[6 * idx_f + 5];
				    }

				    n += 6;
			    }
		    }
		    const_n = n;
		    for (j=1;j<=NY;j++){
		    	n = const_n;
			    for (l=(FDORDER/2);l<=FDORDER/2;l++){
				    // storage of left edge of local volume into buffer
				    for (k=1;k<=NZ;k++) {
				    	idx_f = m * cube + j * slice + l * strip + k;
				    	Flr[j * slice_lr + k * strip_lr + n + 0] = F_Fv[6 * idx_f + 1];
				    	Flr[j * slice_lr + k * strip_lr + n + 1] = F_Fv[6 * idx_f + 2];
				    	Flr[j * slice_lr + k * strip_lr + n + 2] = F_Fv[6 * idx_f + 4];
				    	Flr[j * slice_lr + k * strip_lr + n + 3] = F_Fv[6 * idx_f + 5];
				    }

				    n += 4;
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
				    	idx_f = m * cube + j * slice + l * strip + k;
				    	Frl[j * slice_rl + k * strip_rl + n + 0] = F_Fv[6 * (idx_f + (NX -2 * l  + 1) * strip) + 0];
				    	Frl[j * slice_rl + k * strip_rl + n + 1] = F_Fv[6 * (idx_f + (NX -2 * l  + 1) * strip) + 1];
				    	Frl[j * slice_rl + k * strip_rl + n + 2] = F_Fv[6 * (idx_f + (NX -2 * l  + 1) * strip) + 2];
				    	Frl[j * slice_rl + k * strip_rl + n + 3] = F_Fv[6 * (idx_f + (NX -2 * l  + 1) * strip) + 3];
				    	Frl[j * slice_rl + k * strip_rl + n + 4] = F_Fv[6 * (idx_f + (NX -2 * l  + 1) * strip) + 4];
				    	Frl[j * slice_rl + k * strip_rl + n + 5] = F_Fv[6 * (idx_f + (NX -2 * l  + 1) * strip) + 5];
				    }

				    n += 6;
			    }
		    }
		    const_n = n;
		    for (j=1;j<=NY;j++){
		    	n = const_n;
			    for (l=FDORDER/2;l>=(FDORDER/2);l--){
				    // storage of left edge of local volume into buffer
				    for (k=1;k<=NZ;k++) {
				    	idx_f = m * cube + j * slice + l * strip + k;
				    	Frl[j * slice_rl + k * strip_rl + n + 0] = F_Fv[6 * (idx_f + (NX -2 * l  + 1) * strip) + 0];
				    	Frl[j * slice_rl + k * strip_rl + n + 1] = F_Fv[6 * (idx_f + (NX -2 * l  + 1) * strip) + 3];
				    }

				    n += 2;
			    }
		    }
		}



		MPI_Sendrecv_replace(&Flr[NZ * 2 * nf1 + 2 * nf1 + 1],NY*NZ*2*nf1,MPI_FLOAT,INDEX[1],TAG1,INDEX[2],TAG1,MPI_COMM_WORLD,&status);
		MPI_Sendrecv_replace(&Frl[NZ * 2 * nf2 + 2 * nf2 + 1],NY*NZ*2*nf2,MPI_FLOAT,INDEX[2],TAG2,INDEX[1],TAG2,MPI_COMM_WORLD,&status);


		if ((BOUNDARY) || (POS[1]!=NPROCX-1)) {	// no boundary exchange at left edge of global grid
			int n = 1;
		    for (j=1;j<=NY;j++){
		    	n = 1;
			    for (l=1;l<=(FDORDER/2-1);l++){
				    // storage of left edge of local volume into buffer
				    for (k=1;k<=NZ;k++) {
				    	idx_f = m * cube + j * slice + l * strip + k;
				    	F_Fv[6 * (idx_f + NX * strip) + 0] = Flr[j * slice_lr + k * strip_lr + n + 0];
				    	F_Fv[6 * (idx_f + NX * strip) + 1] = Flr[j * slice_lr + k * strip_lr + n + 1];
				    	F_Fv[6 * (idx_f + NX * strip) + 2] = Flr[j * slice_lr + k * strip_lr + n + 2];
				    	F_Fv[6 * (idx_f + NX * strip) + 3] = Flr[j * slice_lr + k * strip_lr + n + 3];
				    	F_Fv[6 * (idx_f + NX * strip) + 4] = Flr[j * slice_lr + k * strip_lr + n + 4];
				    	F_Fv[6 * (idx_f + NX * strip) + 5] = Flr[j * slice_lr + k * strip_lr + n + 5];
				    }

				    n += 6;
			    }
		    }
		    const_n = n;
		    for (j=1;j<=NY;j++){
		    	n = const_n;
			    for (l=(FDORDER/2);l<=FDORDER/2;l++){
				    // storage of left edge of local volume into buffer
				    for (k=1;k<=NZ;k++) {
				    	idx_f = m * cube + j * slice + l * strip + k;
				    	F_Fv[6 * (idx_f + NX * strip) + 1] = Flr[j * slice_lr + k * strip_lr + n + 0];
				    	F_Fv[6 * (idx_f + NX * strip) + 2] = Flr[j * slice_lr + k * strip_lr + n + 1];
				    	F_Fv[6 * (idx_f + NX * strip) + 4] = Flr[j * slice_lr + k * strip_lr + n + 2];
				    	F_Fv[6 * (idx_f + NX * strip) + 5] = Flr[j * slice_lr + k * strip_lr + n + 3];
				    }

				    n += 4;
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
				    	idx_f = m * cube + j * slice + l * strip + k;
				    	F_Fv[6 * (idx_f + (1 -2 * l) * strip) + 0] = Frl[j * slice_rl + k * strip_rl + n + 0];
				    	F_Fv[6 * (idx_f + (1 -2 * l) * strip) + 1] = Frl[j * slice_rl + k * strip_rl + n + 1];
				    	F_Fv[6 * (idx_f + (1 -2 * l) * strip) + 2] = Frl[j * slice_rl + k * strip_rl + n + 2];
				    	F_Fv[6 * (idx_f + (1 -2 * l) * strip) + 3] = Frl[j * slice_rl + k * strip_rl + n + 3];
				    	F_Fv[6 * (idx_f + (1 -2 * l) * strip) + 4] = Frl[j * slice_rl + k * strip_rl + n + 4];
				    	F_Fv[6 * (idx_f + (1 -2 * l) * strip) + 5] = Frl[j * slice_rl + k * strip_rl + n + 5];
				    }

				    n += 6;
			    }
		    }
		    const_n = n;
		    for (j=1;j<=NY;j++){
		    	n = const_n;
			    for (l=FDORDER/2;l>=(FDORDER/2);l--){
				    // storage of left edge of local volume into buffer
				    for (k=1;k<=NZ;k++) {
				    	idx_f = m * cube + j * slice + l * strip + k;
				    	F_Fv[6 * (idx_f + (1 -2 * l) * strip) + 0] = Frl[j * slice_rl + k * strip_rl + n + 0];
				    	F_Fv[6 * (idx_f + (1 -2 * l) * strip) + 3] = Frl[j * slice_rl + k * strip_rl + n + 1];
				    }

				    n += 2;
			    }
		    }
		}





		/* front-back -----------------------------------------------------------*/



		if ((BOUNDARY) || (POS[3]!=0))	//* no boundary exchange at front side of global grid
			for (j=1;j<=NY;j++) {
				for (i=1;i<=NX;i++) {   //* storage of front side of local volume into buffer
					n=1;
				    for (l=1;l<=FDORDER/2;l++){
				    	idx_f = m * cube + j * slice + i * strip + l;
				    	Ffb[j * slice_fb + i * strip_fb + (n++)] = F_Fv[6 * idx_f + 0];
				    	Ffb[j * slice_fb + i * strip_fb + (n++)] = F_Fv[6 * idx_f + 1];
				    	Ffb[j * slice_fb + i * strip_fb + (n++)] = F_Fv[6 * idx_f + 3];
				    	Ffb[j * slice_fb + i * strip_fb + (n++)] = F_Fv[6 * idx_f + 4];
				    }
				    for (l=1;l<=(FDORDER/2-1);l++) {
				    	idx_f = m * cube + j * slice + i * strip + l;
				    	Ffb[j * slice_fb + i * strip_fb + (n++)] = F_Fv[6 * idx_f + 2];
				    	Ffb[j * slice_fb + i * strip_fb + (n++)] = F_Fv[6 * idx_f + 5];
				    }
				}
			}


		// no exchange if periodic boundary condition is applied
		if ((BOUNDARY) || (POS[3]!=NPROCZ-1))	//* no boundary exchange at back side of global grid
			for (j=1;j<=NY;j++) {
				for (i=1;i<=NX;i++) {
						//* storage of back side of local volume into buffer
					n=1;
					for (l=FDORDER/2;l>=1;l--) {
						idx_f = m * cube + j * slice + i * strip + l;
						Fbf[j * slice_bf + i * strip_bf + (n++)] = F_Fv[6 * (idx_f + (NZ -2 * l  + 1)) + 2];
						Fbf[j * slice_bf + i * strip_bf + (n++)] = F_Fv[6 * (idx_f + (NZ -2 * l  + 1)) + 5];
					}

					for (l=(FDORDER/2-1);l>=1;l--) {
						idx_f = m * cube + j * slice + i * strip + l;
	                    Fbf[j * slice_bf + i * strip_bf + (n++)] = F_Fv[6 * (idx_f + (NZ -2 * l  + 1)) + 0];
	                    Fbf[j * slice_bf + i * strip_bf + (n++)] = F_Fv[6 * (idx_f + (NZ -2 * l  + 1)) + 1];
	                    Fbf[j * slice_bf + i * strip_bf + (n++)] = F_Fv[6 * (idx_f + (NZ -2 * l  + 1)) + 3];
	                    Fbf[j * slice_bf + i * strip_bf + (n++)] = F_Fv[6 * (idx_f + (NZ -2 * l  + 1)) + 4];
					}

				}
			}


		MPI_Sendrecv_replace(&Ffb[NX * 2 * nf1 + 2 * nf1 + 1],NX*NY*2*nf1,MPI_FLOAT,INDEX[5],TAG3,INDEX[6],TAG3,MPI_COMM_WORLD,&status);
		MPI_Sendrecv_replace(&Fbf[NX * 2 * nf2 + 2 * nf2 + 1],NX*NY*2*nf2,MPI_FLOAT,INDEX[6],TAG4,INDEX[5],TAG4,MPI_COMM_WORLD,&status);


		// no exchange if periodic boundary condition is applied
		if ((BOUNDARY) || (POS[3]!=NPROCZ-1))  // no boundary exchange at back side of global grid
			for (j=1;j<=NY;j++) {
				for (i=1;i<=NX;i++) {
					n=1;
					for (l=1;l<=FDORDER/2;l++) {
						idx_f = m * cube + j * slice + i * strip + l;
						F_Fv[6 * (idx_f + NZ) + 0] = Ffb[j * slice_fb + i * strip_fb + (n++)];
						F_Fv[6 * (idx_f + NZ) + 1] = Ffb[j * slice_fb + i * strip_fb + (n++)];
						F_Fv[6 * (idx_f + NZ) + 3] = Ffb[j * slice_fb + i * strip_fb + (n++)];
						F_Fv[6 * (idx_f + NZ) + 4] = Ffb[j * slice_fb + i * strip_fb + (n++)];
					}

					for (l=1;l<=(FDORDER/2-1);l++) {
						idx_f = m * cube + j * slice + i * strip + l;
						F_Fv[6 * (idx_f + NZ) + 2] = Ffb[j * slice_fb + i * strip_fb + (n++)];
						F_Fv[6 * (idx_f + NZ) + 5] = Ffb[j * slice_fb + i * strip_fb + (n++)];
					}
				}
			}



		// no exchange if periodic boundary condition is applied
		if ((BOUNDARY) || (POS[3]!=0))	// no boundary exchange at front side of global grid
			for (j=1;j<=NY;j++) {
				for (i=1;i<=NX;i++) {
					n=1;
					for (l=FDORDER/2;l>=1;l--) {
						idx_f = m * cube + j * slice + i * strip + l;
						F_Fv[6 * (idx_f + (1 -2 * l)) + 2] = Fbf[j * slice_bf + i * strip_bf + (n++)];
						F_Fv[6 * (idx_f + (1 -2 * l)) + 5] = Fbf[j * slice_bf + i * strip_bf + (n++)];
					}

					for (l=(FDORDER/2-1);l>=1;l--) {
						idx_f = m * cube + j * slice + i * strip + l;
						F_Fv[6 * (idx_f + (1 -2 * l)) + 0] = Fbf[j * slice_bf + i * strip_bf + (n++)];
						F_Fv[6 * (idx_f + (1 -2 * l)) + 1] = Fbf[j * slice_bf + i * strip_bf + (n++)];
						F_Fv[6 * (idx_f + (1 -2 * l)) + 3] = Fbf[j * slice_bf + i * strip_bf + (n++)];
						F_Fv[6 * (idx_f + (1 -2 * l)) + 4] = Fbf[j * slice_bf + i * strip_bf + (n++)];
					}
				}
			}


	}




	free_trans(Ftb, 1, NX, 1, NZ, 1, 2 * nf1);
	free_trans(Fbt, 1, NX, 1, NZ, 1, 2 * nf2);
	free_trans(Flr, 1, NY, 1, NZ, 1, 2 * nf1);
	free_trans(Frl, 1, NY, 1, NZ, 1, 2 * nf2);
	free_trans(Ffb, 1, NY, 1, NX, 1, 2 * nf1);
	free_trans(Fbf, 1, NY, 1, NX, 1, 2 * nf2);


	/*if (LOG)
	if (MYID==0){
		time2=MPI_Wtime();
		time=time2-time1;
		fprintf(FP," Real time for particle velocity exchange: \t %4.2f s.\n",time);
	}*/
	return time;

}
