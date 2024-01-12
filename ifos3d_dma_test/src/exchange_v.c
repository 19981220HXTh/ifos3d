/* $Id: exchange_v.c,v 1.1.1.1 2009/12/09 14:29:22 sjetschny Exp $ */
/*------------------------------------------------------------------------
 * exchange of particle velocities at grid boundaries between processors
 * when using the standard staggered grid
 *
 *  ----------------------------------------------------------------------*/

#include "fd.h"

double exchange_v(float *Fv, float *** bufferlef_to_rig, float *** bufferrig_to_lef, float *** buffertop_to_bot, float *** bufferbot_to_top, float *** bufferfro_to_bac, float *** bufferbac_to_fro, MPI_Request * req_send, MPI_Request * req_rec) {

	extern int NX, NY, NZ, POS[4], NPROCX, NPROCY, NPROCZ, BOUNDARY,  FDORDER,  INDEX[7], ABS_TYPE;  extern int MYID; /*MYID,LOG,*/
	extern const int TAG1,TAG2,TAG3,TAG4,TAG5,TAG6;

	MPI_Status status;
	int i, j, k, l, n, nf1, nf2, const_n;
	double time=0.0;  /*, time1=0.0;*/

    nf1=3*FDORDER/2-1;
	nf2=nf1-1;
    //printf("--##-------------------->funcname:%s, lineNum:%d\n", __FUNCTION__, __LINE__);
	int idx_v;
	int ll = 1;
	if(ABS_TYPE==1 && FDORDER==2){ll=2;}
	int nrow = NY + 2 * (ll*FDORDER/2) + 1, ncol = NX + 2 * (ll*FDORDER/2), ndep = NZ + 2 * (ll*FDORDER/2);
	int size = nrow * ncol * ndep;
	float *Ftb = transform3(buffertop_to_bot, 1, NX, 1, NZ,  1, nf1);
	float *Fbt = transform3(bufferbot_to_top, 1, NX, 1, NZ,  1, nf2);
	float *Flr = transform3(bufferlef_to_rig, 1, NY, 1, NZ,  1, nf1);
	float *Frl = transform3(bufferrig_to_lef, 1, NY, 1, NZ,  1, nf2);
	float *Ffb = transform3(bufferfro_to_bac, 1, NY, 1, NX,  1, nf1);
	float *Fbf = transform3(bufferbac_to_fro, 1, NY, 1, NX,  1, nf2);
    int strip = ndep;
    int slice = ncol * ndep;
    int strip_tb = nf1;
    int slice_tb = NZ * nf1; 
    int strip_bt = nf2;
    int slice_bt = NZ * nf2;
    int strip_lr = nf1;
    int slice_lr = NZ * nf1;
    int strip_rl = nf2;
    int slice_rl = NZ * nf2;
    int strip_fb = nf1;
    int slice_fb = NX * nf1;
    int strip_bf = nf2;
    int slice_bf = NX * nf2;

	/*if (LOG){
	if (MYID==0) time1=MPI_Wtime();}*/

	/* top-bottom -----------------------------------------------------------*/


	if (POS[2]!=0) {	// no boundary exchange at top of global grid
		int n = 1;
        for (l=1;l<=FDORDER/2-1;l++) {
		    for (i=1;i<=NX;i++) {
			    // storage of top of local volume into buffer
			    for (k=1;k<=NZ;k++) {
			    	idx_v =  l * slice + i * strip + k;
			    	Ftb[i * slice_tb + k * strip_tb + n + 0] = Fv[3 * idx_v + 0];
			    	Ftb[i * slice_tb + k * strip_tb + n + 1] = Fv[3 * idx_v + 1];
			    	Ftb[i * slice_tb + k * strip_tb + n + 2] = Fv[3 * idx_v + 2];
			    }
		    }
		    n = n + 3;
		}
		for (l = FDORDER/2; l<=(FDORDER/2);l++) {
		    for (i=1;i<=NX;i++) {
			    // storage of top of local volume into buffer
			    for (k=1;k<=NZ;k++) {
			    	idx_v =  l * slice + i * strip + k;
			    	Ftb[i * slice_tb + k * strip_tb + n + 0] = Fv[3 * idx_v + 0];
			    	Ftb[i * slice_tb + k * strip_tb + n + 1] = Fv[3 * idx_v + 2];
			    }
		    }
		    n = n + 2;
		}
	}


    if (POS[2] != NPROCY - 1) {	// no boundary exchange at top of global grid
		int n = 1;
		for (l=FDORDER/2-1;l>=1;l--) {
		    for (i=1;i<=NX;i++) {
				// storage of top of local volume into buffer
				for (k=1;k<=NZ;k++) {
			    	idx_v =  l * slice + i * strip + k;
			    	Fbt[i * slice_bt + k * strip_bt + n + 0] = Fv[3 * (idx_v + (NY -2 * l  + 1) * slice) + 0];
			    	Fbt[i * slice_bt + k * strip_bt + n + 1] = Fv[3 * (idx_v + (NY -2 * l  + 1) * slice) + 1];
			    	Fbt[i * slice_bt + k * strip_bt + n + 2] = Fv[3 * (idx_v + (NY -2 * l  + 1) * slice) + 2];
				}
			}
		    n = n + 3;
		}
		for (l = FDORDER/2; l>=(FDORDER/2);l--) {
			for (i=1;i<=NX;i++) {
			    // storage of top of local volume into buffer
				for (k=1;k<=NZ;k++) {
			    	idx_v =  l * slice + i * strip + k;
			    	Fbt[i * slice_bt + k * strip_bt + n + 0] = Fv[3 * (idx_v + (NY -2 * l  + 1) * slice) + 1];
				}
			}
			n = n + 1;
		}
	}


	//MPI_Sendrecv_replace(&buffertop_to_bot[1][1][1],nf1*NX*NZ,MPI_FLOAT,INDEX[3],TAG5,INDEX[4],TAG5,MPI_COMM_WORLD,&status);
    //	MPI_Sendrecv_replace(&bufferbot_to_top[1][1][1], NX * NZ * nf2, MPI_FLOAT, INDEX[4], TAG6, INDEX[3], TAG6, MPI_COMM_WORLD, &status);
	MPI_Sendrecv_replace(&Ftb[NZ * nf1 + nf1 + 1],nf1*NX*NZ,MPI_FLOAT,INDEX[3],TAG5,INDEX[4],TAG5,MPI_COMM_WORLD,&status);
	MPI_Sendrecv_replace(&Fbt[NZ * nf2 + nf2 + 1],nf2 * NX * NZ, MPI_FLOAT, INDEX[4], TAG6, INDEX[3], TAG6, MPI_COMM_WORLD, &status);



	if (POS[2]!=NPROCY-1) {	// no boundary exchange at top of global grid
		int n = 1;
        for (l=1;l<=FDORDER/2-1;l++) {
		    for (i=1;i<=NX;i++) {
			    // storage of top of local volume into buffer
			    for (k=1;k<=NZ;k++) {
			    	idx_v =  l * slice + i * strip + k;
			    	Fv[3 * (idx_v + NY * slice) + 0] = Ftb[i * slice_tb + k * strip_tb + n + 0];
			    	Fv[3 * (idx_v + NY * slice) + 1] = Ftb[i * slice_tb + k * strip_tb + n + 1];
			    	Fv[3 * (idx_v + NY * slice) + 2] = Ftb[i * slice_tb + k * strip_tb + n + 2];
			    }
		    }
		    n = n + 3;
		}
		for (l = FDORDER/2; l<=(FDORDER/2);l++) {
		    for (i=1;i<=NX;i++) {
			    // storage of top of local volume into buffer
			    for (k=1;k<=NZ;k++) {
			    	idx_v =  l * slice + i * strip + k;
			    	Fv[3 * (idx_v + NY * slice) + 0] = Ftb[i * slice_tb + k * strip_tb + n + 0];
			    	Fv[3 * (idx_v + NY * slice) + 2] = Ftb[i * slice_tb + k * strip_tb + n + 1];
			    	//if(MYID == 0)
			    	//printf("--##-------------------->funcname:%s, lineNum:%d, l:%d, i:%d, k, n:%d:%d\n", __FUNCTION__, __LINE__, l, i, k,n);
			    }
		    }
		    n = n + 2;
		}

	}


	if (POS[2]!=0) {	// no boundary exchange at top of global grid
		int n = 1;
        for (l=FDORDER/2-1;l>=1;l--) {
		    for (i=1;i<=NX;i++) {
			    // storage of top of local volume into buffer
			    for (k=1;k<=NZ;k++) {
			    	idx_v =  l * slice + i * strip + k;
			    	Fv[3 * (idx_v + (1 -2 * l) * slice) + 0] = Fbt[i * slice_bt + k * strip_bt + n + 0];
			    	Fv[3 * (idx_v + (1 -2 * l) * slice) + 1] = Fbt[i * slice_bt + k * strip_bt + n + 1];
			    	Fv[3 * (idx_v + (1 -2 * l) * slice) + 2] = Fbt[i * slice_bt + k * strip_bt + n + 2];
			    }
		    }
		    n = n + 3;
		}
		for (l = FDORDER/2; l>=(FDORDER/2);l--) {
		    for (i=1;i<=NX;i++) {
			    // storage of top of local volume into buffer
			    for (k=1;k<=NZ;k++) {
			    	idx_v =  l * slice + i * strip + k;
			    	Fv[3 * (idx_v + (1 -2 * l) * slice) + 1] = Fbt[i * slice_bt + k * strip_bt + n + 0];
			    }
		    }
		    n = n + 1;
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
			    	idx_v =  j * slice + l * strip + k;
			    	Flr[j * slice_lr + k * strip_lr + n + 0] = Fv[3 * idx_v + 0];
			    	Flr[j * slice_lr + k * strip_lr + n + 1] = Fv[3 * idx_v + 1];
			    	Flr[j * slice_lr + k * strip_lr + n + 2] = Fv[3 * idx_v + 2];
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
			    	idx_v =  j * slice + l * strip + k;
			    	Flr[j * slice_lr + k * strip_lr + n + 0] = Fv[3 * idx_v + 1];
			    	Flr[j * slice_lr + k * strip_lr + n + 1] = Fv[3 * idx_v + 2];
			    }

			    n += 2;
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
			    	idx_v =  j * slice + l * strip + k;
			    	Frl[j * slice_rl + k * strip_rl + n + 0] = Fv[3 * (idx_v + (NX -2 * l  + 1) * strip) + 0];
			    	Frl[j * slice_rl + k * strip_rl + n + 1] = Fv[3 * (idx_v + (NX -2 * l  + 1) * strip) + 1];
			    	Frl[j * slice_rl + k * strip_rl + n + 2] = Fv[3 * (idx_v + (NX -2 * l  + 1) * strip) + 2];
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
			    	idx_v =  j * slice + l * strip + k;
			    	Frl[j * slice_rl + k * strip_rl + n + 0] = Fv[3 * (idx_v + (NX -2 * l  + 1) * strip) + 0];
			    }

			    n += 1;
		    }
	    }
	}


	//MPI_Sendrecv_replace(&bufferlef_to_rig[1][1][1],NY*NZ*nf1,MPI_FLOAT,INDEX[1],TAG1,INDEX[2],TAG1,MPI_COMM_WORLD,&status);
	//MPI_Sendrecv_replace(&bufferrig_to_lef[1][1][1],NY*NZ*nf2,MPI_FLOAT,INDEX[2],TAG2,INDEX[1],TAG2,MPI_COMM_WORLD,&status);
	MPI_Sendrecv_replace(&Flr[NZ * nf1 + nf1 + 1],NY*NZ*nf1,MPI_FLOAT,INDEX[1],TAG1,INDEX[2],TAG1,MPI_COMM_WORLD,&status);
	MPI_Sendrecv_replace(&Frl[NZ * nf2 + nf2 + 1],NY*NZ*nf2,MPI_FLOAT,INDEX[2],TAG2,INDEX[1],TAG2,MPI_COMM_WORLD,&status);






	if ((BOUNDARY) || (POS[1]!=NPROCX-1)) {	// no boundary exchange at left edge of global grid
		int n = 1;
	    for (j=1;j<=NY;j++){
	    	n = 1;
		    for (l=1;l<=(FDORDER/2-1);l++){
			    // storage of left edge of local volume into buffer
			    for (k=1;k<=NZ;k++) {
			    	idx_v =  j * slice + l * strip + k;
			    	Fv[3 * (idx_v + NX * strip) + 0] = Flr[j * slice_lr + k * strip_lr + n + 0];
			    	Fv[3 * (idx_v + NX * strip) + 1] = Flr[j * slice_lr + k * strip_lr + n + 1];
			    	Fv[3 * (idx_v + NX * strip) + 2] = Flr[j * slice_lr + k * strip_lr + n + 2];
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
			    	idx_v =  j * slice + l * strip + k;
			    	Fv[3 * (idx_v + NX * strip) + 1] = Flr[j * slice_lr + k * strip_lr + n + 0];
			    	Fv[3 * (idx_v + NX * strip) + 2] = Flr[j * slice_lr + k * strip_lr + n + 1];
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
			    	idx_v =  j * slice + l * strip + k;
			    	Fv[3 * (idx_v + (1 -2 * l) * strip) + 0] = Frl[j * slice_rl + k * strip_rl + n + 0];
			    	Fv[3 * (idx_v + (1 -2 * l) * strip) + 1] = Frl[j * slice_rl + k * strip_rl + n + 1];
			    	Fv[3 * (idx_v + (1 -2 * l) * strip) + 2] = Frl[j * slice_rl + k * strip_rl + n + 2];
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
			    	idx_v =  j * slice + l * strip + k;
			    	Fv[3 * (idx_v + (1 -2 * l) * strip) + 0] = Frl[j * slice_rl + k * strip_rl + n + 0];
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
			    	idx_v =  j * slice + i * strip + l;
			    	Ffb[j * slice_fb + i * strip_fb + (n++)] = Fv[3 * idx_v + 0];
			    	Ffb[j * slice_fb + i * strip_fb + (n++)] = Fv[3 * idx_v + 1];
			    }
			    for (l=1;l<=(FDORDER/2-1);l++) {
			    	idx_v =  j * slice + i * strip + l;
			    	Ffb[j * slice_fb + i * strip_fb + (n++)] = Fv[3 * idx_v + 2];
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
					idx_v =  j * slice + i * strip + l;
					Fbf[j * slice_bf + i * strip_bf + (n++)] = Fv[3 * (idx_v + (NZ -2 * l  + 1)) + 2];
				}

				for (l=(FDORDER/2-1);l>=1;l--) {
					idx_v =  j * slice + i * strip + l;
                    Fbf[j * slice_bf + i * strip_bf + (n++)] = Fv[3 * (idx_v + (NZ -2 * l  + 1)) + 0];
                    Fbf[j * slice_bf + i * strip_bf + (n++)] = Fv[3 * (idx_v + (NZ -2 * l  + 1)) + 1];
				}

			}
		}



	//MPI_Sendrecv_replace(&bufferfro_to_bac[1][1][1],NX*NY*nf1,MPI_FLOAT,INDEX[5],TAG3,INDEX[6],TAG3,MPI_COMM_WORLD,&status);
	//MPI_Sendrecv_replace(&bufferbac_to_fro[1][1][1],NX*NY*nf2,MPI_FLOAT,INDEX[6],TAG4,INDEX[5],TAG4,MPI_COMM_WORLD,&status);
	MPI_Sendrecv_replace(&Ffb[NX * nf1 + nf1 + 1],NX*NY*nf1,MPI_FLOAT,INDEX[5],TAG3,INDEX[6],TAG3,MPI_COMM_WORLD,&status);
	MPI_Sendrecv_replace(&Fbf[NX * nf2 + nf2 + 1],NX*NY*nf2,MPI_FLOAT,INDEX[6],TAG4,INDEX[5],TAG4,MPI_COMM_WORLD,&status);


	/* no exchange if periodic boundary condition is applied */

	if ((BOUNDARY) || (POS[3]!=NPROCZ-1))  //* no boundary exchange at back side of global grid
		for (j=1;j<=NY;j++) {
			for (i=1;i<=NX;i++) {
				n=1;
				for (l=1;l<=FDORDER/2;l++) {
					idx_v =  j * slice + i * strip + l;
					Fv[3 * (idx_v + NZ) + 0] = Ffb[j * slice_fb + i * strip_fb + (n++)];
					Fv[3 * (idx_v + NZ) + 1] = Ffb[j * slice_fb + i * strip_fb + (n++)];
				}

				for (l=1;l<=(FDORDER/2-1);l++) {
					idx_v =  j * slice + i * strip + l;
					Fv[3 * (idx_v + NZ) + 2] = Ffb[j * slice_fb + i * strip_fb + (n++)];
				}
			}
		}



	/* no exchange if periodic boundary condition is applied */
	if ((BOUNDARY) || (POS[3]!=0))	/* no boundary exchange at front side of global grid */
		for (j=1;j<=NY;j++) {
			for (i=1;i<=NX;i++) {
				n=1;
				for (l=FDORDER/2;l>=1;l--) {
					idx_v =  j * slice + i * strip + l;
					Fv[3 * (idx_v + (1 -2 * l)) + 2] = Fbf[j * slice_bf + i * strip_bf + (n++)];
				}

				for (l=(FDORDER/2-1);l>=1;l--) {
					idx_v =  j * slice + i * strip + l;
					Fv[3 * (idx_v + (1 -2 * l)) + 0] = Fbf[j * slice_bf + i * strip_bf + (n++)];
					Fv[3 * (idx_v + (1 -2 * l)) + 1] = Fbf[j * slice_bf + i * strip_bf + (n++)];
				}
			}
		}



	free_trans(Ftb, 1, NX, 1, NZ, 1, nf1);
	free_trans(Fbt, 1, NX, 1, NZ, 1, nf2);
	free_trans(Flr, 1, NY, 1, NZ, 1, nf1);
	free_trans(Frl, 1, NY, 1, NZ, 1, nf2);
	free_trans(Ffb, 1, NY, 1, NX, 1, nf1);
	free_trans(Fbf, 1, NY, 1, NX, 1, nf2);

	return time;

}


double exchange_v2_1(float *Fv, float *** buffertop_to_bot,  MPI_Request * req_send, MPI_Request * req_rec) {

	extern int NX, NY, NZ, POS[4], NPROCX, NPROCY, NPROCZ, BOUNDARY,  FDORDER,  INDEX[7], ABS_TYPE;  extern int MYID; /*MYID,LOG,*/
	extern const int TAG1,TAG2,TAG3,TAG4,TAG5,TAG6;

	MPI_Status status;
	int i, j, k, l, n, nf1, nf2, const_n;
	double time=0.0;  /*, time1=0.0;*/

    nf1=3*FDORDER/2-1;
	nf2=nf1-1;
    //printf("--##-------------------->funcname:%s, lineNum:%d\n", __FUNCTION__, __LINE__);
	int idx_v;
	int ll = 1;
	if(ABS_TYPE==1 && FDORDER==2){ll=2;}
	int nrow = NY + 2 * (ll*FDORDER/2) + 1, ncol = NX + 2 * (ll*FDORDER/2), ndep = NZ + 2 * (ll*FDORDER/2);
	int size = nrow * ncol * ndep;
	float *Ftb = transform3(buffertop_to_bot, 1, NX, 1, NZ,  1, nf1);
	int ny = ((NY % 2 == 0) ? (NY/2):(NY/2+1));
	float Flr[(ny+1)*(NZ+1)*(nf1+1)],Frl[(ny+1)*(NZ+1)*(nf2+1)],Ffb[(ny+1)*(NX+1)*(nf1+1)],Fbf[(ny+1)*(NX+1)*(nf2+1)];
    int strip = ndep;
    int slice = ncol * ndep;
    int strip_tb = nf1;
    int slice_tb = NZ * nf1; 
    int strip_bt = nf2;
    int slice_bt = NZ * nf2;
    int strip_lr = nf1;
    int slice_lr = NZ * nf1;
    int strip_rl = nf2;
    int slice_rl = NZ * nf2;
    int strip_fb = nf1;
    int slice_fb = NX * nf1;
    int strip_bf = nf2;
    int slice_bf = NX * nf2;

	/*if (LOG){
	if (MYID==0) time1=MPI_Wtime();}*/

	/* top-bottom -----------------------------------------------------------*/


	if (POS[2]!=0) {	// no boundary exchange at top of global grid
		int n = 1;
        for (l=1;l<=FDORDER/2-1;l++) {
		    for (i=1;i<=NX;i++) {
			    // storage of top of local volume into buffer
			    for (k=1;k<=NZ;k++) {
			    	idx_v =  l * slice + i * strip + k;
			    	Ftb[i * slice_tb + k * strip_tb + n + 0] = Fv[3 * idx_v + 0];
			    	Ftb[i * slice_tb + k * strip_tb + n + 1] = Fv[3 * idx_v + 1];
			    	Ftb[i * slice_tb + k * strip_tb + n + 2] = Fv[3 * idx_v + 2];
			    }
		    }
		    n = n + 3;
		}
		for (l = FDORDER/2; l<=(FDORDER/2);l++) {
		    for (i=1;i<=NX;i++) {
			    // storage of top of local volume into buffer
			    for (k=1;k<=NZ;k++) {
			    	idx_v =  l * slice + i * strip + k;
			    	Ftb[i * slice_tb + k * strip_tb + n + 0] = Fv[3 * idx_v + 0];
			    	Ftb[i * slice_tb + k * strip_tb + n + 1] = Fv[3 * idx_v + 2];
			    }
		    }
		    n = n + 2;
		}
	}


    // if (POS[2] != NPROCY - 1) {	// no boundary exchange at top of global grid
	// 	int n = 1;
	// 	for (l=FDORDER/2-1;l>=1;l--) {
	// 	    for (i=1;i<=NX;i++) {
	// 			// storage of top of local volume into buffer
	// 			for (k=1;k<=NZ;k++) {
	// 		    	idx_v =  l * slice + i * strip + k;
	// 		    	Fbt[i * slice_bt + k * strip_bt + n + 0] = Fv[3 * (idx_v + (NY -2 * l  + 1) * slice) + 0];
	// 		    	Fbt[i * slice_bt + k * strip_bt + n + 1] = Fv[3 * (idx_v + (NY -2 * l  + 1) * slice) + 1];
	// 		    	Fbt[i * slice_bt + k * strip_bt + n + 2] = Fv[3 * (idx_v + (NY -2 * l  + 1) * slice) + 2];
	// 			}
	// 		}
	// 	    n = n + 3;
	// 	}
	// 	for (l = FDORDER/2; l>=(FDORDER/2);l--) {
	// 		for (i=1;i<=NX;i++) {
	// 		    // storage of top of local volume into buffer
	// 			for (k=1;k<=NZ;k++) {
	// 		    	idx_v =  l * slice + i * strip + k;
	// 		    	Fbt[i * slice_bt + k * strip_bt + n + 0] = Fv[3 * (idx_v + (NY -2 * l  + 1) * slice) + 1];
	// 			}
	// 		}
	// 		n = n + 1;
	// 	}
	// }


	//MPI_Sendrecv_replace(&buffertop_to_bot[1][1][1],nf1*NX*NZ,MPI_FLOAT,INDEX[3],TAG5,INDEX[4],TAG5,MPI_COMM_WORLD,&status);
    //	MPI_Sendrecv_replace(&bufferbot_to_top[1][1][1], NX * NZ * nf2, MPI_FLOAT, INDEX[4], TAG6, INDEX[3], TAG6, MPI_COMM_WORLD, &status);
	MPI_Sendrecv_replace(&Ftb[NZ * nf1 + nf1 + 1],nf1*NX*NZ,MPI_FLOAT,INDEX[3],TAG5,INDEX[4],TAG5,MPI_COMM_WORLD,&status);
	// MPI_Sendrecv_replace(&Fbt[NZ * nf2 + nf2 + 1],nf2 * NX * NZ, MPI_FLOAT, INDEX[4], TAG6, INDEX[3], TAG6, MPI_COMM_WORLD, &status);



	if (POS[2]!=NPROCY-1) {	// no boundary exchange at top of global grid
		int n = 1;
        for (l=1;l<=FDORDER/2-1;l++) {
		    for (i=1;i<=NX;i++) {
			    // storage of top of local volume into buffer
			    for (k=1;k<=NZ;k++) {
			    	idx_v =  l * slice + i * strip + k;
			    	Fv[3 * (idx_v + NY * slice) + 0] = Ftb[i * slice_tb + k * strip_tb + n + 0];
			    	Fv[3 * (idx_v + NY * slice) + 1] = Ftb[i * slice_tb + k * strip_tb + n + 1];
			    	Fv[3 * (idx_v + NY * slice) + 2] = Ftb[i * slice_tb + k * strip_tb + n + 2];
			    }
		    }
		    n = n + 3;
		}
		for (l = FDORDER/2; l<=(FDORDER/2);l++) {
		    for (i=1;i<=NX;i++) {
			    // storage of top of local volume into buffer
			    for (k=1;k<=NZ;k++) {
			    	idx_v =  l * slice + i * strip + k;
			    	Fv[3 * (idx_v + NY * slice) + 0] = Ftb[i * slice_tb + k * strip_tb + n + 0];
			    	Fv[3 * (idx_v + NY * slice) + 2] = Ftb[i * slice_tb + k * strip_tb + n + 1];
			    	//if(MYID == 0)
			    	//printf("--##-------------------->funcname:%s, lineNum:%d, l:%d, i:%d, k, n:%d:%d\n", __FUNCTION__, __LINE__, l, i, k,n);
			    }
		    }
		    n = n + 2;
		}

	}


	// if (POS[2]!=0) {	// no boundary exchange at top of global grid
	// 	int n = 1;
    //     for (l=FDORDER/2-1;l>=1;l--) {
	// 	    for (i=1;i<=NX;i++) {
	// 		    // storage of top of local volume into buffer
	// 		    for (k=1;k<=NZ;k++) {
	// 		    	idx_v =  l * slice + i * strip + k;
	// 		    	Fv[3 * (idx_v + (1 -2 * l) * slice) + 0] = Fbt[i * slice_bt + k * strip_bt + n + 0];
	// 		    	Fv[3 * (idx_v + (1 -2 * l) * slice) + 1] = Fbt[i * slice_bt + k * strip_bt + n + 1];
	// 		    	Fv[3 * (idx_v + (1 -2 * l) * slice) + 2] = Fbt[i * slice_bt + k * strip_bt + n + 2];
	// 		    }
	// 	    }
	// 	    n = n + 3;
	// 	}
	// 	for (l = FDORDER/2; l>=(FDORDER/2);l--) {
	// 	    for (i=1;i<=NX;i++) {
	// 		    // storage of top of local volume into buffer
	// 		    for (k=1;k<=NZ;k++) {
	// 		    	idx_v =  l * slice + i * strip + k;
	// 		    	Fv[3 * (idx_v + (1 -2 * l) * slice) + 1] = Fbt[i * slice_bt + k * strip_bt + n + 0];
	// 		    }
	// 	    }
	// 	    n = n + 1;
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
			    	idx_v =  j * slice + l * strip + k;
			    	Flr[j * slice_lr + k * strip_lr + n + 0] = Fv[3 * idx_v + 0];
			    	Flr[j * slice_lr + k * strip_lr + n + 1] = Fv[3 * idx_v + 1];
			    	Flr[j * slice_lr + k * strip_lr + n + 2] = Fv[3 * idx_v + 2];
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
			    	idx_v =  j * slice + l * strip + k;
			    	Flr[j * slice_lr + k * strip_lr + n + 0] = Fv[3 * idx_v + 1];
			    	Flr[j * slice_lr + k * strip_lr + n + 1] = Fv[3 * idx_v + 2];
			    }

			    n += 2;
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
			    	idx_v =  j * slice + l * strip + k;
			    	Frl[j * slice_rl + k * strip_rl + n + 0] = Fv[3 * (idx_v + (NX -2 * l  + 1) * strip) + 0];
			    	Frl[j * slice_rl + k * strip_rl + n + 1] = Fv[3 * (idx_v + (NX -2 * l  + 1) * strip) + 1];
			    	Frl[j * slice_rl + k * strip_rl + n + 2] = Fv[3 * (idx_v + (NX -2 * l  + 1) * strip) + 2];
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
			    	idx_v =  j * slice + l * strip + k;
			    	Frl[j * slice_rl + k * strip_rl + n + 0] = Fv[3 * (idx_v + (NX -2 * l  + 1) * strip) + 0];
			    }

			    n += 1;
		    }
	    }
	}


	//MPI_Sendrecv_replace(&bufferlef_to_rig[1][1][1],NY*NZ*nf1,MPI_FLOAT,INDEX[1],TAG1,INDEX[2],TAG1,MPI_COMM_WORLD,&status);
	//MPI_Sendrecv_replace(&bufferrig_to_lef[1][1][1],NY*NZ*nf2,MPI_FLOAT,INDEX[2],TAG2,INDEX[1],TAG2,MPI_COMM_WORLD,&status);
	MPI_Sendrecv_replace(&Flr[NZ * nf1 + nf1 + 1],ny*NZ*nf1,MPI_FLOAT,INDEX[1],TAG1,INDEX[2],TAG1,MPI_COMM_WORLD,&status);
	MPI_Sendrecv_replace(&Frl[NZ * nf2 + nf2 + 1],ny*NZ*nf2,MPI_FLOAT,INDEX[2],TAG2,INDEX[1],TAG2,MPI_COMM_WORLD,&status);






	if ((BOUNDARY) || (POS[1]!=NPROCX-1)) {	// no boundary exchange at left edge of global grid
		int n = 1;
	    for (j=1;j<=ny;j++){
	    	n = 1;
		    for (l=1;l<=(FDORDER/2-1);l++){
			    // storage of left edge of local volume into buffer
			    for (k=1;k<=NZ;k++) {
			    	idx_v =  j * slice + l * strip + k;
			    	Fv[3 * (idx_v + NX * strip) + 0] = Flr[j * slice_lr + k * strip_lr + n + 0];
			    	Fv[3 * (idx_v + NX * strip) + 1] = Flr[j * slice_lr + k * strip_lr + n + 1];
			    	Fv[3 * (idx_v + NX * strip) + 2] = Flr[j * slice_lr + k * strip_lr + n + 2];
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
			    	idx_v =  j * slice + l * strip + k;
			    	Fv[3 * (idx_v + NX * strip) + 1] = Flr[j * slice_lr + k * strip_lr + n + 0];
			    	Fv[3 * (idx_v + NX * strip) + 2] = Flr[j * slice_lr + k * strip_lr + n + 1];
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
			    	idx_v =  j * slice + l * strip + k;
			    	Fv[3 * (idx_v + (1 -2 * l) * strip) + 0] = Frl[j * slice_rl + k * strip_rl + n + 0];
			    	Fv[3 * (idx_v + (1 -2 * l) * strip) + 1] = Frl[j * slice_rl + k * strip_rl + n + 1];
			    	Fv[3 * (idx_v + (1 -2 * l) * strip) + 2] = Frl[j * slice_rl + k * strip_rl + n + 2];
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
			    	idx_v =  j * slice + l * strip + k;
			    	Fv[3 * (idx_v + (1 -2 * l) * strip) + 0] = Frl[j * slice_rl + k * strip_rl + n + 0];
			    }

			    n += 1;
		    }
	    }
	}







	/* front-back -----------------------------------------------------------*/


	if ((BOUNDARY) || (POS[3]!=0))	//* no boundary exchange at front side of global grid
		for (j=1;j<=ny;j++) {
			for (i=1;i<=NX;i++) {   //* storage of front side of local volume into buffer
				n=1;
			    for (l=1;l<=FDORDER/2;l++){
			    	idx_v =  j * slice + i * strip + l;
			    	Ffb[j * slice_fb + i * strip_fb + (n++)] = Fv[3 * idx_v + 0];
			    	Ffb[j * slice_fb + i * strip_fb + (n++)] = Fv[3 * idx_v + 1];
			    }
			    for (l=1;l<=(FDORDER/2-1);l++) {
			    	idx_v =  j * slice + i * strip + l;
			    	Ffb[j * slice_fb + i * strip_fb + (n++)] = Fv[3 * idx_v + 2];
			    }
			}
		}


	// no exchange if periodic boundary condition is applied


	if ((BOUNDARY) || (POS[3]!=NPROCZ-1))	//* no boundary exchange at back side of global grid
		for (j=1;j<=ny;j++) {
			for (i=1;i<=NX;i++) {
					//* storage of back side of local volume into buffer
				n=1;
				for (l=FDORDER/2;l>=1;l--) {
					idx_v =  j * slice + i * strip + l;
					Fbf[j * slice_bf + i * strip_bf + (n++)] = Fv[3 * (idx_v + (NZ -2 * l  + 1)) + 2];
				}

				for (l=(FDORDER/2-1);l>=1;l--) {
					idx_v =  j * slice + i * strip + l;
                    Fbf[j * slice_bf + i * strip_bf + (n++)] = Fv[3 * (idx_v + (NZ -2 * l  + 1)) + 0];
                    Fbf[j * slice_bf + i * strip_bf + (n++)] = Fv[3 * (idx_v + (NZ -2 * l  + 1)) + 1];
				}

			}
		}



	//MPI_Sendrecv_replace(&bufferfro_to_bac[1][1][1],NX*NY*nf1,MPI_FLOAT,INDEX[5],TAG3,INDEX[6],TAG3,MPI_COMM_WORLD,&status);
	//MPI_Sendrecv_replace(&bufferbac_to_fro[1][1][1],NX*NY*nf2,MPI_FLOAT,INDEX[6],TAG4,INDEX[5],TAG4,MPI_COMM_WORLD,&status);
	MPI_Sendrecv_replace(&Ffb[NX * nf1 + nf1 + 1],NX*ny*nf1,MPI_FLOAT,INDEX[5],TAG3,INDEX[6],TAG3,MPI_COMM_WORLD,&status);
	MPI_Sendrecv_replace(&Fbf[NX * nf2 + nf2 + 1],NX*ny*nf2,MPI_FLOAT,INDEX[6],TAG4,INDEX[5],TAG4,MPI_COMM_WORLD,&status);


	/* no exchange if periodic boundary condition is applied */

	if ((BOUNDARY) || (POS[3]!=NPROCZ-1))  //* no boundary exchange at back side of global grid
		for (j=1;j<=ny;j++) {
			for (i=1;i<=NX;i++) {
				n=1;
				for (l=1;l<=FDORDER/2;l++) {
					idx_v =  j * slice + i * strip + l;
					Fv[3 * (idx_v + NZ) + 0] = Ffb[j * slice_fb + i * strip_fb + (n++)];
					Fv[3 * (idx_v + NZ) + 1] = Ffb[j * slice_fb + i * strip_fb + (n++)];
				}

				for (l=1;l<=(FDORDER/2-1);l++) {
					idx_v =  j * slice + i * strip + l;
					Fv[3 * (idx_v + NZ) + 2] = Ffb[j * slice_fb + i * strip_fb + (n++)];
				}
			}
		}



	/* no exchange if periodic boundary condition is applied */
	if ((BOUNDARY) || (POS[3]!=0))	/* no boundary exchange at front side of global grid */
		for (j=1;j<=ny;j++) {
			for (i=1;i<=NX;i++) {
				n=1;
				for (l=FDORDER/2;l>=1;l--) {
					idx_v =  j * slice + i * strip + l;
					Fv[3 * (idx_v + (1 -2 * l)) + 2] = Fbf[j * slice_bf + i * strip_bf + (n++)];
				}

				for (l=(FDORDER/2-1);l>=1;l--) {
					idx_v =  j * slice + i * strip + l;
					Fv[3 * (idx_v + (1 -2 * l)) + 0] = Fbf[j * slice_bf + i * strip_bf + (n++)];
					Fv[3 * (idx_v + (1 -2 * l)) + 1] = Fbf[j * slice_bf + i * strip_bf + (n++)];
				}
			}
		}



	free_trans(Ftb, 1, NX, 1, NZ, 1, nf1);
	// free_trans(Fbt, 1, NX, 1, NZ, 1, nf2);
	// free_trans(Flr, 1, NY, 1, NZ, 1, nf1);
	// free_trans(Frl, 1, NY, 1, NZ, 1, nf2);
	// free_trans(Ffb, 1, NY, 1, NX, 1, nf1);
	// free_trans(Fbf, 1, NY, 1, NX, 1, nf2);

	return time;

}




double exchange_v2_2(float *Fv, float *** bufferbot_to_top, MPI_Request * req_send, MPI_Request * req_rec) {

	extern int NX, NY, NZ, POS[4], NPROCX, NPROCY, NPROCZ, BOUNDARY,  FDORDER,  INDEX[7], ABS_TYPE;  extern int MYID; /*MYID,LOG,*/
	extern const int TAG1,TAG2,TAG3,TAG4,TAG5,TAG6;

	MPI_Status status;
	int i, j, k, l, n, nf1, nf2, const_n;
	double time=0.0;  /*, time1=0.0;*/

    nf1=3*FDORDER/2-1;
	nf2=nf1-1;
    //printf("--##-------------------->funcname:%s, lineNum:%d\n", __FUNCTION__, __LINE__);
	int idx_v;
	int ll = 1;
	if(ABS_TYPE==1 && FDORDER==2){ll=2;}
	int nrow = NY + 2 * (ll*FDORDER/2) + 1, ncol = NX + 2 * (ll*FDORDER/2), ndep = NZ + 2 * (ll*FDORDER/2);
	int size = nrow * ncol * ndep;
	// float *Ftb = transform3(buffertop_to_bot, 1, NX, 1, NZ,  1, nf1);
	float *Fbt = transform3(bufferbot_to_top, 1, NX, 1, NZ,  1, nf2);
	int ny = NY/2;
	int ny1 = ((NY % 2 == 0) ? (NY/2):(NY/2+1));
	float Flr[(ny+1)*(NZ+1)*(nf1+1)],Frl[(ny+1)*(NZ+1)*(nf2+1)],Ffb[(ny+1)*(NX+1)*(nf1+1)],Fbf[(ny+1)*(NX+1)*(nf2+1)];
    int strip = ndep;
    int slice = ncol * ndep;
    int strip_tb = nf1;
    int slice_tb = NZ * nf1; 
    int strip_bt = nf2;
    int slice_bt = NZ * nf2;
    int strip_lr = nf1;
    int slice_lr = NZ * nf1;
    int strip_rl = nf2;
    int slice_rl = NZ * nf2;
    int strip_fb = nf1;
    int slice_fb = NX * nf1;
    int strip_bf = nf2;
    int slice_bf = NX * nf2;

	/*if (LOG){
	if (MYID==0) time1=MPI_Wtime();}*/

	/* top-bottom -----------------------------------------------------------*/




    if (POS[2] != NPROCY - 1) {	// no boundary exchange at top of global grid
		int n = 1;
		for (l=FDORDER/2-1;l>=1;l--) {
		    for (i=1;i<=NX;i++) {
				// storage of top of local volume into buffer
				for (k=1;k<=NZ;k++) {
			    	idx_v =  l * slice + i * strip + k;
			    	Fbt[i * slice_bt + k * strip_bt + n + 0] = Fv[3 * (idx_v + (NY -2 * l  + 1) * slice) + 0];
			    	Fbt[i * slice_bt + k * strip_bt + n + 1] = Fv[3 * (idx_v + (NY -2 * l  + 1) * slice) + 1];
			    	Fbt[i * slice_bt + k * strip_bt + n + 2] = Fv[3 * (idx_v + (NY -2 * l  + 1) * slice) + 2];
				}
			}
		    n = n + 3;
		}
		for (l = FDORDER/2; l>=(FDORDER/2);l--) {
			for (i=1;i<=NX;i++) {
			    // storage of top of local volume into buffer
				for (k=1;k<=NZ;k++) {
			    	idx_v =  l * slice + i * strip + k;
			    	Fbt[i * slice_bt + k * strip_bt + n + 0] = Fv[3 * (idx_v + (NY -2 * l  + 1) * slice) + 1];
				}
			}
			n = n + 1;
		}
	}


	//MPI_Sendrecv_replace(&buffertop_to_bot[1][1][1],nf1*NX*NZ,MPI_FLOAT,INDEX[3],TAG5,INDEX[4],TAG5,MPI_COMM_WORLD,&status);
    //	MPI_Sendrecv_replace(&bufferbot_to_top[1][1][1], NX * NZ * nf2, MPI_FLOAT, INDEX[4], TAG6, INDEX[3], TAG6, MPI_COMM_WORLD, &status);
	//MPI_Sendrecv_replace(&Ftb[NZ * nf1 + nf1 + 1],nf1*NX*NZ,MPI_FLOAT,INDEX[3],TAG5,INDEX[4],TAG5,MPI_COMM_WORLD,&status);
	MPI_Sendrecv_replace(&Fbt[NZ * nf2 + nf2 + 1],nf2 * NX * NZ, MPI_FLOAT, INDEX[4], TAG6, INDEX[3], TAG6, MPI_COMM_WORLD, &status);




	if (POS[2]!=0) {	// no boundary exchange at top of global grid
		int n = 1;
        for (l=FDORDER/2-1;l>=1;l--) {
		    for (i=1;i<=NX;i++) {
			    // storage of top of local volume into buffer
			    for (k=1;k<=NZ;k++) {
			    	idx_v =  l * slice + i * strip + k;
			    	Fv[3 * (idx_v + (1 -2 * l) * slice) + 0] = Fbt[i * slice_bt + k * strip_bt + n + 0];
			    	Fv[3 * (idx_v + (1 -2 * l) * slice) + 1] = Fbt[i * slice_bt + k * strip_bt + n + 1];
			    	Fv[3 * (idx_v + (1 -2 * l) * slice) + 2] = Fbt[i * slice_bt + k * strip_bt + n + 2];
			    }
		    }
		    n = n + 3;
		}
		for (l = FDORDER/2; l>=(FDORDER/2);l--) {
		    for (i=1;i<=NX;i++) {
			    // storage of top of local volume into buffer
			    for (k=1;k<=NZ;k++) {
			    	idx_v =  l * slice + i * strip + k;
			    	Fv[3 * (idx_v + (1 -2 * l) * slice) + 1] = Fbt[i * slice_bt + k * strip_bt + n + 0];
			    }
		    }
		    n = n + 1;
		}

	}





	/* left-right -----------------------------------------------------------*/



	if ((BOUNDARY) || (POS[1]!=0)) {	// no boundary exchange at left edge of global grid
		int n = 1;
	    for (j=ny1+1;j<=NY;j++){
	    	n = 1;
		    for (l=1;l<=(FDORDER/2-1);l++){
			    // storage of left edge of local volume into buffer
			    for (k=1;k<=NZ;k++) {
			    	idx_v =  j * slice + l * strip + k;
			    	Flr[(j - ny1) * slice_lr + k * strip_lr + n + 0] = Fv[3 * idx_v + 0];
			    	Flr[(j - ny1) * slice_lr + k * strip_lr + n + 1] = Fv[3 * idx_v + 1];
			    	Flr[(j - ny1) * slice_lr + k * strip_lr + n + 2] = Fv[3 * idx_v + 2];
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
			    	idx_v =  j * slice + l * strip + k;
			    	Flr[(j - ny1) * slice_lr + k * strip_lr + n + 0] = Fv[3 * idx_v + 1];
			    	Flr[(j - ny1) * slice_lr + k * strip_lr + n + 1] = Fv[3 * idx_v + 2];
			    }

			    n += 2;
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
			    	idx_v =  j * slice + l * strip + k;
			    	Frl[(j - ny1) * slice_rl + k * strip_rl + n + 0] = Fv[3 * (idx_v + (NX -2 * l  + 1) * strip) + 0];
			    	Frl[(j - ny1) * slice_rl + k * strip_rl + n + 1] = Fv[3 * (idx_v + (NX -2 * l  + 1) * strip) + 1];
			    	Frl[(j - ny1) * slice_rl + k * strip_rl + n + 2] = Fv[3 * (idx_v + (NX -2 * l  + 1) * strip) + 2];
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
			    	idx_v =  j * slice + l * strip + k;
			    	Frl[(j - ny1)* slice_rl + k * strip_rl + n + 0] = Fv[3 * (idx_v + (NX -2 * l  + 1) * strip) + 0];
			    }

			    n += 1;
		    }
	    }
	}


	//MPI_Sendrecv_replace(&bufferlef_to_rig[1][1][1],NY*NZ*nf1,MPI_FLOAT,INDEX[1],TAG1,INDEX[2],TAG1,MPI_COMM_WORLD,&status);
	//MPI_Sendrecv_replace(&bufferrig_to_lef[1][1][1],NY*NZ*nf2,MPI_FLOAT,INDEX[2],TAG2,INDEX[1],TAG2,MPI_COMM_WORLD,&status);
	MPI_Sendrecv_replace(&Flr[NZ * nf1 + nf1 + 1],ny*NZ*nf1,MPI_FLOAT,INDEX[1],TAG1,INDEX[2],TAG1,MPI_COMM_WORLD,&status);
	MPI_Sendrecv_replace(&Frl[NZ * nf2 + nf2 + 1],ny*NZ*nf2,MPI_FLOAT,INDEX[2],TAG2,INDEX[1],TAG2,MPI_COMM_WORLD,&status);






	if ((BOUNDARY) || (POS[1]!=NPROCX-1)) {	// no boundary exchange at left edge of global grid
		int n = 1;
	    for (j=ny1+1;j<=NY;j++){
	    	n = 1;
		    for (l=1;l<=(FDORDER/2-1);l++){
			    // storage of left edge of local volume into buffer
			    for (k=1;k<=NZ;k++) {
			    	idx_v =  j * slice + l * strip + k;
			    	Fv[3 * (idx_v + NX * strip) + 0] = Flr[(j - ny1) * slice_lr + k * strip_lr + n + 0];
			    	Fv[3 * (idx_v + NX * strip) + 1] = Flr[(j - ny1) * slice_lr + k * strip_lr + n + 1];
			    	Fv[3 * (idx_v + NX * strip) + 2] = Flr[(j - ny1) * slice_lr + k * strip_lr + n + 2];
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
			    	idx_v =  j * slice + l * strip + k;
			    	Fv[3 * (idx_v + NX * strip) + 1] = Flr[(j - ny1) * slice_lr + k * strip_lr + n + 0];
			    	Fv[3 * (idx_v + NX * strip) + 2] = Flr[(j - ny1) * slice_lr + k * strip_lr + n + 1];
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
			    	idx_v =  j * slice + l * strip + k;
			    	Fv[3 * (idx_v + (1 -2 * l) * strip) + 0] = Frl[(j - ny1) * slice_rl + k * strip_rl + n + 0];
			    	Fv[3 * (idx_v + (1 -2 * l) * strip) + 1] = Frl[(j - ny1) * slice_rl + k * strip_rl + n + 1];
			    	Fv[3 * (idx_v + (1 -2 * l) * strip) + 2] = Frl[(j - ny1) * slice_rl + k * strip_rl + n + 2];
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
			    	idx_v =  j * slice + l * strip + k;
			    	Fv[3 * (idx_v + (1 -2 * l) * strip) + 0] = Frl[(j - ny1) * slice_rl + k * strip_rl + n + 0];
			    }

			    n += 1;
		    }
	    }
	}







	/* front-back -----------------------------------------------------------*/


	if ((BOUNDARY) || (POS[3]!=0))	//* no boundary exchange at front side of global grid
		for (j=ny1+1;j<=NY;j++) {
			for (i=1;i<=NX;i++) {   //* storage of front side of local volume into buffer
				n=1;
			    for (l=1;l<=FDORDER/2;l++){
			    	idx_v =  j * slice + i * strip + l;
			    	Ffb[(j - ny1) * slice_fb + i * strip_fb + (n++)] = Fv[3 * idx_v + 0];
			    	Ffb[(j - ny1) * slice_fb + i * strip_fb + (n++)] = Fv[3 * idx_v + 1];
			    }
			    for (l=1;l<=(FDORDER/2-1);l++) {
			    	idx_v =  j * slice + i * strip + l;
			    	Ffb[(j - ny1) * slice_fb + i * strip_fb + (n++)] = Fv[3 * idx_v + 2];
			    }
			}
		}


	// no exchange if periodic boundary condition is applied


	if ((BOUNDARY) || (POS[3]!=NPROCZ-1))	//* no boundary exchange at back side of global grid
		for (j=ny1+1;j<=NY;j++) {
			for (i=1;i<=NX;i++) {
					//* storage of back side of local volume into buffer
				n=1;
				for (l=FDORDER/2;l>=1;l--) {
					idx_v =  j * slice + i * strip + l;
					Fbf[(j - ny1) * slice_bf + i * strip_bf + (n++)] = Fv[3 * (idx_v + (NZ -2 * l  + 1)) + 2];
				}

				for (l=(FDORDER/2-1);l>=1;l--) {
					idx_v =  j * slice + i * strip + l;
                    Fbf[(j - ny1) * slice_bf + i * strip_bf + (n++)] = Fv[3 * (idx_v + (NZ -2 * l  + 1)) + 0];
                    Fbf[(j - ny1) * slice_bf + i * strip_bf + (n++)] = Fv[3 * (idx_v + (NZ -2 * l  + 1)) + 1];
				}

			}
		}



	//MPI_Sendrecv_replace(&bufferfro_to_bac[1][1][1],NX*NY*nf1,MPI_FLOAT,INDEX[5],TAG3,INDEX[6],TAG3,MPI_COMM_WORLD,&status);
	//MPI_Sendrecv_replace(&bufferbac_to_fro[1][1][1],NX*NY*nf2,MPI_FLOAT,INDEX[6],TAG4,INDEX[5],TAG4,MPI_COMM_WORLD,&status);
	MPI_Sendrecv_replace(&Ffb[NX * nf1 + nf1 + 1],NX*ny*nf1,MPI_FLOAT,INDEX[5],TAG3,INDEX[6],TAG3,MPI_COMM_WORLD,&status);
	MPI_Sendrecv_replace(&Fbf[NX * nf2 + nf2 + 1],NX*ny*nf2,MPI_FLOAT,INDEX[6],TAG4,INDEX[5],TAG4,MPI_COMM_WORLD,&status);


	/* no exchange if periodic boundary condition is applied */

	if ((BOUNDARY) || (POS[3]!=NPROCZ-1))  //* no boundary exchange at back side of global grid
		for (j=ny1+1;j<=NY;j++) {
			for (i=1;i<=NX;i++) {
				n=1;
				for (l=1;l<=FDORDER/2;l++) {
					idx_v =  j * slice + i * strip + l;
					Fv[3 * (idx_v + NZ) + 0] = Ffb[(j - ny1) * slice_fb + i * strip_fb + (n++)];
					Fv[3 * (idx_v + NZ) + 1] = Ffb[(j - ny1) * slice_fb + i * strip_fb + (n++)];
				}

				for (l=1;l<=(FDORDER/2-1);l++) {
					idx_v =  j * slice + i * strip + l;
					Fv[3 * (idx_v + NZ) + 2] = Ffb[(j - ny1) * slice_fb + i * strip_fb + (n++)];
				}
			}
		}



	/* no exchange if periodic boundary condition is applied */
	if ((BOUNDARY) || (POS[3]!=0))	/* no boundary exchange at front side of global grid */
		for (j=ny1+1;j<=NY;j++) {
			for (i=1;i<=NX;i++) {
				n=1;
				for (l=FDORDER/2;l>=1;l--) {
					idx_v =  j * slice + i * strip + l;
					Fv[3 * (idx_v + (1 -2 * l)) + 2] = Fbf[(j - ny1) * slice_bf + i * strip_bf + (n++)];
				}

				for (l=(FDORDER/2-1);l>=1;l--) {
					idx_v =  j * slice + i * strip + l;
					Fv[3 * (idx_v + (1 -2 * l)) + 0] = Fbf[(j - ny1) * slice_bf + i * strip_bf + (n++)];
					Fv[3 * (idx_v + (1 -2 * l)) + 1] = Fbf[(j - ny1) * slice_bf + i * strip_bf + (n++)];
				}
			}
		}



	free_trans(Fbt, 1, NX, 1, NZ, 1, nf2);


	return time;

}



