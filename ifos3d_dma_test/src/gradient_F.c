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
 * gradient calculation in frequency domain:
 * gradient as multiplication of forward and conjugate backpropagated wavefield 
 * spatial derivatives are calculated by 4th order finite differences
 * S. Butzer 2013
--------------------------------------------------------------------------- */


#include "fd.h"
#include <signal.h>
#include <slave.h>
extern SLAVE_FUN(gradient_kernel_slave)(Param_grad *);


inline long get_cpe_pc(int cpe_id) {
	athread_task_info();
    char cgn = athread_get_cgn();
    int _CGN =(int) cgn;
	return *(long*)(549806153728L + _CGN * 68719476736L + cpe_id * 65536L);
}

void pc_on_sig(int sig) {
	printf("writing cpe PCs on signal %d\n", sig);
	int i, j;
	for (i = 0; i < 8; i ++) {
		for (j = 0; j < 8; j ++)
			printf("%3d: %12llx ", i * 8 + j, get_cpe_pc(i * 8 + j));
		puts("");

	}
}




void gradient_F(int nx, int ny, int nz, float *F_fv, float *F_bv, float ***grad1, float ***grad2,float ***grad3,int nt, float  ***  rho, float ***  pi, float ***  u, float * finv, int nf, int ntr_hess, int iteration) {

	extern float DX, DY, DZ, DT, MYID;
	extern int POS[4], FDCOEFF;
	extern char  MFILE[STRING_SIZE];


	/*float vp0=6200.0, vs0=3600.0, rho0=2800.0;*/

	char gradfile1[STRING_SIZE],gradfile2[STRING_SIZE],gradfile3[STRING_SIZE],gradfile4[STRING_SIZE],gradfile5[STRING_SIZE],gradfile6[STRING_SIZE],gradfile7[STRING_SIZE];
	
	/*FILE *fpmod1, *fpmod2, *fpmod3,*fpmod4, *fpmod5, *fpmod6,*fpmod7;*/
	



	extern int NX, NY, NZ, FDORDER, NFMAX, ABS_TYPE;
	int ll = 1;
	if(ABS_TYPE==1 && FDORDER==2){ll=2;}


	//int nval = NFMAX;
	int nrow = NY + 2 * (ll*FDORDER/2) + 1, ncol = NX + 2 * (ll*FDORDER/2), ndep = NZ + 2 * (ll*FDORDER/2);
	//int nval_b = NFMAX*(ntr_hess+1);
    float *Fpi = transform3(pi, 1, NY + 1, 1, NX + 1, 1, NZ + 1);
    float *Fu = transform3(u, 1, NY + 1, 1, NX + 1, 1, NZ + 1);
    float *Frho = transform3(rho, 1, NY + 1, 1, NX + 1, 1, NZ + 1);
    float *Fgrad = transform_3(grad1, grad2, grad3, 1, NY, 1, NX, 1, NZ);
    int strip = ndep;
    int slice = ncol * ndep;
    int cube = nrow * ncol * ndep;
    int strip_grad = NZ;
    int slice_grad = NZ * NX;
	int strip_pi = NZ + 1;
	int slice_pi = (NX + 1) * (NZ + 1);


		
	sprintf(gradfile1,"%s.grad1.%i%i%i",MFILE,POS[1],POS[2],POS[3]);
	sprintf(gradfile2,"%s.grad2.%i%i%i",MFILE,POS[1],POS[2],POS[3]);
	sprintf(gradfile3,"%s.grad3.%i%i%i",MFILE,POS[1],POS[2],POS[3]);
	sprintf(gradfile4,"%s.grad4.%i%i%i",MFILE,POS[1],POS[2],POS[3]);
	sprintf(gradfile5,"%s.grad5.%i%i%i",MFILE,POS[1],POS[2],POS[3]);
	sprintf(gradfile6,"%s.grad6.%i%i%i",MFILE,POS[1],POS[2],POS[3]);
	sprintf(gradfile7,"%s.grad7.%i%i%i",MFILE,POS[1],POS[2],POS[3]);
	
	/*fpmod1=fopen(gradfile1,"w");
	fpmod2=fopen(gradfile2,"w");
	fpmod3=fopen(gradfile3,"w");
	fpmod4=fopen(gradfile4,"w");
	fpmod5=fopen(gradfile5,"w");
	fpmod6=fopen(gradfile6,"w");
	fpmod7=fopen(gradfile7,"w");*/
	
	Param_grad param;

	float *fvel = F_fv + 6 * (slice + strip + 1);
	float *bvel = F_bv + 6 * (slice + strip + 1);
	float *grad = Fgrad + 3 * (slice_grad + strip_grad + 1);
	float *Pi = Fpi + (slice_pi + strip_pi + 1);
	float *U = Fu + (slice_pi + strip_pi + 1);
	float *Rho = Frho + (slice_pi + strip_pi + 1);

	param.dim_x = nx;
	param.dim_y = ny;
	param.dim_z = nz;
	param.cube = cube;
	param.slice = slice;
	param.strip = strip;
	param.slice_grad = slice_grad;
	param.strip_grad = strip_grad;
	param.slice_pi = slice_pi;
	param.strip_pi = strip_pi;
	param.DT = DT;
	param.DX = DX;
	param.DY = DY;
	param.DZ = DZ;
	param.FDCOEFF = FDCOEFF;
	param.nf = nf;
	param.Fgrad = Fgrad;
	param.F_fv = F_fv;
	param.F_bv = F_bv;
	param.Fpi = Fpi;
	param.Fu = Fu;
	param.Frho = Frho;
	param.finv = finv;
	param.MYID = MYID;
	param.Fgrad = Fgrad;


	static int init = 1;

	//signal(SIGUSR1, pc_on_sig);
	athread_spawn(gradient_kernel_slave, &param);
	athread_join();



	inverse_3(Fgrad, grad1, grad2, grad3, 1, NY, 1, NX, 1, NZ);

	free_trans_3(Fgrad, 1, NY, 1, NX, 1, NZ);
    free_trans(Fpi, 1, NY + 1, 1, NX + 1, 1, NZ + 1);
    free_trans(Fu, 1, NY + 1, 1, NX + 1, 1, NZ + 1);
    free_trans(Frho, 1, NY + 1, 1, NX + 1, 1, NZ + 1);
	

}
