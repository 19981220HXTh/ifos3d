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
 * calculation of diagonal Hessian approximation in frequency domain:
 * spatial derivatives are calculated by 4th order finite differences
 * S. Butzer 2013
--------------------------------------------------------------------------- */

#include "fd.h"


void hess_F(int nx,int ny,int nz,float *F_fv, float *F_bv, float ***hess1, float ***hess2,float ***hess3,int nt, float  ***  rho, float ***  pi, float ***  u, float * finv, int nf, int ntr_hess) {

	extern float DX, DY, DZ, MYID;
	extern int  FDCOEFF;
	/*extern char  MFILE[STRING_SIZE];*/

	printf("---###------------------------------------------------------> %s,  %d,  MYID: %d\n", __FUNCTION__, __LINE__, MYID);
		
	float fvxx=0.0,fvxy=0.0,fvxz=0.0,fvyx=0.0,fvyy=0.0,fvyz=0.0,fvzx=0.0,fvzy=0.0,fvzz=0.0;
	float bvxx=0.0,bvxy=0.0,bvxz=0.0,bvyx=0.0,bvyy=0.0,bvyz=0.0,bvzx=0.0,bvzy=0.0,bvzz=0.0;
	float fivxx=0.0,fivxy=0.0,fivxz=0.0,fivyx=0.0,fivyy=0.0,fivyz=0.0,fivzx=0.0,fivzy=0.0,fivzz=0.0;
	float bivxx=0.0,bivxy=0.0,bivxz=0.0,bivyx=0.0,bivyy=0.0,bivyz=0.0,bivzx=0.0,bivzy=0.0,bivzz=0.0;
	float relam, remu, rerho,imlam, immu, imrho,revs,imvs, rerho1, imrho1;
	/*float hessmu=0.0, hesslam=0.0, hessrho=0.0;*/
	float b1,b2,fdummy;
	/*float vp0=6200.0, vs0=3600.0, rho0=2800.0;*/
	
	int i,j,k,l,m,n,pos_f,pos_b,pos_pi;
	int pos_f_jp1,pos_f_jp2,pos_f_jm1,pos_f_jm2,pos_b_jp1,pos_b_jp2,pos_b_jm1,pos_b_jm2;
	int pos_f_ip1,pos_f_ip2,pos_f_im1,pos_f_im2,pos_b_ip1,pos_b_ip2,pos_b_im1,pos_b_im2;
	int pos_f_kp1,pos_f_kp2,pos_f_km1,pos_f_km2,pos_b_kp1,pos_b_kp2,pos_b_km1,pos_b_km2;

	extern int NX, NY, NZ, FDORDER, NFMAX, ABS_TYPE;
	int ll = 1;
	if(ABS_TYPE==1 && FDORDER==2){ll=2;}
	int nval = NFMAX, nrow = NY + 2 * (ll*FDORDER/2) + 1, ncol = NX + 2 * (ll*FDORDER/2), ndep = NZ + 2 * (ll*FDORDER/2);
	int nval_b = NFMAX*(ntr_hess+1);
    float *Fpi = transform3(pi, 1, NY + 1, 1, NX + 1, 1, NZ + 1);
    float *Fu = transform3(u, 1, NY + 1, 1, NX + 1, 1, NZ + 1);
    float *Frho = transform3(rho, 1, NY + 1, 1, NX + 1, 1, NZ + 1);
    int strip = ndep;
    int slice = ncol * ndep;
    int cube = nrow * ncol * ndep;
	int strip_pi = NZ + 1;
	int slice_pi = (NX + 1) * (NZ + 1);



	/*t=NT-nt+1;*/
	b1=9.0/8.0; b2=-1.0/24.0; /* Taylor coefficients (4th order)*/
		if(FDCOEFF==2){
		b1=1.1382; b2=-0.046414;} /* Holberg coefficients E=0.1 %*/	
	for(n=1;n<=ntr_hess;n++)
	for(l=1;l<=nf;l++){
	  m=(n)*nf+l;
	  fdummy=0.0;
	  fdummy=finv[l-1]*M_PI*2;
	  fdummy=1/fdummy;
		for (j=1;j<=ny;j++){
			for (i=1;i<=nx;i++){
				for (k=1;k<=nz;k++){

				// spatial derivatives of the components of the velocities are computed
				pos_f = l * cube +  j * slice + i * strip + k;
				pos_b = m * cube +  j * slice + i * strip + k;
				pos_pi = j * slice_pi + i * strip_pi + k;
				pos_b_jm2 = pos_b - 2 * slice;
				pos_b_jm1 = pos_b - slice;
				pos_b_jp2 = pos_b + 2 * slice;
				pos_b_jp1 = pos_b + slice;
				pos_f_jm2 = pos_f - 2 * slice;
				pos_f_jm1 = pos_f - slice;
				pos_f_jp2 = pos_f + 2 * slice;
				pos_f_jp1 = pos_f + slice;
				pos_b_im2 = pos_b - 2 * strip;
				pos_b_im1 = pos_b - strip;
				pos_b_ip2 = pos_b + 2 * strip;
				pos_b_ip1 = pos_b + strip;
				pos_f_im2 = pos_f - 2 * strip;
				pos_f_im1 = pos_f - strip;
				pos_f_ip2 = pos_f + 2 * strip;
				pos_f_ip1 = pos_f + strip;
				pos_b_km2 = pos_b - 2;
				pos_b_km1 = pos_b - 1;
				pos_b_kp2 = pos_b + 2;
				pos_b_kp1 = pos_b + 1;
				pos_f_km2 = pos_f - 2;
				pos_f_km1 = pos_f - 1;
				pos_f_kp2 = pos_f + 2;
				pos_f_kp1 = pos_f + 1;

				fvxx = (b1 * (F_fv[pos_f * 6] - F_fv[pos_f_im1 * 6]) + b2 * (F_fv[pos_f_ip1 * 6] - F_fv[pos_f_im2 * 6])) / DX;
				fvxy = (b1 * (F_fv[pos_f_jp1 * 6] - F_fv[pos_f * 6]) + b2 * (F_fv[pos_f_jp2 * 6] - F_fv[pos_f_jm1 * 6])) / DY;
				fvxz = (b1 * (F_fv[pos_f_kp1 * 6] - F_fv[pos_f * 6]) + b2 * (F_fv[pos_f_kp2 * 6] - F_fv[pos_f_km1 * 6])) / DZ;
				fvyx = (b1 * (F_fv[pos_f_ip1 * 6 + 1] - F_fv[pos_f * 6 + 1]) + b2 * (F_fv[pos_f_ip2 * 6 + 1] - F_fv[pos_f_im1 * 6 + 1])) / DX;
				fvyy = (b1 * (F_fv[pos_f * 6 + 1] - F_fv[pos_f_jm1 * 6 + 1]) + b2 * (F_fv[pos_f_jp1 * 6 + 1] - F_fv[pos_f_jm2 * 6 + 1])) / DY;
				fvyz = (b1 * (F_fv[pos_f_kp1 * 6 + 1] - F_fv[pos_f * 6 + 1]) + b2 * (F_fv[pos_f_kp2 * 6 + 1] - F_fv[pos_f_km1 * 6 + 1])) / DZ;
				fvzx = (b1 * (F_fv[pos_f_ip1 * 6 + 2] - F_fv[pos_f * 6 + 2]) + b2 * (F_fv[pos_f_ip2 * 6 + 2] - F_fv[pos_f_im1 * 6 + 2])) / DX;
				fvzy = (b1 * (F_fv[pos_f_jp1 * 6 + 2] - F_fv[pos_f * 6 + 2]) + b2 * (F_fv[pos_f_jp2 * 6 + 2] - F_fv[pos_f_jm1 * 6 + 2])) / DY;
				fvzz = (b1 * (F_fv[pos_f * 6 + 2] - F_fv[pos_f_km1 * 6 + 2]) + b2 * (F_fv[pos_f_kp1 * 6 + 2] - F_fv[pos_f_km2 * 6 + 2])) / DY;

				fivxx = (b1 * (F_fv[pos_f * 6 + 3] - F_fv[pos_f_im1 * 6 + 3]) + b2 * (F_fv[pos_f_ip1 * 6 + 3] - F_fv[pos_f_im2 * 6 + 3])) / DX;
				fivxy = (b1 * (F_fv[pos_f_jp1 * 6 + 3] - F_fv[pos_f * 6 + 3]) + b2 * (F_fv[pos_f_jp2 * 6 + 3] - F_fv[pos_f_jm1 * 6 + 3])) / DY;
				fivxz = (b1 * (F_fv[pos_f_kp1 * 6 + 3] - F_fv[pos_f * 6 + 3]) + b2 * (F_fv[pos_f_kp2 * 6 + 3] - F_fv[pos_f_km1 * 6 + 3])) / DZ;
				fivyx = (b1 * (F_fv[pos_f_ip1 * 6 + 4] - F_fv[pos_f * 6 + 4]) + b2 * (F_fv[pos_f_ip2 * 6 + 4] - F_fv[pos_f_im1 * 6 + 4])) / DX;
				fivyy = (b1 * (F_fv[pos_f * 6 + 4] - F_fv[pos_f_jm1 * 6 + 4]) + b2 * (F_fv[pos_f_jp1 * 6 + 4] - F_fv[pos_f_jm2 * 6 + 4])) / DY;
				fivyz = (b1 * (F_fv[pos_f_kp1 * 6 + 4] - F_fv[pos_f * 6 + 4]) + b2 * (F_fv[pos_f_kp2 * 6 + 4] - F_fv[pos_f_km1 * 6 + 4])) / DZ;
				fivzx = (b1 * (F_fv[pos_f_ip1 * 6 + 5] - F_fv[pos_f * 6 + 5]) + b2 * (F_fv[pos_f_ip2 * 6 + 5] - F_fv[pos_f_im1 * 6 + 5])) / DX;
				fivzy = (b1 * (F_fv[pos_f_jp1 * 6 + 5] - F_fv[pos_f * 6 + 5]) + b2 * (F_fv[pos_f_jp2 * 6 + 5] - F_fv[pos_f_jm1 * 6 + 5])) / DY;
				fivzz = (b1 * (F_fv[pos_f * 6 + 5] - F_fv[pos_f_km1 * 6 + 5]) + b2 * (F_fv[pos_f_kp1 * 6 + 5] - F_fv[pos_f_km2 * 6 + 5])) / DY;

				bvxx = (b1 * (F_bv[pos_b * 6] - F_bv[pos_b_im1 * 6]) + b2 * (F_bv[pos_b_ip1 * 6] - F_bv[pos_b_im2 * 6])) / DX;
				bvxy = (b1 * (F_bv[pos_b_jp1 * 6] - F_bv[pos_b * 6]) + b2 * (F_bv[pos_b_jp2 * 6] - F_bv[pos_b_jm1 * 6])) / DY;
				bvxz = (b1 * (F_bv[pos_b_kp1 * 6] - F_bv[pos_b * 6]) + b2 * (F_bv[pos_b_kp2 * 6] - F_bv[pos_b_km1 * 6])) / DZ;
				bvyx = (b1 * (F_bv[pos_b_ip1 * 6 + 1] - F_bv[pos_b * 6 + 1]) + b2 * (F_bv[pos_b_ip2 * 6 + 1] - F_bv[pos_b_im1 * 6 + 1])) / DX;
				bvyy = (b1 * (F_bv[pos_b * 6 + 1] - F_bv[pos_b_jm1 * 6 + 1]) + b2 * (F_bv[pos_b_jp1 * 6 + 1] - F_bv[pos_b_jm2 * 6 + 1])) / DY;
				bvyz = (b1 * (F_bv[pos_b_kp1 * 6 + 1] - F_bv[pos_b * 6 + 1]) + b2 * (F_bv[pos_b_kp2 * 6 + 1] - F_bv[pos_b_km1 * 6 + 1])) / DZ;
				bvzx = (b1 * (F_bv[pos_b_ip1 * 6 + 2] - F_bv[pos_b * 6 + 2]) + b2 * (F_bv[pos_b_ip2 * 6 + 2] - F_bv[pos_b_im1 * 6 + 2])) / DX;
				bvzy = (b1 * (F_bv[pos_b_jp1 * 6 + 2] - F_bv[pos_b * 6 + 2]) + b2 * (F_bv[pos_b_jp2 * 6 + 2] - F_bv[pos_b_jm1 * 6 + 2])) / DY;
				bvzz = (b1 * (F_bv[pos_b * 6 + 2] - F_bv[pos_b_km1 * 6 + 2]) + b2 * (F_bv[pos_b_kp1 * 6 + 2] - F_bv[pos_b_km2 * 6 + 2])) / DY;

				bivxx = (b1 * (F_bv[pos_b * 6 + 3] - F_bv[pos_b_im1 * 6 + 3]) + b2 * (F_bv[pos_b_ip1 * 6 + 3] - F_bv[pos_b_im2 * 6 + 3])) / DX;
				bivxy = (b1 * (F_bv[pos_b_jp1 * 6 + 3] - F_bv[pos_b * 6 + 3]) + b2 * (F_bv[pos_b_jp2 * 6 + 3] - F_bv[pos_b_jm1 * 6 + 3])) / DY;
				bivxz = (b1 * (F_bv[pos_b_kp1 * 6 + 3] - F_bv[pos_b * 6 + 3]) + b2 * (F_bv[pos_b_kp2 * 6 + 3] - F_bv[pos_b_km1 * 6 + 3])) / DZ;
				bivyx = (b1 * (F_bv[pos_b_ip1 * 6 + 4] - F_bv[pos_b * 6 + 4]) + b2 * (F_bv[pos_b_ip2 * 6 + 4] - F_bv[pos_b_im1 * 6 + 4])) / DX;
				bivyy = (b1 * (F_bv[pos_b * 6 + 4] - F_bv[pos_b_jm1 * 6 + 4]) + b2 * (F_bv[pos_b_jp1 * 6 + 4] - F_bv[pos_b_jm2 * 6 + 4])) / DY;
				bivyz = (b1 * (F_bv[pos_b_kp1 * 6 + 4] - F_bv[pos_b * 6 + 4]) + b2 * (F_bv[pos_b_kp2 * 6 + 4] - F_bv[pos_b_km1 * 6 + 4])) / DZ;
				bivzx = (b1 * (F_bv[pos_b_ip1 * 6 + 5] - F_bv[pos_b * 6 + 5]) + b2 * (F_bv[pos_b_ip2 * 6 + 5] - F_bv[pos_b_im1 * 6 + 5])) / DX;
				bivzy = (b1 * (F_bv[pos_b_jp1 * 6 + 5] - F_bv[pos_b * 6 + 5]) + b2 * (F_bv[pos_b_jp2 * 6 + 5] - F_bv[pos_b_jm1 * 6 + 5])) / DY;
				bivzz = (b1 * (F_bv[pos_b * 6 + 5] - F_bv[pos_b_km1 * 6 + 5]) + b2 * (F_bv[pos_b_kp1 * 6 + 5] - F_bv[pos_b_km2 * 6 + 5])) / DY;


			

				relam=0.0; imlam=0.0;   /*relam and imlam correponds to partial derivative wavefields with respect to lambda*/ /*Geändert wg. f als Geschwindigkeit, b als displacement*/
				relam=(fivxx+fivyy+fivzz)*(bvxx+bvyy+bvzz)+(fvxx+fvyy+fvzz)*(bivxx+bivyy+bivzz);
				relam=fdummy*relam;
				imlam=(fivxx+fivyy+fivzz)*(bivxx+bivyy+bivzz)-(fvxx+fvyy+fvzz)*(bvxx+bvyy+bvzz);
				imlam=fdummy*imlam;
				
				remu=0.0; immu=0.0;		
				remu=2*fivxx*bvxx+2*fivyy*bvyy+2*fivzz*bvzz+2*fvxx*bivxx+2*fvyy*bivyy+2*fvzz*bivzz+(fivxy+fivyx)*(bvxy+bvyx)+(fvxy+fvyx)*(bivxy+bivyx)+(fivxz+fivzx)*(bvxz+bvzx)+(fvxz+fvzx)*(bivxz+bivzx)+(fivyz+fivzy)*(bvyz+bvzy)+(fvyz+fvzy)*(bivyz+bivzy);
				remu=fdummy*remu;	
				immu=2*fivxx*bivxx+2*fivyy*bivyy+2*fivzz*bivzz-2*fvxx*bvxx-2*fvyy*bvyy-2*fvzz*bvzz+(fivxy+fivyx)*(bivxy+bivyx)-(fvxy+fvyx)*(bvxy+bvyx)+(fivxz+fivzx)*(bivxz+bivzx)-(fvxz+fvzx)*(bvxz+bvzx)+(fivyz+fivzy)*(bivyz+bivzy)-(fvyz+fvzy)*(bvyz+bvzy);
				immu=fdummy*immu;
				
				rerho=0.0; imrho=0.0;/*rho noch nicht abgeändert!!*/
				rerho = (F_fv[pos_f * 6] * F_bv[pos_b * 6] + F_fv[pos_f * 6 + 1] * F_bv[pos_b * 6 + 1] + F_fv[pos_f * 6 + 2] * F_bv[pos_b * 6 + 2]) - (F_fv[pos_f * 6 + 3] * F_bv[pos_b * 6 + 3] + F_fv[pos_f * 6 + 4] * F_bv[pos_b * 6 + 4] + F_fv[pos_f * 6 + 5] * F_bv[pos_b * 6 + 5]);
				imrho = (F_fv[pos_f * 6] * F_bv[pos_b * 6 + 3] + F_fv[pos_f * 6 + 1] * F_bv[pos_b * 6 + 4] + F_fv[pos_f * 6 + 2] * F_bv[pos_b * 6 + 5]) + (F_fv[pos_f * 6 + 3] * F_bv[pos_b * 6] + F_fv[pos_f * 6 + 4] * F_bv[pos_b * 6 + 1] + F_fv[pos_f * 6 + 5] * F_bv[pos_b * 6 + 2]);
				
				revs=0.0; imvs=0.0;
				revs=-2*relam+remu;
				imvs=-2*imlam+immu;
				
				rerho1=0.0; imrho1=0.0;


				rerho1=rerho+relam*(Fpi[pos_pi]-2*Fu[pos_pi])/Frho[pos_pi]+remu*Fu[pos_pi]/Frho[pos_pi];
				imrho1=imrho+imlam*(Fpi[pos_pi]-2*Fu[pos_pi])/Frho[pos_pi]+immu*Fu[pos_pi]/Frho[pos_pi];

				
				
				hess1[j][i][k]+=(relam*relam+imlam*imlam)*Frho[pos_pi]*Fpi[pos_pi]*4;
				hess2[j][i][k]+=(revs*revs+imvs*imvs)*4*Frho[pos_pi]*Fu[pos_pi];
				hess3[j][i][k]+=rerho1*rerho1+imrho1*imrho1;
				
							
				

	
				}
			}
		}
	}
	
    free_trans(Fpi, 1, NY + 1, 1, NX + 1, 1, NZ + 1);
    free_trans(Fu, 1, NY + 1, 1, NX + 1, 1, NZ + 1);
    free_trans(Frho, 1, NY + 1, 1, NX + 1, 1, NZ + 1);
	
}
