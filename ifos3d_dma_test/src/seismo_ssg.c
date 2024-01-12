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
 *   store amplitudes (particle velocities or pressure) at receiver positions
     in arrays
 *  ----------------------------------------------------------------------*/

#include "fd.h"

void seismo(int lsamp, int ntr, int **recpos, float **sectionvx, float **sectionvy, 
float **sectionvz, float **sectiondiv, float **sectioncurl, float **sectionp, 
float *Fv, float *Fs, float ***pi, float ***u){
		
	extern int SEISMO; 
	int itr, ins, nxrec, nyrec, nzrec, pos, pos_pi, pos_jm1, pos_jp1, pos_im1, pos_ip1, pos_km1, pos_kp1;
	float amp, dh24x, dh24y, dh24z, vyx, vxy, vxx, vyy, vzx, vyz, vxz, vzy, vzz;
	extern float DX, DY, DZ;
	extern int NX, NY, NZ, MYID, FDORDER,  ABS_TYPE;
	int ll = 1;
	if(ABS_TYPE==1 && FDORDER==2) {ll=2; printf("ABS:%d, FDOR:%d\n", ABS_TYPE, FDORDER);}
    float *Fpi = transform3(pi, 1, NY + 1, 1, NX + 1, 1, NZ + 1);
    float *Fu = transform3(u, 1, NY + 1, 1, NX + 1, 1, NZ + 1);
	int ncol = NX + 2 * (ll*FDORDER/2), ndep = NZ + 2 * (ll*FDORDER/2);
	int strip = ndep;
	int slice = ncol * ndep;
	int strip_pi = NZ + 1;
	int slice_pi = (NX + 1) * (NZ + 1);
	//printf("debugging MYID:%d %s: ------------------------------------------------------------->\n", MYID,  __FUNCTION__);

	ins=lsamp; /* changed from "ins=lsamp/NDT;" (neccessary after correction of the buggy ns in ifos.c) */
	dh24x=1.0/DX;
	dh24y=1.0/DY;
	dh24z=1.0/DZ;
	
	for (itr=1;itr<=ntr;itr++){
		nxrec=recpos[1][itr];
		nyrec=recpos[2][itr];
		nzrec=recpos[3][itr];
		sectionvx[itr][ins]=0.0;
		sectionvy[itr][ins]=0.0;
		sectionvz[itr][ins]=0.0;
		pos = nyrec * slice + nxrec * strip + nzrec;
		pos_pi = nyrec * slice_pi + nxrec * strip_pi + nzrec;
		pos_jm1 = pos - slice;
		pos_jp1 = pos + slice;
		pos_im1 = pos - strip;
		pos_ip1 = pos + strip;
		pos_km1 = pos - 1;
		pos_kp1 = pos + 1;
		switch (SEISMO){
		case 1 :
			sectionvx[itr][ins]= Fv[pos * 3 + 0];
			sectionvy[itr][ins]= Fv[pos * 3 + 1];
			sectionvz[itr][ins]= Fv[pos * 3 + 2];
			 break;
		case 2 :
			sectionp[itr][ins] = - Fs[pos * 6 + 0] - Fs[pos * 6 + 1] - Fs[pos * 6 + 2];
			 break;
		case 3 :
			vxy=(Fv[pos_jp1 * 3 + 0] - Fv[pos * 3 + 0])*(dh24y);
		    vxz=(Fv[pos_kp1 * 3 + 0] - Fv[pos * 3 + 0])*(dh24z);
			vyx=(Fv[pos_ip1 * 3 + 1] - Fv[pos * 3 + 1])*(dh24x);
			vyz=(Fv[pos_kp1 * 3 + 1] - Fv[pos * 3 + 1])*(dh24z);
			vzx=(Fv[pos_ip1 * 3 + 2] - Fv[pos * 3 + 2])*(dh24x);
			vzy=(Fv[pos_jp1 * 3 + 2] - Fv[pos * 3 + 2])*(dh24y);
			amp=Fu[pos_pi]*((vyz-vzy)*fabs(vyz-vzy) + (vzx-vxz)*fabs(vzx-vxz)+(vxy-vyx)*fabs(vxy-vyx));
			sectioncurl[itr][ins]=fsign(amp)*sqrt(fabs(amp));

			/*vxx=(-vx[j][i+1][k]+27.0*(vx[j][i][k]-vx[j][i-1][k])+vx[j][i-2][k])*(24.0*DX);
			vyy=(-vy[j+1][i][k]+27.0*(vy[j][i][k]-vy[j-1][i][k])+vy[j-2][i][k])*(24.0*DY);
			vzz=(-vz[j][i][k+1]+27.0*(vz[j][i][k]-vz[j][i][k-1])+vz[j][i][k-2])*(24.0*DZ);*/
			
			vxx=(Fv[pos * 3 + 0] - Fv[pos_im1 * 3 + 0])*(dh24x);
			vyy=(Fv[pos * 3 + 1] - Fv[pos_jm1 * 3 + 1])*(dh24y);
			vzz=(Fv[pos * 3 + 2] - Fv[pos_km1 * 3 + 2])*(dh24z);
			
			sectiondiv[itr][ins]=(vxx+vyy+vzz)*sqrt(Fpi[pos_pi]);

			break;
		case 4 :				

			sectionvx[itr][ins] = Fv[pos * 3 + 0];
			sectionvy[itr][ins] = Fv[pos * 3 + 1];
			sectionvz[itr][ins] = Fv[pos * 3 + 2];
			sectionp[itr][ins] = - Fs[pos * 6 + 0] - Fs[pos * 6 + 1] - Fs[pos * 6 + 2];
			
			
			vxy=(Fv[pos_jp1 * 3 + 0] - Fv[pos * 3 + 0])*(dh24y);
		    vxz=(Fv[pos_kp1 * 3 + 0] - Fv[pos * 3 + 0])*(dh24z);
			vyx=(Fv[pos_ip1 * 3 + 1] - Fv[pos * 3 + 1])*(dh24x);
			vyz=(Fv[pos_kp1 * 3 + 1] - Fv[pos * 3 + 1])*(dh24z);
			vzx=(Fv[pos_ip1 * 3 + 2] - Fv[pos * 3 + 2])*(dh24x);
			vzy=(Fv[pos_jp1 * 3 + 2] - Fv[pos * 3 + 2])*(dh24y);
			
			amp=Fu[pos_pi]*((vyz-vzy)*fabs(vyz-vzy)+(vzx-vxz)*fabs(vzx-vxz)+(vxy-vyx)*fabs(vxy-vyx));
			sectioncurl[itr][ins]=fsign(amp)*sqrt(fabs(amp));

			
			vxx=(Fv[pos * 3 + 0] - Fv[pos_im1 * 3 + 0])*(dh24x);
			vyy=(Fv[pos * 3 + 1] - Fv[pos_jm1 * 3 + 1])*(dh24y);
			vzz=(Fv[pos * 3 + 2] - Fv[pos_km1 * 3 + 2])*(dh24z);
			
			sectiondiv[itr][ins]=(vxx+vyy+vzz)*sqrt(Fpi[pos_pi]);
			break;
			 		 
		}
	}
    free_trans(Fpi, 1, NY + 1, 1, NX + 1, 1, NZ + 1);
    free_trans(Fu, 1, NY + 1, 1, NX + 1, 1, NZ + 1);
}
