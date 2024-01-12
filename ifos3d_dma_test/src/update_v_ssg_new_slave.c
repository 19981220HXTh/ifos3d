#include <stdio.h>
#include <slave.h>
#include <dma.h>
#include "sw.h"
//#include "fd.h"

#define wy 1
#define wx 10
#define wz 4
//#define DEBUG

typedef struct Param_vel {
	int dim_x;
	int dim_y;
	int dim_z;
	int slice;
	int strip;
	int slice_rp;
	int strip_rp;
	float  dx;
	float  dy;
	float dz;
	int FDCOEFF;
	float *vel;
	float *stress;
	float * rp;
} Param_vel;

#define get_row_id(rid) asm volatile ("rcsr %0, 1" : "=r"(rid))
#define get_col_id(cid) asm volatile ("rcsr %0, 2" : "=r"(cid))

inline  addn(int x, int n, int x1) {
	return x + n < x1 ? x + n : x + n - x1;
}

inline int id_map(int id) {
	int i = id / 8;
	if (i % 2 == 1) {
		int beg = i * 8;
		int end = (i + 1) * 8 - 1;
		return beg + end - id;
	}
	return id;
}

void update_v_kernel_fusion_new_slave(Param_vel *param) 
{
    volatile int id , cid, rid;
    volatile int get_reply;
    volatile int put_reply;
	id = athread_get_id(-1);
    get_row_id(rid);
    get_col_id(cid);
	int dim_x = param->dim_x;
	int dim_y = param->dim_y;
	int dim_z = param->dim_z;
	int xstep = param->strip;
	int ystep = param->slice;
	int xstep_rp = param->strip_rp;
	int ystep_rp = param->slice_rp;
	float dx = param->dx;
	float dy = param->dy;
	float dz = param->dz;
	int FDCOEFF = param->FDCOEFF;

	float *vel_s0 = param->vel;
	float *strs_s0 = param->stress;
	float *rp_s0 = param->rp;

	const int x0 = 2;
	const int x1 = 2;
	const int y0 = 2;
	const int y1 = 2;
	const int z0 = 2;
	const int z1 = 2;
	const int dummy = 4;

	int i, j, k, ix, iy, iz, iix, iiz;
	int izbeg, izend, izn;
	int ixbeg, ixend, ixn;

	int xstep_6 = xstep*6;
	int ystep_6 = ystep*6;
	int xstep_3 = xstep*3;
	int ystep_3 = ystep*3;
	int ystep_rp3 = ystep_rp*3;
	int xstep_rp3 = xstep_rp*3;

	int xstep_l = wz + z0 + z1;
	int ystep_l = (wx + x0 + x1)*(wz + z0 + z1);
	int xstep_l2 = xstep_l*2;
	int ystep_l6 = ystep_l*6;
	int xstep_l6 = xstep_l*6;


    // nx, nz, respect the number dma which relations with x-axis and y-axis
	int nx = (dim_x + wx - 1)/wx;
	int nz = (dim_z + 8*wz - 1)/(8*wz);
    int loop_x = wx/8;

    float strs[(wy + y0 + y1)*(wx + x0 + x1 + dummy)*(wz + z0 + z1)*6];
	float vel[wx*wz*3];
	float rp[wx*wz*3];
    float buffer[(wz*8 + z0 + z1)*6];
	float *vel_s, *strs_s, *rp_s;
	float *strs_ls, *str_s1_temp;
    float *vel_ls, *rp_ls;
	float sxx_x, sxy_y, sxz_z, syy_y, sxy_x, syz_z, szz_z, sxz_x, syz_y;
	float *strs_c = strs + x0 * xstep_l6 + z0 * 6;

	int plane_comp = 0;
	float b1, b2;
	b1=9.0/8.0; b2=-1.0/24.0; /* Taylor coefficients*/
	if(FDCOEFF==2){
	b1=1.1382; b2=-0.046414;} /* Holberg coefficients E=0.1 %*/


#define DEBUG
#ifdef DEBUG
	int ldm_vel = sizeof(float)*wx*wz*3;
	int ldm_strs = sizeof(float)*(wx + x0 + x1 + dummy)*(wy + y0 + y1)*(wz + z0 + z1)*6;
	int ldm_rp = sizeof(float)*wx*wz*3;
    int ldm_buffer = sizeof(float)*(wz + z0 + z1)*8*6;
	int ldm_total = ldm_vel + ldm_strs + ldm_rp + ldm_buffer;
	if(id == 0) {
		//printf("LDM vel = %d Bytes\n", ldm_vel);
		//printf("LDM stress = %d Bytes\n", ldm_strs);
		//printf("LDM up = %d Bytes\n", ldm_rp);
		//printf("LDM total of update_v = %d Bytes\n", ldm_total);
	}
#endif

    for (iiz = 0; iiz < nz; iiz++) {
        for (iix = 0; iix < nx; iix++) {
            // obtain halo data for y-axis
            for (iy = -y0; iy < wy + y1; iy++) {
                // obtain data from mpe
                for (i = 0; i < loop_x + 1; i++) {
                    strs_s = strs_s0 + (iix*wx + i*8 + cid)*xstep_6 + rid*8*wz*6 + iy*ystep_6 - x0*xstep_6 - z0*6; 
                    strs_ls = strs + (y0 + iy)*ystep_l6 + i*8*xstep_l6; 
			        get_reply = 0;
                    // get wz*8 + z0 + z1 for buffer
			        athread_get(PE_MODE, strs_s0, buffer, (wz*8 + z0 + z1)*6*sizeof(float), &get_reply, 0, 0, 0);
			        while(get_reply!=1);
                    regshffle_strs(buffer, strs_ls, (int)((wz + z0 + z1)*6/4));
                }
            }
			plane_comp = 0;
            for (iy = 0; iy < dim_y; iy++) {
                for (i = 0; i < loop_x; i++) {
                    vel_s = vel_s0 + (iix*wx + i*8 + cid)*xstep_3 + rid*8*wz*3 + iy*ystep_3;
                    // why xstep_rp, becase rp has no halo
                    rp_s = rp_s0 + (iix*wx + i*8 + cid)*xstep_rp3 + rid*8*wz*3 + iy*ystep_rp3;
                    vel_ls = vel + i*8*wz*3;
                    rp_ls = rp + i*8*wz*3;
			        get_reply = 0;
			        athread_get(PE_MODE, vel_s, buffer, wz*8*3*sizeof(float), &get_reply, 0, 0, 0);
			        while(get_reply!=1);
                    regshffle_vel(buffer, vel_ls, (int)(wz*3/4));
			        get_reply = 0;
			        athread_get(PE_MODE, rp_s, buffer, wz*8*3*sizeof(float), &get_reply, 0, 0, 0);
			        while(get_reply!=1);
                    regshffle_vel(buffer, rp_ls, (int)(wz*3/4));
                }
		        athread_syn(ARRAY_SCOPE, 0xffff);\

				// computing part 
				for (ix = 0; ix < 1; ix++) {
					for (iz = 0; iz < 1; iz++) {
						int ym2 = plane_comp;
						int ym1 = addn(plane_comp, 1, wy + y0 + y1);
						int y = addn(plane_comp, y0 , wy + y0 + y1);
						int yp1 = addn(plane_comp, y0 + 1, wy + y0 + y1);
						int yp2 = addn(plane_comp, y0 + 2, wy + y0 + y1);

						int pos = y * ystep_l + ix * xstep_l + iz;

						int pos_km2_l = pos - 2;
						int pos_km1_l = pos - 1;
						int pos_kp1_l = pos + 1;
						int pos_kp2_l = pos + 2;

						int pos_im2_l = pos - xstep_l2;
						int pos_im1_l = pos - xstep_l;
						int pos_ip1_l = pos + xstep_l;
						int pos_ip2_l = pos + xstep_l2;

						int pos_jm2_l = ym2 * ystep_l + ix * xstep_l + iz;
						int pos_jm1_l = ym1 * ystep_l + ix * xstep_l + iz;
						int pos_jp1_l = yp1 * ystep_l + ix * xstep_l + iz;
						int pos_jp2_l = yp2 * ystep_l + ix * xstep_l + iz;

						int pos_v = ix * wz + iz;
						sxx_x = dx * (b1 * (strs_c[pos_ip1_l * 6 + 0] - strs_c[pos * 6 + 0]) + b2 * (strs_c[pos_ip2_l * 6 + 0] - strs_c[pos_im1_l * 6 + 0]));
						sxy_y = dy * (b1 * (strs_c[pos * 6 + 3] - strs_c[pos_jm1_l * 6 + 3]) + b2 * (strs_c[pos_jp1_l * 6 + 3] - strs_c[pos_jm2_l * 6 + 3]));
						sxz_z = dz * (b1 * (strs_c[pos * 6 + 5] - strs_c[pos_km1_l * 6 + 5]) + b2 * (strs_c[pos_kp1_l * 6 + 5] - strs_c[pos_km2_l * 6 + 5]));

						vel[pos_v * 3 + 0] += (sxx_x + sxy_y + sxz_z) / rp[pos_v * 3 + 0];

						syy_y = dy * (b1 * (strs_c[pos_jp1_l * 6 + 1] - strs_c[pos * 6 + 1]) + b2 * (strs_c[pos_jp2_l * 6 + 1] - strs_c[pos_jm1_l * 6 + 1]));
						sxy_x = dx * (b1 * (strs_c[pos * 6 + 3] - strs_c[pos_im1_l * 6 + 3]) + b2 * (strs_c[pos_ip1_l * 6 + 3] - strs_c[pos_im2_l * 6 + 3]));
						syz_z = dz * (b1 * (strs_c[pos * 6 + 4] - strs_c[pos_km1_l * 6 + 4]) + b2 * (strs_c[pos_kp1_l * 6 + 4] - strs_c[pos_km2_l * 6 + 4]));

						vel[pos_v * 3 + 1] += (syy_y + sxy_x + syz_z) / rp[pos_v * 3 + 1];

						szz_z = dz * (b1 * (strs_c[pos_kp1_l * 6 + 2] - strs_c[pos * 6 + 2]) + b2 * (strs_c[pos_kp2_l * 6 + 2] - strs_c[pos_km1_l * 6 + 2]));
						sxz_x = dx * (b1 * (strs_c[pos * 6 + 5] - strs_c[pos_im1_l * 6 + 5]) + b2 * (strs_c[pos_ip1_l * 6 + 5] - strs_c[pos_im2_l * 6 + 5]));
						syz_y = dy * (b1 * (strs_c[pos * 6 + 4] - strs_c[pos_jm1_l * 6 + 4]) + b2 * (strs_c[pos_jp1_l * 6 + 4] - strs_c[pos_jm2_l * 6 + 4]));

						vel[pos_v * 3 + 2] += (szz_z + sxz_x + syz_y) / rp[pos_v * 3 + 2];
                    }
                }
				
				
				for (i = 0; i < loop_x; i++) {
                    vel_s = vel_s0 + (iix*wx + i*8 + cid)*xstep_3 + rid*8*wz*3 + iy*ystep_3;
                    vel_ls = vel + i*8*wz*3;
					regshffle_vel(vel_ls, buffer , (int)(wz*3/4));
			        put_reply = 0;
			        athread_put(PE_MODE, buffer, vel_s, wz*8*3*sizeof(float), &put_reply, 0, 0);
			        while(put_reply!=1);
                }
				
				//get next plane stress data
                for (i = 0; i < loop_x + 1; i++) {
                    strs_s = strs_s0 + (iix*wx + i*8 + cid)*xstep_6 + rid*8*wz*6 + (iy + wy + y1)*ystep_6 - x0*xstep_6 - z0*6; 
                    strs_ls = strs + plane_comp*ystep_l6 + i*8*xstep_l6; 
                    // get wz*8 + z0 + z1 for buffer
			        get_reply = 0;
			        athread_get(PE_MODE, strs_s0, buffer, (wz*8 + z0 + z1)*6*sizeof(float), &get_reply, 0, 0, 0);
			        while(get_reply!=1);
                    regshffle_strs(buffer, strs_ls, (int)((wz + z0 + z1)*6/4));
                }
				plane_comp = addn(plane_comp, 1, wy + y0 + y1);
            }
        }
    }
}

