/*
 * gradient_F_slave.c
 *
 *  Created on: Dec 19, 2017
 *      Author: bingo
 */



#include <slave.h>
// #include <dma.h>
#include <math.h>

#define wx 10
#define wy 1
#define MX 8
#define MZ 8
//#define DEBUG
#define RMA
//#define RMA0
#define RMAP

typedef struct {
	int dim_x;
	int dim_y;
	int dim_z;
	int cube;
	int slice;
	int strip;
	int slice_pi;
	int strip_pi;
	int slice_grad;
	int strip_grad;
	float DT;
	float DX;
	float DY;
	float DZ;
	int FDCOEFF;
	int nf;
	float *Fgrad;
	float *F_fv;
	float *F_bv;
	float *Fpi;
	float *Fu;
	float *Frho;
	float *finv;
	int MYID;
} Param_grad;


inline  addn(int x, int n, int x1) {
	return x + n < x1 ? x + n : x + n - x1;
}

volatile __thread_local int reply,reply_l,reply_r,l_reply,r_reply,reply_mask;

void readData_dma_with_rma4(float *p,float *data,float *data1,int wy_tmp,int wz,int y0,int y1,int z0,int z1,int ystep_w,
                                int len76,int len75,int len74,int len32,int len31,int len30){     
    int dz = z0+z1;
    int dy = y0+y1;
    float *data_ldm = ldm_malloc(((wy_tmp+dy)/4+1)*(wz*4+dz)*sizeof(float));
    reply=0;
    int row_4aver_remainder=0;
    int flag=0;
    int rem_row=0;
    float *r_saddr = data;//7号从核的首地址-包括Halo
    float *r_saddr1 = data1;//3号从核的首地址-包括Halo
    int len7 = len74 + wz + dz; 
    int len3 = len30 + wz + dz; 
    float *l_saddr = data_ldm;
    l_reply=0;r_reply=0;
    athread_ssync_array();
    if(_ROW==0){
        row_4aver_remainder = (wy_tmp+dy) % 4; //need均分(wy+d)%4=0/1/2/3
        flag = ((row_4aver_remainder-(_COL%4))>=1)? 1:0;
        rem_row = (flag==1)? (_COL%4):row_4aver_remainder;
        if(_COL/4==0){//0-3
            r_saddr1 = r_saddr1 + (((wy_tmp+dy)/4)*(_COL%4)+rem_row)*ystep_w;
            for(int i=0;i<((wy_tmp+dy)/4+flag);i++){
                athread_dma_iget(l_saddr,r_saddr1, sizeof(float)*len3, &reply);
                l_saddr+=len3;
                r_saddr1+=ystep_w;
                athread_dma_wait_value(&reply,i+1); 
                //rma
                athread_rma_iput(data_ldm+i*len3,&l_reply,sizeof(float)*(wz+dz),_PEN+(3-0-_COL%4),p+((_COL%4)*((wy_tmp+dy)/4)+rem_row)*(wz+dz)+i*(wz+dz),&r_reply); 
                athread_rma_iput(data_ldm+i*len3+len32,&l_reply,sizeof(float)*(wz+dz),_PEN+(3-1-_COL%4),p+((_COL%4)*((wy_tmp+dy)/4)+rem_row)*(wz+dz)+i*(wz+dz),&r_reply);
                athread_rma_iput(data_ldm+i*len3+len31,&l_reply,sizeof(float)*(wz+dz),_PEN+(3-2-_COL%4),p+((_COL%4)*((wy_tmp+dy)/4)+rem_row)*(wz+dz)+i*(wz+dz),&r_reply);
                athread_rma_iput(data_ldm+i*len3+len30,&l_reply,sizeof(float)*(wz+dz),_PEN+(3-3-_COL%4),p+((_COL%4)*((wy_tmp+dy)/4)+rem_row)*(wz+dz)+i*(wz+dz),&r_reply);
                athread_rma_wait_value(&l_reply,(i+1)*4);
            }
        }else{//4-7
            r_saddr = r_saddr + (((wy_tmp+dy)/4)*(_COL%4)+rem_row)*ystep_w;
            for(int i=0;i<((wy_tmp+dy)/4+flag);i++){
                athread_dma_iget(l_saddr,r_saddr, sizeof(float)*len7, &reply);
                l_saddr+=len7;
                r_saddr+=ystep_w;
                athread_dma_wait_value(&reply,i+1);     
                //rma
                athread_rma_iput(data_ldm+i*len7,&l_reply,sizeof(float)*(wz+dz),_PEN+(3-0-_COL%4),p+((_COL%4)*((wy_tmp+dy)/4)+rem_row)*(wz+dz)+i*(wz+dz),&r_reply);
                athread_rma_iput(data_ldm+i*len7+len76,&l_reply,sizeof(float)*(wz+dz),_PEN+(3-1-_COL%4),p+((_COL%4)*((wy_tmp+dy)/4)+rem_row)*(wz+dz)+i*(wz+dz),&r_reply);
                athread_rma_iput(data_ldm+i*len7+len75,&l_reply,sizeof(float)*(wz+dz),_PEN+(3-2-_COL%4),p+((_COL%4)*((wy_tmp+dy)/4)+rem_row)*(wz+dz)+i*(wz+dz),&r_reply);
                athread_rma_iput(data_ldm+i*len7+len74,&l_reply,sizeof(float)*(wz+dz),_PEN+(3-3-_COL%4),p+((_COL%4)*((wy_tmp+dy)/4)+rem_row)*(wz+dz)+i*(wz+dz),&r_reply);
                athread_rma_wait_value(&l_reply,(i+1)*4);
            }
        }
    }else{//_ROW=1234567
        row_4aver_remainder = (wy_tmp) % 4; //need均分(wy+d)%4=0/1/2/3
        flag = ((row_4aver_remainder-(_COL%4))>=1)? 1:0;
        rem_row = (flag==1)? (_COL%4):row_4aver_remainder;
        if(_COL/4==0){
            r_saddr1 = r_saddr1 + ((wy_tmp+dy)+(wy_tmp*(_ROW-1))+(wy_tmp/4)*(_COL%4)+rem_row)*ystep_w;
            for(int i=0;i<(wy_tmp/4+flag);i++){
                athread_dma_iget(l_saddr,r_saddr1, sizeof(float)*len3, &reply);
                l_saddr+=len3;
                r_saddr1+=ystep_w;
                athread_dma_wait_value(&reply,i+1); 
                //rma
                athread_rma_iput(data_ldm+i*len3,&l_reply,sizeof(float)*(wz+dz),_PEN+(3-0-_COL%4),p+(dy+(_COL%4)*(wy_tmp/4)+rem_row)*(wz+dz)+i*(wz+dz),&r_reply); 
                athread_rma_iput(data_ldm+i*len3+len32,&l_reply,sizeof(float)*(wz+dz),_PEN+(3-1-_COL%4),p+(dy+(_COL%4)*(wy_tmp/4)+rem_row)*(wz+dz)+i*(wz+dz),&r_reply); 
                athread_rma_iput(data_ldm+i*len3+len31,&l_reply,sizeof(float)*(wz+dz),_PEN+(3-2-_COL%4),p+(dy+(_COL%4)*(wy_tmp/4)+rem_row)*(wz+dz)+i*(wz+dz),&r_reply); 
                athread_rma_iput(data_ldm+i*len3+len30,&l_reply,sizeof(float)*(wz+dz),_PEN+(3-3-_COL%4),p+(dy+(_COL%4)*(wy_tmp/4)+rem_row)*(wz+dz)+i*(wz+dz),&r_reply); 
                athread_rma_wait_value(&l_reply,(i+1)*4);
            }
        }else{
            r_saddr = r_saddr + ((wy_tmp+dy)+(wy_tmp*(_ROW-1))+(wy_tmp/4)*(_COL%4)+rem_row)*ystep_w;
            for(int i=0;i<(wy_tmp/4+flag);i++){
                athread_dma_iget(l_saddr,r_saddr, sizeof(float)*len7, &reply);
                l_saddr+=len7;
                r_saddr+=ystep_w;
                athread_dma_wait_value(&reply,i+1); 
                // rma
                athread_rma_iput(data_ldm+i*len7,&l_reply,sizeof(float)*(wz+dz),_PEN+(3-0-_COL%4),p+(dy+(_COL%4)*(wy_tmp/4)+rem_row)*(wz+dz)+i*(wz+dz),&r_reply); 
                athread_rma_iput(data_ldm+i*len7+len76,&l_reply,sizeof(float)*(wz+dz),_PEN+(3-1-_COL%4),p+(dy+(_COL%4)*(wy_tmp/4)+rem_row)*(wz+dz)+i*(wz+dz),&r_reply); 
                athread_rma_iput(data_ldm+i*len7+len75,&l_reply,sizeof(float)*(wz+dz),_PEN+(3-2-_COL%4),p+(dy+(_COL%4)*(wy_tmp/4)+rem_row)*(wz+dz)+i*(wz+dz),&r_reply); 
                athread_rma_iput(data_ldm+i*len7+len74,&l_reply,sizeof(float)*(wz+dz),_PEN+(3-3-_COL%4),p+(dy+(_COL%4)*(wy_tmp/4)+rem_row)*(wz+dz)+i*(wz+dz),&r_reply); 
                athread_rma_wait_value(&l_reply,(i+1)*4);
            }
        }
    }
    // -group组间 上下 rma通信
    //RMA: p->p
    l_reply=0;r_reply=0;
    athread_ssync_array();//同步上述DMA、RMA操作
    if(_ROW!=7){
        athread_rma_iput(p+wy_tmp*(wz+dz),&l_reply,sizeof(float)*(dy*(wz+dz)),_PEN+8,p,&r_reply); 
        athread_rma_wait_value(&l_reply,1);
    }
    if(_ROW!=0){
        athread_rma_wait_value(&r_reply,1);
    }
    ldm_free(data_ldm,((wy_tmp+dy)/4+1)*(wz*4+dz)*sizeof(float));
}
//dma rma try2
void readData_rma4_1(float *data_ldm_tmp, float *data,float *data1,int wy_tmp,int wz,int y0,int y1,int z0,int z1,int ystep_w,
                                    int len74, int len30){  
    int dz = z0+z1;
    int dy = y0+y1;
    reply_mask=0;
    int row_4aver_remainder=0;
    int flag=0;
    int rem_row=0;
    float *r_saddr = data;//7号从核的首地址-包括Halo
    float *r_saddr1 = data1;//3号从核的首地址-包括Halo
    int len7 = len74 + wz + dz; 
    int len3 = len30 + wz + dz; 
    float *l_saddr = data_ldm_tmp;
    reply_l=0;reply_r=0;
    athread_ssync_array();

    if(_ROW==0){
        row_4aver_remainder = (wy_tmp+dy) % 4; 
        flag = ((row_4aver_remainder-(_COL%4))>=1)? 1:0;
        rem_row = (flag==1)? (_COL%4):row_4aver_remainder;
        if(_COL/4==0){//0-3
            r_saddr1 = r_saddr1 + (((wy_tmp+dy)/4)*(_COL%4)+rem_row)*ystep_w;
            for(int i=0;i<((wy_tmp+dy)/4+flag);i++){
                athread_dma_iget(l_saddr,r_saddr1, sizeof(float)*len3, &reply_mask);
                l_saddr+=len3;
                r_saddr1+=ystep_w;
            }
        }else{//4-7
            r_saddr = r_saddr + (((wy_tmp+dy)/4)*(_COL%4)+rem_row)*ystep_w;
            for(int i=0;i<((wy_tmp+dy)/4+flag);i++){
                athread_dma_iget(l_saddr,r_saddr, sizeof(float)*len7, &reply_mask);
                l_saddr+=len7;
                r_saddr+=ystep_w;   
            }
        }
    }else{//_ROW=1234567
        row_4aver_remainder = (wy_tmp) % 4; 
        flag = ((row_4aver_remainder-(_COL%4))>=1)? 1:0;
        rem_row = (flag==1)? (_COL%4):row_4aver_remainder;
        if(_COL/4==0){
            r_saddr1 = r_saddr1 + ((wy_tmp+dy)+(wy_tmp*(_ROW-1))+(wy_tmp/4)*(_COL%4)+rem_row)*ystep_w;
            for(int i=0;i<(wy_tmp/4+flag);i++){
                athread_dma_iget(l_saddr,r_saddr1, sizeof(float)*len3, &reply_mask);
                l_saddr+=len3;
                r_saddr1+=ystep_w;
            }
        }else{
            r_saddr = r_saddr + ((wy_tmp+dy)+(wy_tmp*(_ROW-1))+(wy_tmp/4)*(_COL%4)+rem_row)*ystep_w;
            for(int i=0;i<(wy_tmp/4+flag);i++){
                athread_dma_iget(l_saddr,r_saddr, sizeof(float)*len7, &reply_mask);
                l_saddr+=len7;
                r_saddr+=ystep_w;
            }
        }
    }
}
    
void readData_rma4_2(float *data_ldm_tmp, float *wave_l,int wy_tmp,int wz,int y0,int y1,int z0,int z1,int ystep_w,
                                int len76,int len75,int len74,int len32,int len31,int len30){    
    int dz = z0+z1;
    int dy = y0+y1;
    int row_4aver_remainder=0;
    int flag=0;
    int rem_row=0;
    int len7 = len74 + wz + dz; 
    int len3 = len30 + wz + dz; 
    float *l_saddr = data_ldm_tmp;

    if(_ROW==0){
        row_4aver_remainder = (wy_tmp+dy) % 4; 
        flag = ((row_4aver_remainder-(_COL%4))>=1)? 1:0;
        rem_row = (flag==1)? (_COL%4):row_4aver_remainder;
        if(_COL/4==0){//0-3
            for(int i=0;i<((wy_tmp+dy)/4+flag);i++){
                athread_dma_wait_value(&reply_mask,i+1); 
                //rma
                athread_rma_iput(data_ldm_tmp+i*len3,&reply_l,sizeof(float)*(wz+dz),_PEN+(3-0-_COL%4),wave_l+((_COL%4)*((wy_tmp+dy)/4)+rem_row)*(wz+dz)+i*(wz+dz),&reply_r); 
                athread_rma_iput(data_ldm_tmp+i*len3+len32,&reply_l,sizeof(float)*(wz+dz),_PEN+(3-1-_COL%4),wave_l+((_COL%4)*((wy_tmp+dy)/4)+rem_row)*(wz+dz)+i*(wz+dz),&reply_r);
                athread_rma_iput(data_ldm_tmp+i*len3+len31,&reply_l,sizeof(float)*(wz+dz),_PEN+(3-2-_COL%4),wave_l+((_COL%4)*((wy_tmp+dy)/4)+rem_row)*(wz+dz)+i*(wz+dz),&reply_r);
                athread_rma_iput(data_ldm_tmp+i*len3+len30,&reply_l,sizeof(float)*(wz+dz),_PEN+(3-3-_COL%4),wave_l+((_COL%4)*((wy_tmp+dy)/4)+rem_row)*(wz+dz)+i*(wz+dz),&reply_r);
                athread_rma_wait_value(&reply_l,(i+1)*4);
            }
        }else{//4-7
            for(int i=0;i<((wy_tmp+dy)/4+flag);i++){
                athread_dma_wait_value(&reply_mask,i+1);     
                //rma
                athread_rma_iput(data_ldm_tmp+i*len7,&reply_l,sizeof(float)*(wz+dz),_PEN+(3-0-_COL%4),wave_l+((_COL%4)*((wy_tmp+dy)/4)+rem_row)*(wz+dz)+i*(wz+dz),&reply_r);
                athread_rma_iput(data_ldm_tmp+i*len7+len76,&reply_l,sizeof(float)*(wz+dz),_PEN+(3-1-_COL%4),wave_l+((_COL%4)*((wy_tmp+dy)/4)+rem_row)*(wz+dz)+i*(wz+dz),&reply_r);
                athread_rma_iput(data_ldm_tmp+i*len7+len75,&reply_l,sizeof(float)*(wz+dz),_PEN+(3-2-_COL%4),wave_l+((_COL%4)*((wy_tmp+dy)/4)+rem_row)*(wz+dz)+i*(wz+dz),&reply_r);
                athread_rma_iput(data_ldm_tmp+i*len7+len74,&reply_l,sizeof(float)*(wz+dz),_PEN+(3-3-_COL%4),wave_l+((_COL%4)*((wy_tmp+dy)/4)+rem_row)*(wz+dz)+i*(wz+dz),&reply_r);
                athread_rma_wait_value(&reply_l,(i+1)*4);
            }
        }
    }else{//_ROW=1234567
        row_4aver_remainder = (wy_tmp) % 4; 
        flag = ((row_4aver_remainder-(_COL%4))>=1)? 1:0;
        rem_row = (flag==1)? (_COL%4):row_4aver_remainder;
        if(_COL/4==0){
            for(int i=0;i<(wy_tmp/4+flag);i++){
                athread_dma_wait_value(&reply_mask,i+1); 
                //rma
                athread_rma_iput(data_ldm_tmp+i*len3,&reply_l,sizeof(float)*(wz+dz),_PEN+(3-0-_COL%4),wave_l+(dy+(_COL%4)*(wy_tmp/4)+rem_row)*(wz+dz)+i*(wz+dz),&reply_r); 
                athread_rma_iput(data_ldm_tmp+i*len3+len32,&reply_l,sizeof(float)*(wz+dz),_PEN+(3-1-_COL%4),wave_l+(dy+(_COL%4)*(wy_tmp/4)+rem_row)*(wz+dz)+i*(wz+dz),&reply_r); 
                athread_rma_iput(data_ldm_tmp+i*len3+len31,&reply_l,sizeof(float)*(wz+dz),_PEN+(3-2-_COL%4),wave_l+(dy+(_COL%4)*(wy_tmp/4)+rem_row)*(wz+dz)+i*(wz+dz),&reply_r); 
                athread_rma_iput(data_ldm_tmp+i*len3+len30,&reply_l,sizeof(float)*(wz+dz),_PEN+(3-3-_COL%4),wave_l+(dy+(_COL%4)*(wy_tmp/4)+rem_row)*(wz+dz)+i*(wz+dz),&reply_r); 
                athread_rma_wait_value(&reply_l,(i+1)*4);
            }
        }else{
            for(int i=0;i<(wy_tmp/4+flag);i++){
                athread_dma_wait_value(&reply_mask,i+1); 
                // rma
                athread_rma_iput(data_ldm_tmp+i*len7,&reply_l,sizeof(float)*(wz+dz),_PEN+(3-0-_COL%4),wave_l+(dy+(_COL%4)*(wy_tmp/4)+rem_row)*(wz+dz)+i*(wz+dz),&reply_r); 
                athread_rma_iput(data_ldm_tmp+i*len7+len76,&reply_l,sizeof(float)*(wz+dz),_PEN+(3-1-_COL%4),wave_l+(dy+(_COL%4)*(wy_tmp/4)+rem_row)*(wz+dz)+i*(wz+dz),&reply_r); 
                athread_rma_iput(data_ldm_tmp+i*len7+len75,&reply_l,sizeof(float)*(wz+dz),_PEN+(3-2-_COL%4),wave_l+(dy+(_COL%4)*(wy_tmp/4)+rem_row)*(wz+dz)+i*(wz+dz),&reply_r); 
                athread_rma_iput(data_ldm_tmp+i*len7+len74,&reply_l,sizeof(float)*(wz+dz),_PEN+(3-3-_COL%4),wave_l+(dy+(_COL%4)*(wy_tmp/4)+rem_row)*(wz+dz)+i*(wz+dz),&reply_r); 
                athread_rma_wait_value(&reply_l,(i+1)*4);
            }
        }
    }
    // -group组间 上下 rma通信
    reply_l=0;reply_r=0;
    athread_ssync_array();//同步上述DMA、RMA操作
    if(_ROW!=7){
        athread_rma_iput(wave_l+wy_tmp*(wz+dz),&reply_l,sizeof(float)*(dy*(wz+dz)),_PEN+8,wave_l,&reply_r); 
        athread_rma_wait_value(&reply_l,1);
    }
    if(_ROW!=0){
        athread_rma_wait_value(&reply_r,1);
    }
}
//dma rma try2
void gradient_kernel_slave(Param_grad *param) {

	int id = athread_get_id(-1);

	int dim_x = param->dim_x;
	int dim_y = param->dim_y;
	int dim_z = param->dim_z;
	int cube = param->cube;
	int slice = param->slice;
	int strip = param->strip;
	int slice_pi = param->slice_pi;
	int strip_pi = param->strip_pi;
	int slice_grad = param->slice_grad;
	int strip_grad = param->strip_grad;
	float DT = param->DT;
	float DX = param->DX;
	float DY = param->DY;
	float DZ = param->DZ;
	int FDCOEFF = param->FDCOEFF;
	int nf = param->nf;

	float *Fgrad = param->Fgrad;
	float *F_fv = param->F_fv;
	float *F_bv = param->F_bv;
	float *Fpi = param->Fpi;
	float *Fu = param->Fu;
	float *Frho = param->Frho;
	float *finv = param->finv;
	int MYID = param->MYID;

	int wz = 5;
	wz = (dim_z <= wz) ? dim_z : wz; 


	// int NX = (dim_x + wx * MX - 1) / (wx * MX);
	// int NZ = (dim_z + wz - 1) / wz;
	int yid = id / (MX * MZ);
    int xid = (id - yid * (MX * MZ)) / MZ;
    int zid = (id - yid * (MX * MZ)) % MZ;
	int nx = (dim_x + wx * MX  - 1) / (wx * MX);
	int nz = (dim_z + wz * MZ - 1) / (wz * MZ);



	const int x0 = 2;
	const int x1 = 2;
	const int y0 = 2;
	const int y1 = 2;
	const int z0 = 2;
	const int z1 = 2;

	int l, iix, iiz, ix, iy, iz;
	int izbeg, izend, izn, ixbeg, ixend, ixn;
	int izbeg0,izbeg1,izbeg2,izbeg3,izbeg4,izbeg5,izbeg6,izbeg7,ixbeg0;
	int len76,len75,len74,len32,len31,len30;
	int izend0,sum_izn;//put change


	// int xstep_6 = xstep * 6;
	// int ystep_6 = ystep * 6;
	// int xstep_3 = xstep * 3;
	// int ystep_3 = ystep * 3;
	// int ystep_grad_3 = ystep_grad * 3;
	// int xstep_grad_3 = xstep_grad * 3;
	// int xstep_l = wz + z0 + z1;
	// int ystep_l = (wx + x0 + x1) * (wz + z0 + z1);
	// int xstep_l2 = xstep_l * 2;
	// int ystep_l3 = ystep_l * 3;
	// int xstep_l3 = xstep_l * 3;
	// int ystep_l6 = ystep_l * 6;
	// int xstep_l6 = xstep_l * 6;


	// float *fvel_s_tmp, *fvel_s1, *bvel_s1, * bvel_s_tmp, *grad_s1, *grad_s1_tmp;
    // float *pi_s1, *pi_s1_tmp, *u_s1, *u_s1_tmp, *rho_s1, *rho_s1_tmp;
    float fvel[(wy + y0 + y1 + 1) * (wx + x0 + x1) * (wz + z0 + z1) * 6];
    float bvel[(wy + y0 + y1 + 1) * (wx + x0 + x1) * (wz + z0 + z1) * 6];
    float grad[wx * wz * 3];
    float u[wx * wz], pi[wx * wz], rho[wx * wz];
    // float *fvel_l_tmp, *fvel_l_c, *bvel_l_tmp, *bvel_l_c, *grad_l_tmp;
    // float *u_l_tmp, *pi_l_tmp, *rho_l_tmp;
    // fvel_l_c = fvel_l + x0 * xstep_l6 + z0 * 6;
    // bvel_l_c = bvel_l + x0 * xstep_l6 + z0 * 6;




	float finv_l[nf];
	float fvxx=0.0,fvxy=0.0,fvxz=0.0,fvyx=0.0,fvyy=0.0,fvyz=0.0,fvzx=0.0,fvzy=0.0,fvzz=0.0;
	float bvxx=0.0,bvxy=0.0,bvxz=0.0,bvyx=0.0,bvyy=0.0,bvyz=0.0,bvzx=0.0,bvzy=0.0,bvzz=0.0;
	float fivxx=0.0,fivxy=0.0,fivxz=0.0,fivyx=0.0,fivyy=0.0,fivyz=0.0,fivzx=0.0,fivzy=0.0,fivzz=0.0;
	float bivxx=0.0,bivxy=0.0,bivxz=0.0,bivyx=0.0,bivyy=0.0,bivyz=0.0,bivzx=0.0,bivzy=0.0,bivzz=0.0;
	float gradlam, gradmu, gradrho;
	float b1,b2,fdummy;
	int i,j,k,m,n,pos_f,pos_b,pos_pi,ll;
	int pos_f_jp1,pos_f_jp2,pos_f_jm1,pos_f_jm2,pos_b_jp1,pos_b_jp2,pos_b_jm1,pos_b_jm2;
	int pos_f_ip1,pos_f_ip2,pos_f_im1,pos_f_im2,pos_b_ip1,pos_b_ip2,pos_b_im1,pos_b_im2;
	int pos_f_kp1,pos_f_kp2,pos_f_km1,pos_f_km2,pos_b_kp1,pos_b_kp2,pos_b_km1,pos_b_km2;
	int slice_v,strip_v;
	int plane,plane_comp,plane_last,plane1,plane2,plane3,plane4,plane5;



	//To obtain data of finv
	reply = 0;
	athread_get(PE_MODE, finv, finv_l, nf * sizeof(float), &reply, 0, 0, 0);
	while(reply!=1);





	b1 = 9.0 / 8.0; b2 = -1.0 / 24.0; /* Taylor coefficients (4th order)*/
	if(FDCOEFF==2){
		b1=1.1382; b2=-0.046414;
	} /* Holberg coefficients E=0.1 %*/
    slice_v = (wx+x0+x1)*(wz+z0+z1);
	strip_v = (wz+z0+z1);


	for(l=1;l<=nf;l++){
		fdummy=0.0;
		fdummy=2.0*2.0*finv_l[l-1]*finv_l[l-1]*M_PI*M_PI;
		fdummy=1/fdummy;
		for (iiz = 0; iiz < nz; iiz++) {
			izbeg = dim_z - iiz*MZ*wz - wz*(zid + 1);
			izend = dim_z - iiz*MZ*wz - wz*zid;
			izbeg = izbeg < 0 ? 0 : izbeg;
			izn= izend - izbeg;
			izbeg7 = dim_z - iiz*MZ*wz - wz*(7 + 1);//第7列从核
			izbeg7 = izbeg7 < 0 ? 0 : izbeg7;
			izbeg6 = dim_z - iiz*MZ*wz - wz*(6 + 1);//第6列从核
			izbeg6 = izbeg6 < 0 ? 0 : izbeg6;
			izbeg5 = dim_z - iiz*MZ*wz - wz*(5 + 1);//第5列从核
			izbeg5 = izbeg5 < 0 ? 0 : izbeg5;
			izbeg4 = dim_z - iiz*MZ*wz - wz*(4 + 1);//第4列从核
			izbeg4 = izbeg4 < 0 ? 0 : izbeg4;
			izbeg3 = dim_z - iiz*MZ*wz - wz*(3 + 1);//第3列从核
			izbeg3 = izbeg3 < 0 ? 0 : izbeg3;
			izbeg2 = dim_z - iiz*MZ*wz - wz*(2 + 1);//第2列从核
			izbeg2 = izbeg2 < 0 ? 0 : izbeg2;
			izbeg1 = dim_z - iiz*MZ*wz - wz*(1 + 1);//第1列从核
			izbeg1 = izbeg1 < 0 ? 0 : izbeg1;
			izbeg0 = dim_z - iiz*MZ*wz - wz*(0 + 1);//第0列从核
			izbeg0 = izbeg0 < 0 ? 0 : izbeg0;
			len76 = izbeg6 - izbeg7;//第7、6列从核在z方向的差值
			len75 = izbeg5 - izbeg7; 
			len74 = izbeg4 - izbeg7;
			len32 = izbeg2 - izbeg3;//第3、2列从核在z方向的差值
			len31 = izbeg1 - izbeg3;
			len30 = izbeg0 - izbeg3;

			izend0 = dim_z - iiz*MZ*wz ;
			sum_izn = izend0 - izbeg7; //put change
			float grad_change[sum_izn*3];

			for (iix = 0; iix < nx; iix++) {
				ixbeg = wx * MX * iix + wx * xid;
				ixend = wx * MX * iix + wx * (xid + 1);
				ixend = ixend < dim_x ? ixend : dim_x;
				ixn = ixend - ixbeg;
				ixbeg0 = wx * MX * iix;
#ifdef RMA 
				for(ll = 0;ll<wy+y0+y1;ll++){ 
					// eight_coop(F_fv+(l*cube+(j+ll-y0)*slice+(ixbeg0+1-x0)*strip+izbeg7+1-z0)*6,F_fv+((j+ll-y0)*slice+(ixbeg0+1-x0)*strip+izbeg0+1-z0)*6,
					// fvel+ll*(wx+x0+x1)*(wz+z0+z1)*6,strip*6,wz*6,x0*6,x1*6,z0*6,z1*6);
					// eight_coop(F_bv+(l*cube+(j+ll-y0)*slice+(ixbeg0+1-x0)*strip+izbeg7+1-z0)*6,F_bv+((j+ll-y0)*slice+(ixbeg0+1-x0)*strip+izbeg0+1-z0)*6,
					// bvel+ll*(wx+x0+x1)*(wz+z0+z1)*6,strip*6,wz*6,x0,x1,z0*6,z1*6);
					readData_dma_with_rma4(fvel+ll*(wx+x0+x1)*(wz+z0+z1)*6,F_fv+(l*cube+(1+ll-y0)*slice+(ixbeg0+1-x0)*strip+izbeg7+1-z0)*6,
					F_fv+(l*cube+(1+ll-y0)*slice+(ixbeg0+1-x0)*strip+izbeg3+1-z0)*6, wx, wz*6, x0, x1, z0*6, z1*6, strip*6,
						len76*6, len75*6, len74*6, len32*6, len31*6, len30*6);
					readData_dma_with_rma4(bvel+ll*(wx+x0+x1)*(wz+z0+z1)*6,F_bv+(l*cube+(1+ll-y0)*slice+(ixbeg0+1-x0)*strip+izbeg7+1-z0)*6,
					F_bv+(l*cube+(1+ll-y0)*slice+(ixbeg0+1-x0)*strip+izbeg3+1-z0)*6, wx, wz*6, x0, x1, z0*6, z1*6, strip*6,
						len76*6, len75*6, len74*6, len32*6, len31*6, len30*6);
				}

#else
				if(izn > 0){
					for(ll = 0;ll<wy+y0+y1;ll++){
						for (i = ixbeg+1; i <= ixend+x0+x1; i++) {
							reply = 0;
							athread_get(PE_MODE, F_fv+(l*cube+(1+ll-y0)*slice+(i-x0)*strip+izbeg+1-z0)*6, fvel+ll*(wx+x0+x1)*(wz+z0+z1)*6+(i-ixbeg-1)*(wz+z0+z1)*6, (izn+z0+z1)*6*sizeof(float), &reply, 0, 0, 0);
							athread_get(PE_MODE, F_bv+(l*cube+(1+ll-y0)*slice+(i-x0)*strip+izbeg+1-z0)*6, bvel+ll*(wx+x0+x1)*(wz+z0+z1)*6+(i-ixbeg-1)*(wz+z0+z1)*6, (izn+z0+z1)*6*sizeof(float), &reply, 0, 0, 0);
							while(reply!=2);
						}
					}						
				}
#endif
                plane_comp = 0;
				plane_last = y0 + wy + y1;
				for (j = 1; j <= dim_y; j++) {

					reply_mask = 0;
					if(izn > 0 && j < dim_y && ixn > 0){
						for (i = ixbeg+1; i <= ixend+x0+x1; i++) {
							// reply_mask = 0;
							athread_get(PE_MODE, F_fv+(l*cube+(j+wy+y1)*slice+(i-x0)*strip+izbeg+1-z0)*6, fvel+plane_last*(wx+x0+x1)*(wz+z0+z1)*6+(i-ixbeg-1)*(wz+z0+z1)*6, (izn+z0+z1)*6*sizeof(float), &reply_mask, 0, 0, 0);
							athread_get(PE_MODE, F_bv+(l*cube+(j+wy+y1)*slice+(i-x0)*strip+izbeg+1-z0)*6, bvel+plane_last*(wx+x0+x1)*(wz+z0+z1)*6+(i-ixbeg-1)*(wz+z0+z1)*6, (izn+z0+z1)*6*sizeof(float), &reply_mask, 0, 0, 0);
							// while(reply_mask!=1);
						}
					}
					plane1 = plane_comp;
					plane2 = addn(plane_comp, 1, wy + y0 + y1 + 1);
					plane3 = addn(plane_comp, 2, wy + y0 + y1 + 1);
					plane4 = addn(plane_comp, 3, wy + y0 + y1 + 1);
					plane5 = addn(plane_comp, 4, wy + y0 + y1 + 1);
					plane =  plane3;
#ifdef RMA0
					readData_dma_with_rma4(grad,Fgrad+(j*slice_grad+(ixbeg0+1)*strip_grad+izbeg7+1)*3,
					Fgrad+(j*slice_grad+(ixbeg0+1)*strip_grad+izbeg3+1)*3, wx, wz*3, 0, 0, 0, 0, strip_grad*3,
						len76*3, len75*3, len74*3, len32*3, len31*3, len30*3);
					readData_dma_with_rma4(pi,Fpi+(j*slice_pi+(ixbeg0+1)*strip_pi+izbeg7+1),
					Fpi+(j*slice_pi+(ixbeg0+1)*strip_pi+izbeg3+1), wx, wz, 0, 0, 0, 0, strip_pi,
						len76, len75, len74, len32, len31, len30);
					readData_dma_with_rma4(u,Fu+(j*slice_pi+(ixbeg0+1)*strip_pi+izbeg7+1),
					Fu+(j*slice_pi+(ixbeg0+1)*strip_pi+izbeg3+1), wx, wz, 0, 0, 0, 0, strip_pi,
						len76, len75, len74, len32, len31, len30);
					readData_dma_with_rma4(rho,Frho+(j*slice_pi+(ixbeg0+1)*strip_pi+izbeg7+1),
					Frho+(j*slice_pi+(ixbeg0+1)*strip_pi+izbeg3+1), wx, wz, 0, 0, 0, 0, strip_pi,
						len76, len75, len74, len32, len31, len30);
#else
					if(izn > 0){
						for (i = ixbeg+1; i <= ixend; i++) {
							reply = 0;
							athread_get(PE_MODE, Fgrad+(j*slice_grad+i*strip_grad+izbeg+1)*3, grad+(i-ixbeg-1)*wz*3, izn*3*sizeof(float), &reply, 0, 0, 0);
							athread_get(PE_MODE, Fpi+(j*slice_pi+i*strip_pi+izbeg+1), pi+(i-ixbeg-1)*wz, izn*sizeof(float), &reply, 0, 0, 0);
							athread_get(PE_MODE, Fu+(j*slice_pi+i*strip_pi+izbeg+1), u+(i-ixbeg-1)*wz, izn*sizeof(float), &reply, 0, 0, 0);
							athread_get(PE_MODE, Frho+(j*slice_pi+i*strip_pi+izbeg+1), rho+(i-ixbeg-1)*wz, izn*sizeof(float), &reply, 0, 0, 0);
							while(reply!=4);
						}
				    }
#endif


					for (i = ixbeg+1; i <= ixend; i++) {
						for (k = izbeg+1; k <= izend; k++) {

					/* spatial derivatives of the components of the velocities are computed */
										
							// pos_f = l * cube +  j * slice + i * strip + k;
							// pos_b = l * cube +  j * slice + i * strip + k;
							// pos_pi = j * slice_pi + i * strip_pi + k;
							// pos_b_jm2 = pos_b - 2 * slice;
							// pos_b_jm1 = pos_b - slice;
							// pos_b_jp2 = pos_b + 2 * slice;
							// pos_b_jp1 = pos_b + slice;
							// pos_f_jm2 = pos_f - 2 * slice;
							// pos_f_jm1 = pos_f - slice;
							// pos_f_jp2 = pos_f + 2 * slice;
							// pos_f_jp1 = pos_f + slice;
							// pos_b_im2 = pos_b - 2 * strip;
							// pos_b_im1 = pos_b - strip;
							// pos_b_ip2 = pos_b + 2 * strip;
							// pos_b_ip1 = pos_b + strip;
							// pos_f_im2 = pos_f - 2 * strip;
							// pos_f_im1 = pos_f - strip;
							// pos_f_ip2 = pos_f + 2 * strip;
							// pos_f_ip1 = pos_f + strip;
							// pos_b_km2 = pos_b - 2;
							// pos_b_km1 = pos_b - 1;
							// pos_b_kp2 = pos_b + 2;
							// pos_b_kp1 = pos_b + 1;
							// pos_f_km2 = pos_f - 2;
							// pos_f_km1 = pos_f - 1;
							// pos_f_kp2 = pos_f + 2;
							// pos_f_kp1 = pos_f + 1;

							// fvxx = (b1 * (F_fv[pos_f * 6] - F_fv[pos_f_im1 * 6]) + b2 * (F_fv[pos_f_ip1 * 6] - F_fv[pos_f_im2 * 6])) / DX;
							// fvxy = (b1 * (F_fv[pos_f_jp1 * 6] - F_fv[pos_f * 6]) + b2 * (F_fv[pos_f_jp2 * 6] - F_fv[pos_f_jm1 * 6])) / DY;
							// fvxz = (b1 * (F_fv[pos_f_kp1 * 6] - F_fv[pos_f * 6]) + b2 * (F_fv[pos_f_kp2 * 6] - F_fv[pos_f_km1 * 6])) / DZ;
							// fvyx = (b1 * (F_fv[pos_f_ip1 * 6 + 1] - F_fv[pos_f * 6 + 1]) + b2 * (F_fv[pos_f_ip2 * 6 + 1] - F_fv[pos_f_im1 * 6 + 1])) / DX;
							// fvyy = (b1 * (F_fv[pos_f * 6 + 1] - F_fv[pos_f_jm1 * 6 + 1]) + b2 * (F_fv[pos_f_jp1 * 6 + 1] - F_fv[pos_f_jm2 * 6 + 1])) / DY;
							// fvyz = (b1 * (F_fv[pos_f_kp1 * 6 + 1] - F_fv[pos_f * 6 + 1]) + b2 * (F_fv[pos_f_kp2 * 6 + 1] - F_fv[pos_f_km1 * 6 + 1])) / DZ;
							// fvzx = (b1 * (F_fv[pos_f_ip1 * 6 + 2] - F_fv[pos_f * 6 + 2]) + b2 * (F_fv[pos_f_ip2 * 6 + 2] - F_fv[pos_f_im1 * 6 + 2])) / DX;
							// fvzy = (b1 * (F_fv[pos_f_jp1 * 6 + 2] - F_fv[pos_f * 6 + 2]) + b2 * (F_fv[pos_f_jp2 * 6 + 2] - F_fv[pos_f_jm1 * 6 + 2])) / DY;
							// fvzz = (b1 * (F_fv[pos_f * 6 + 2] - F_fv[pos_f_km1 * 6 + 2]) + b2 * (F_fv[pos_f_kp1 * 6 + 2] - F_fv[pos_f_km2 * 6 + 2])) / DY;

							// fivxx = (b1 * (F_fv[pos_f * 6 + 3] - F_fv[pos_f_im1 * 6 + 3]) + b2 * (F_fv[pos_f_ip1 * 6 + 3] - F_fv[pos_f_im2 * 6 + 3])) / DX;
							// fivxy = (b1 * (F_fv[pos_f_jp1 * 6 + 3] - F_fv[pos_f * 6 + 3]) + b2 * (F_fv[pos_f_jp2 * 6 + 3] - F_fv[pos_f_jm1 * 6 + 3])) / DY;
							// fivxz = (b1 * (F_fv[pos_f_kp1 * 6 + 3] - F_fv[pos_f * 6 + 3]) + b2 * (F_fv[pos_f_kp2 * 6 + 3] - F_fv[pos_f_km1 * 6 + 3])) / DZ;
							// fivyx = (b1 * (F_fv[pos_f_ip1 * 6 + 4] - F_fv[pos_f * 6 + 4]) + b2 * (F_fv[pos_f_ip2 * 6 + 4] - F_fv[pos_f_im1 * 6 + 4])) / DX;
							// fivyy = (b1 * (F_fv[pos_f * 6 + 4] - F_fv[pos_f_jm1 * 6 + 4]) + b2 * (F_fv[pos_f_jp1 * 6 + 4] - F_fv[pos_f_jm2 * 6 + 4])) / DY;
							// fivyz = (b1 * (F_fv[pos_f_kp1 * 6 + 4] - F_fv[pos_f * 6 + 4]) + b2 * (F_fv[pos_f_kp2 * 6 + 4] - F_fv[pos_f_km1 * 6 + 4])) / DZ;
							// fivzx = (b1 * (F_fv[pos_f_ip1 * 6 + 5] - F_fv[pos_f * 6 + 5]) + b2 * (F_fv[pos_f_ip2 * 6 + 5] - F_fv[pos_f_im1 * 6 + 5])) / DX;
							// fivzy = (b1 * (F_fv[pos_f_jp1 * 6 + 5] - F_fv[pos_f * 6 + 5]) + b2 * (F_fv[pos_f_jp2 * 6 + 5] - F_fv[pos_f_jm1 * 6 + 5])) / DY;
							// fivzz = (b1 * (F_fv[pos_f * 6 + 5] - F_fv[pos_f_km1 * 6 + 5]) + b2 * (F_fv[pos_f_kp1 * 6 + 5] - F_fv[pos_f_km2 * 6 + 5])) / DY;

							// bvxx = (b1 * (F_bv[pos_b * 6] - F_bv[pos_b_im1 * 6]) + b2 * (F_bv[pos_b_ip1 * 6] - F_bv[pos_b_im2 * 6])) / DX;
							// bvxy = (b1 * (F_bv[pos_b_jp1 * 6] - F_bv[pos_b * 6]) + b2 * (F_bv[pos_b_jp2 * 6] - F_bv[pos_b_jm1 * 6])) / DY;
							// bvxz = (b1 * (F_bv[pos_b_kp1 * 6] - F_bv[pos_b * 6]) + b2 * (F_bv[pos_b_kp2 * 6] - F_bv[pos_b_km1 * 6])) / DZ;
							// bvyx = (b1 * (F_bv[pos_b_ip1 * 6 + 1] - F_bv[pos_b * 6 + 1]) + b2 * (F_bv[pos_b_ip2 * 6 + 1] - F_bv[pos_b_im1 * 6 + 1])) / DX;
							// bvyy = (b1 * (F_bv[pos_b * 6 + 1] - F_bv[pos_b_jm1 * 6 + 1]) + b2 * (F_bv[pos_b_jp1 * 6 + 1] - F_bv[pos_b_jm2 * 6 + 1])) / DY;
							// bvyz = (b1 * (F_bv[pos_b_kp1 * 6 + 1] - F_bv[pos_b * 6 + 1]) + b2 * (F_bv[pos_b_kp2 * 6 + 1] - F_bv[pos_b_km1 * 6 + 1])) / DZ;
							// bvzx = (b1 * (F_bv[pos_b_ip1 * 6 + 2] - F_bv[pos_b * 6 + 2]) + b2 * (F_bv[pos_b_ip2 * 6 + 2] - F_bv[pos_b_im1 * 6 + 2])) / DX;
							// bvzy = (b1 * (F_bv[pos_b_jp1 * 6 + 2] - F_bv[pos_b * 6 + 2]) + b2 * (F_bv[pos_b_jp2 * 6 + 2] - F_bv[pos_b_jm1 * 6 + 2])) / DY;
							// bvzz = (b1 * (F_bv[pos_b * 6 + 2] - F_bv[pos_b_km1 * 6 + 2]) + b2 * (F_bv[pos_b_kp1 * 6 + 2] - F_bv[pos_b_km2 * 6 + 2])) / DY;

							// bivxx = (b1 * (F_bv[pos_b * 6 + 3] - F_bv[pos_b_im1 * 6 + 3]) + b2 * (F_bv[pos_b_ip1 * 6 + 3] - F_bv[pos_b_im2 * 6 + 3])) / DX;
							// bivxy = (b1 * (F_bv[pos_b_jp1 * 6 + 3] - F_bv[pos_b * 6 + 3]) + b2 * (F_bv[pos_b_jp2 * 6 + 3] - F_bv[pos_b_jm1 * 6 + 3])) / DY;
							// bivxz = (b1 * (F_bv[pos_b_kp1 * 6 + 3] - F_bv[pos_b * 6 + 3]) + b2 * (F_bv[pos_b_kp2 * 6 + 3] - F_bv[pos_b_km1 * 6 + 3])) / DZ;
							// bivyx = (b1 * (F_bv[pos_b_ip1 * 6 + 4] - F_bv[pos_b * 6 + 4]) + b2 * (F_bv[pos_b_ip2 * 6 + 4] - F_bv[pos_b_im1 * 6 + 4])) / DX;
							// bivyy = (b1 * (F_bv[pos_b * 6 + 4] - F_bv[pos_b_jm1 * 6 + 4]) + b2 * (F_bv[pos_b_jp1 * 6 + 4] - F_bv[pos_b_jm2 * 6 + 4])) / DY;
							// bivyz = (b1 * (F_bv[pos_b_kp1 * 6 + 4] - F_bv[pos_b * 6 + 4]) + b2 * (F_bv[pos_b_kp2 * 6 + 4] - F_bv[pos_b_km1 * 6 + 4])) / DZ;
							// bivzx = (b1 * (F_bv[pos_b_ip1 * 6 + 5] - F_bv[pos_b * 6 + 5]) + b2 * (F_bv[pos_b_ip2 * 6 + 5] - F_bv[pos_b_im1 * 6 + 5])) / DX;
							// bivzy = (b1 * (F_bv[pos_b_jp1 * 6 + 5] - F_bv[pos_b * 6 + 5]) + b2 * (F_bv[pos_b_jp2 * 6 + 5] - F_bv[pos_b_jm1 * 6 + 5])) / DY;
							// bivzz = (b1 * (F_bv[pos_b * 6 + 5] - F_bv[pos_b_km1 * 6 + 5]) + b2 * (F_bv[pos_b_kp1 * 6 + 5] - F_bv[pos_b_km2 * 6 + 5])) / DY;
						

							// gradlam=0.0;
							// gradlam=(fvxx+fvyy+fvzz)*(bvxx+bvyy+bvzz)+(fivxx+fivyy+fivzz)*(bivxx+bivyy+bivzz);
							// gradlam=-gradlam*fdummy;
							
							// gradmu=0.0;
							// gradmu= 2*fvxx*bvxx+2*fvyy*bvyy+2*fvzz*bvzz+2*fivxx*bivxx+2*fivyy*bivyy+2*fivzz*bivzz+(fvxy+fvyx)*(bvxy+bvyx)+(fivxy+fivyx)*(bivxy+bivyx)+(fvxz+fvzx)*(bvxz+bvzx)+(fivxz+fivzx)*(bivxz+bivzx)+(fvyz+fvzy)*(bvyz+bvzy)+(fivyz+fivzy)*(bivyz+bivzy);
							// gradmu=-gradmu*fdummy;
							
							// gradrho=0.0;
							// gradrho=(F_fv[pos_f * 6] * F_bv[pos_b * 6] + F_fv[pos_f * 6 + 1] * F_bv[pos_b * 6 + 1] + F_fv[pos_f * 6 + 2] * F_bv[pos_b * 6 + 2]) + (F_fv[pos_f * 6 + 3] * F_bv[pos_b * 6 + 3] + F_fv[pos_f * 6 + 4] * F_bv[pos_b * 6 + 4] + F_fv[pos_f * 6 + 5] * F_bv[pos_b * 6 + 5]);

							// //gradrho=gradrho;
							
							// /*parametrisation vp, vs, rho*/
							// Fgrad[3 * (j*slice_grad + i*strip_grad + k) + 0]+=sqrt(Frho[pos_pi]*Fpi[pos_pi])*2*gradlam; /*gradient vp*/
							
							// Fgrad[3 * (j*slice_grad + i*strip_grad + k) + 1]+=-4*sqrt(Frho[pos_pi]*Fu[pos_pi])*gradlam+2*sqrt(Frho[pos_pi]*Fu[pos_pi])*gradmu; /*gradient vs*/
					
							// Fgrad[3 * (j*slice_grad + i*strip_grad + k) + 2]+=gradrho+Fu[pos_pi]/Frho[pos_pi]*gradmu+(Fpi[pos_pi]-2*Fu[pos_pi])/Frho[pos_pi]*gradlam; /*gradient rho*/
					        
							pos_f = plane3 *slice_v + (i-ixbeg-1+x0)*strip_v + (k-izbeg-1+z0);
							pos_b = plane3 *slice_v + (i-ixbeg-1+x0)*strip_v + (k-izbeg-1+z0);
							pos_pi = (i-ixbeg-1) * wz + k-izbeg-1;
							pos_b_jm2 = plane1 *slice_v + (i-ixbeg-1+x0)*strip_v + (k-izbeg-1+z0);
							pos_b_jm1 = plane2 *slice_v + (i-ixbeg-1+x0)*strip_v + (k-izbeg-1+z0);
							pos_b_jp2 = plane5 *slice_v + (i-ixbeg-1+x0)*strip_v + (k-izbeg-1+z0);
							pos_b_jp1 = plane4 *slice_v + (i-ixbeg-1+x0)*strip_v + (k-izbeg-1+z0);
							pos_f_jm2 = plane1 *slice_v + (i-ixbeg-1+x0)*strip_v + (k-izbeg-1+z0);
							pos_f_jm1 = plane2 *slice_v + (i-ixbeg-1+x0)*strip_v + (k-izbeg-1+z0);
							pos_f_jp2 = plane5 *slice_v + (i-ixbeg-1+x0)*strip_v + (k-izbeg-1+z0);
							pos_f_jp1 = plane4 *slice_v + (i-ixbeg-1+x0)*strip_v + (k-izbeg-1+z0);
							pos_b_im2 = pos_b - 2 * strip_v;
							pos_b_im1 = pos_b - strip_v;
							pos_b_ip2 = pos_b + 2 * strip_v;
							pos_b_ip1 = pos_b + strip_v;
							pos_f_im2 = pos_f - 2 * strip_v;
							pos_f_im1 = pos_f - strip_v;
							pos_f_ip2 = pos_f + 2 * strip_v;
							pos_f_ip1 = pos_f + strip_v;
							pos_b_km2 = pos_b - 2;
							pos_b_km1 = pos_b - 1;
							pos_b_kp2 = pos_b + 2;
							pos_b_kp1 = pos_b + 1;
							pos_f_km2 = pos_f - 2;
							pos_f_km1 = pos_f - 1;
							pos_f_kp2 = pos_f + 2;
							pos_f_kp1 = pos_f + 1;
							// if(Fgrad[3*(j*slice_grad+i*strip_grad+k)+0] != grad[3 * pos_pi + 0]){
							// 	printf("grad 0 dma error \n");
							// }
							// if(Fgrad[3*(j*slice_grad+i*strip_grad+k)+1] != grad[3 * pos_pi + 1]){
							// 	printf("grad 1 dma error \n");
							// }
							// if(Fgrad[3*(j*slice_grad+i*strip_grad+k)+2] != grad[3 * pos_pi + 2]){
							// 	printf("grad 2 dma error \n");
							// }
							// if(F_bv[6*(l*cube+j*slice+i*strip+k)+1] != bvel[6*pos_b+1]){
							// 	printf("bv 0 dma error \n");
							// }
							// if(F_bv[6*(l*cube+(j+1)*slice+i*strip+k)+4] != bvel[6*pos_b_jp1+4]){
							// 	printf("bv 1 dma error \n");
							// }
							// if(F_bv[6*(l*cube+j*slice+(i-1)*strip+k)+2] != bvel[6*pos_b_im1+2]){
							// 	printf("bv 2 dma error \n");
							// }
							// if(F_fv[6*(l*cube+j*slice+i*strip+k)+1] != fvel[6*pos_f+1]){
							// 	printf("fv 0 dma error \n");
							// }
							// if(F_fv[6*(l*cube+(j+1)*slice+i*strip+k)+4] != fvel[6*pos_f_jp1+4]){
							// 	printf("fv 1 dma error \n");
							// }
							// if(F_fv[6*(l*cube+j*slice+(i-1)*strip+k)+2] != fvel[6*pos_f_im1+2]){
							// 	printf("fv 2 dma error \n");
							// }

							fvxx = (b1 * (fvel[pos_f * 6] - fvel[pos_f_im1 * 6]) + b2 * (fvel[pos_f_ip1 * 6] - fvel[pos_f_im2 * 6])) / DX;
							fvxy = (b1 * (fvel[pos_f_jp1 * 6] - fvel[pos_f * 6]) + b2 * (fvel[pos_f_jp2 * 6] - fvel[pos_f_jm1 * 6])) / DY;
							fvxz = (b1 * (fvel[pos_f_kp1 * 6] - fvel[pos_f * 6]) + b2 * (fvel[pos_f_kp2 * 6] - fvel[pos_f_km1 * 6])) / DZ;
							fvyx = (b1 * (fvel[pos_f_ip1 * 6 + 1] - fvel[pos_f * 6 + 1]) + b2 * (fvel[pos_f_ip2 * 6 + 1] - fvel[pos_f_im1 * 6 + 1])) / DX;
							fvyy = (b1 * (fvel[pos_f * 6 + 1] - fvel[pos_f_jm1 * 6 + 1]) + b2 * (fvel[pos_f_jp1 * 6 + 1] - fvel[pos_f_jm2 * 6 + 1])) / DY;
							fvyz = (b1 * (fvel[pos_f_kp1 * 6 + 1] - fvel[pos_f * 6 + 1]) + b2 * (fvel[pos_f_kp2 * 6 + 1] - fvel[pos_f_km1 * 6 + 1])) / DZ;
							fvzx = (b1 * (fvel[pos_f_ip1 * 6 + 2] - fvel[pos_f * 6 + 2]) + b2 * (fvel[pos_f_ip2 * 6 + 2] - fvel[pos_f_im1 * 6 + 2])) / DX;
							fvzy = (b1 * (fvel[pos_f_jp1 * 6 + 2] - fvel[pos_f * 6 + 2]) + b2 * (fvel[pos_f_jp2 * 6 + 2] - fvel[pos_f_jm1 * 6 + 2])) / DY;
							fvzz = (b1 * (fvel[pos_f * 6 + 2] - fvel[pos_f_km1 * 6 + 2]) + b2 * (fvel[pos_f_kp1 * 6 + 2] - fvel[pos_f_km2 * 6 + 2])) / DY;

							fivxx = (b1 * (fvel[pos_f * 6 + 3] - fvel[pos_f_im1 * 6 + 3]) + b2 * (fvel[pos_f_ip1 * 6 + 3] - fvel[pos_f_im2 * 6 + 3])) / DX;
							fivxy = (b1 * (fvel[pos_f_jp1 * 6 + 3] - fvel[pos_f * 6 + 3]) + b2 * (fvel[pos_f_jp2 * 6 + 3] - fvel[pos_f_jm1 * 6 + 3])) / DY;
							fivxz = (b1 * (fvel[pos_f_kp1 * 6 + 3] - fvel[pos_f * 6 + 3]) + b2 * (fvel[pos_f_kp2 * 6 + 3] - fvel[pos_f_km1 * 6 + 3])) / DZ;
							fivyx = (b1 * (fvel[pos_f_ip1 * 6 + 4] - fvel[pos_f * 6 + 4]) + b2 * (fvel[pos_f_ip2 * 6 + 4] - fvel[pos_f_im1 * 6 + 4])) / DX;
							fivyy = (b1 * (fvel[pos_f * 6 + 4] - fvel[pos_f_jm1 * 6 + 4]) + b2 * (fvel[pos_f_jp1 * 6 + 4] - fvel[pos_f_jm2 * 6 + 4])) / DY;
							fivyz = (b1 * (fvel[pos_f_kp1 * 6 + 4] - fvel[pos_f * 6 + 4]) + b2 * (fvel[pos_f_kp2 * 6 + 4] - fvel[pos_f_km1 * 6 + 4])) / DZ;
							fivzx = (b1 * (fvel[pos_f_ip1 * 6 + 5] - fvel[pos_f * 6 + 5]) + b2 * (fvel[pos_f_ip2 * 6 + 5] - fvel[pos_f_im1 * 6 + 5])) / DX;
							fivzy = (b1 * (fvel[pos_f_jp1 * 6 + 5] - fvel[pos_f * 6 + 5]) + b2 * (fvel[pos_f_jp2 * 6 + 5] - fvel[pos_f_jm1 * 6 + 5])) / DY;
							fivzz = (b1 * (fvel[pos_f * 6 + 5] - fvel[pos_f_km1 * 6 + 5]) + b2 * (fvel[pos_f_kp1 * 6 + 5] - fvel[pos_f_km2 * 6 + 5])) / DY;

							bvxx = (b1 * (bvel[pos_b * 6] - bvel[pos_b_im1 * 6]) + b2 * (bvel[pos_b_ip1 * 6] - bvel[pos_b_im2 * 6])) / DX;
							bvxy = (b1 * (bvel[pos_b_jp1 * 6] - bvel[pos_b * 6]) + b2 * (bvel[pos_b_jp2 * 6] - bvel[pos_b_jm1 * 6])) / DY;
							bvxz = (b1 * (bvel[pos_b_kp1 * 6] - bvel[pos_b * 6]) + b2 * (bvel[pos_b_kp2 * 6] - bvel[pos_b_km1 * 6])) / DZ;
							bvyx = (b1 * (bvel[pos_b_ip1 * 6 + 1] - bvel[pos_b * 6 + 1]) + b2 * (bvel[pos_b_ip2 * 6 + 1] - bvel[pos_b_im1 * 6 + 1])) / DX;
							bvyy = (b1 * (bvel[pos_b * 6 + 1] - bvel[pos_b_jm1 * 6 + 1]) + b2 * (bvel[pos_b_jp1 * 6 + 1] - bvel[pos_b_jm2 * 6 + 1])) / DY;
							bvyz = (b1 * (bvel[pos_b_kp1 * 6 + 1] - bvel[pos_b * 6 + 1]) + b2 * (bvel[pos_b_kp2 * 6 + 1] - bvel[pos_b_km1 * 6 + 1])) / DZ;
							bvzx = (b1 * (bvel[pos_b_ip1 * 6 + 2] - bvel[pos_b * 6 + 2]) + b2 * (bvel[pos_b_ip2 * 6 + 2] - bvel[pos_b_im1 * 6 + 2])) / DX;
							bvzy = (b1 * (bvel[pos_b_jp1 * 6 + 2] - bvel[pos_b * 6 + 2]) + b2 * (bvel[pos_b_jp2 * 6 + 2] - bvel[pos_b_jm1 * 6 + 2])) / DY;
							bvzz = (b1 * (bvel[pos_b * 6 + 2] - bvel[pos_b_km1 * 6 + 2]) + b2 * (bvel[pos_b_kp1 * 6 + 2] - bvel[pos_b_km2 * 6 + 2])) / DY;

							bivxx = (b1 * (bvel[pos_b * 6 + 3] - bvel[pos_b_im1 * 6 + 3]) + b2 * (bvel[pos_b_ip1 * 6 + 3] - bvel[pos_b_im2 * 6 + 3])) / DX;
							bivxy = (b1 * (bvel[pos_b_jp1 * 6 + 3] - bvel[pos_b * 6 + 3]) + b2 * (bvel[pos_b_jp2 * 6 + 3] - bvel[pos_b_jm1 * 6 + 3])) / DY;
							bivxz = (b1 * (bvel[pos_b_kp1 * 6 + 3] - bvel[pos_b * 6 + 3]) + b2 * (bvel[pos_b_kp2 * 6 + 3] - bvel[pos_b_km1 * 6 + 3])) / DZ;
							bivyx = (b1 * (bvel[pos_b_ip1 * 6 + 4] - bvel[pos_b * 6 + 4]) + b2 * (bvel[pos_b_ip2 * 6 + 4] - bvel[pos_b_im1 * 6 + 4])) / DX;
							bivyy = (b1 * (bvel[pos_b * 6 + 4] - bvel[pos_b_jm1 * 6 + 4]) + b2 * (bvel[pos_b_jp1 * 6 + 4] - bvel[pos_b_jm2 * 6 + 4])) / DY;
							bivyz = (b1 * (bvel[pos_b_kp1 * 6 + 4] - bvel[pos_b * 6 + 4]) + b2 * (bvel[pos_b_kp2 * 6 + 4] - bvel[pos_b_km1 * 6 + 4])) / DZ;
							bivzx = (b1 * (bvel[pos_b_ip1 * 6 + 5] - bvel[pos_b * 6 + 5]) + b2 * (bvel[pos_b_ip2 * 6 + 5] - bvel[pos_b_im1 * 6 + 5])) / DX;
							bivzy = (b1 * (bvel[pos_b_jp1 * 6 + 5] - bvel[pos_b * 6 + 5]) + b2 * (bvel[pos_b_jp2 * 6 + 5] - bvel[pos_b_jm1 * 6 + 5])) / DY;
							bivzz = (b1 * (bvel[pos_b * 6 + 5] - bvel[pos_b_km1 * 6 + 5]) + b2 * (bvel[pos_b_kp1 * 6 + 5] - bvel[pos_b_km2 * 6 + 5])) / DY;
						

							gradlam=0.0;
							gradlam=(fvxx+fvyy+fvzz)*(bvxx+bvyy+bvzz)+(fivxx+fivyy+fivzz)*(bivxx+bivyy+bivzz);
							gradlam=-gradlam*fdummy;
							
							gradmu=0.0;
							gradmu= 2*fvxx*bvxx+2*fvyy*bvyy+2*fvzz*bvzz+2*fivxx*bivxx+2*fivyy*bivyy+2*fivzz*bivzz+(fvxy+fvyx)*(bvxy+bvyx)+(fivxy+fivyx)*(bivxy+bivyx)+(fvxz+fvzx)*(bvxz+bvzx)+(fivxz+fivzx)*(bivxz+bivzx)+(fvyz+fvzy)*(bvyz+bvzy)+(fivyz+fivzy)*(bivyz+bivzy);
							gradmu=-gradmu*fdummy;
							
							gradrho=0.0;
							gradrho=(F_fv[pos_f * 6] * F_bv[pos_b * 6] + F_fv[pos_f * 6 + 1] * F_bv[pos_b * 6 + 1] + F_fv[pos_f * 6 + 2] * F_bv[pos_b * 6 + 2]) + (F_fv[pos_f * 6 + 3] * F_bv[pos_b * 6 + 3] + F_fv[pos_f * 6 + 4] * F_bv[pos_b * 6 + 4] + F_fv[pos_f * 6 + 5] * F_bv[pos_b * 6 + 5]);

							//gradrho=gradrho;
							
							/*parametrisation vp, vs, rho*/
							grad[3 * pos_pi + 0]+=sqrt(rho[pos_pi]*pi[pos_pi])*2*gradlam; /*gradient vp*/
							
							grad[3 * pos_pi + 1]+=-4*sqrt(rho[pos_pi]*u[pos_pi])*gradlam+2*sqrt(rho[pos_pi]*u[pos_pi])*gradmu; /*gradient vs*/
					
							grad[3 * pos_pi + 2]+=gradrho+u[pos_pi]/rho[pos_pi]*gradmu+(pi[pos_pi]-2*u[pos_pi])/rho[pos_pi]*gradlam; /*gradient rho*/
					        // if(Fpi[j*slice_pi+i*strip_pi+k] != pi[pos_pi]){
							// 	printf("pi dma error");
							// }
							// if(Fu[j*slice_pi+i*strip_pi+k] != u[pos_pi]){
							// 	printf("u dma error");
							// }
							// if(Frho[j*slice_pi+i*strip_pi+k] != rho[pos_pi]){
							// 	printf("rho dma error");
							// }
	
						}//k
#ifdef RMAP 
                    athread_ssync (ROW_SCOPE, 0xff);
					if(izn > 0){  
						l_reply = 0; 
						r_reply = 0;                
						athread_rma_iput(grad+(i-ixbeg-1)*wz*3,&l_reply,sizeof(float)*izn*3,_ROW*8+7,grad_change+(sum_izn-_COL*wz-izn)*3,&r_reply);
					    while(l_reply!=1);
					}
					athread_ssync (ROW_SCOPE, 0xff);
                    if(_COL == 7){
						reply = 0;
						athread_put(PE_MODE, grad_change, Fgrad+(j*slice_grad+i*strip_grad+izbeg+1)*3,sizeof(float)*sum_izn*3,&reply, 0, 0);
						while(reply!=1);
					}
#endif
					}//i
#ifdef RMAP
#else
					if(izn > 0){
						for (i = ixbeg+1; i <= ixend; i++) {
							reply = 0;
							athread_put(PE_MODE, grad+(i-ixbeg-1)*wz*3, Fgrad+(j*slice_grad+i*strip_grad+izbeg+1)*3, izn*3*sizeof(float),&reply, 0, 0);
							while(reply!=1);
						}
					}
#endif

                if(izn > 0 && j < dim_y && ixn > 0){
                    while(reply_mask!=2*(ixn+x1+x0));
                }

				plane_last = addn(plane_last, 1, wy + y0 + y1 + 1); 
                plane_comp = addn(plane_comp, 1, wy + y0 + y1 + 1); 
				}//j
			}//iix
		}//iiz
	}//nl



}















