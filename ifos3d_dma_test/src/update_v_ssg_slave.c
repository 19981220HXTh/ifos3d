#include <stdio.h>
#include <slave.h>
#define wx 10
#define wy 1
#define MX 8
#define MZ 8
//#define DEBUG
#define pe_get(mem, ldm, size, reply) athread_get(PE_MODE, mem, ldm, size, (void*)(reply), 0, 0, 0) 
#define pe_put(mem, ldm, size, reply) athread_put(PE_MODE, ldm, mem, size, (void*)(reply), 0, 0)
#define RMA
//#define RMA0
#define RMAP

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
	// float *vel;
	// float *stress;
	// float * rp;
	float *Fv;
	float *Fs;
	float *Frp;
} Param_vel;


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
void update_v_kernel_fusion_slave(Param_vel *param) {

	int id = athread_get_id(-1);
	int dim_x = param->dim_x;
	int dim_y = param->dim_y;
	int dim_z = param->dim_z;
	int strip = param->strip;
	int slice = param->slice;
	int strip_rp = param->strip_rp;
	int slice_rp = param->slice_rp;
	float dx = param->dx;
	float dy = param->dy;
	float dz = param->dz;
	int FDCOEFF = param->FDCOEFF;

	float *Fv = param->Fv;
	float *Fs = param->Fs;
	float *Frp = param->Frp;

	int wz = 5;
	wz = (dim_z <= wz) ? dim_z : wz; 

	const int x0 = 2;
	const int x1 = 2;
	const int y0 = 2;
	const int y1 = 2;
	const int z0 = 2;
	const int z1 = 2;

	int ix, iy, iz, iix, iiz;
	int izbeg, izend, izn;
	int ixbeg, ixend, ixn;



	int yid = id / (MX * MZ);
    int xid = (id - yid * (MX * MZ)) / MZ;
    int zid = (id - yid * (MX * MZ)) % MZ;
	int nx = (dim_x + wx * MX  - 1) / (wx * MX);
	int nz = (dim_z + wz * MZ - 1) / (wz * MZ);

	float strs[(wy + y0 + y1 + 1) * (wx + x0 + x1) * (wz + z0 + z1) * 6];
	float vel[wx * wz * 3];
	float rp[wx * wz * 3];
	// float *vel_s1, * str_s1, *rp_s1;
	// float *str_temp, *str_s1_temp;



	float b1, b2;
	b1=9.0/8.0; b2=-1.0/24.0; /* Taylor coefficients*/
	if(FDCOEFF==2){
	b1=1.1382; b2=-0.046414;} /* Holberg coefficients E=0.1 %*/
    



	int i, j, k, l;
	float sxx_x, sxy_y, sxz_z, syy_y, sxy_x, syz_z;
	float szz_z, sxz_x, syz_y;

	int ll = 1;
	int idx,idx_rp,slice_str,strip_str;
	int izbeg0,izbeg1,izbeg2,izbeg3,izbeg4,izbeg5,izbeg6,izbeg7,ixbeg0;
	int len76,len75,len74,len32,len31,len30;
	int izend0,sum_izn;//put change
    int pos_jp1,pos_jp2,pos_jm1,pos_jm2,pos_ip1,pos_ip2,pos_im1,pos_im2,pos_kp1,pos_kp2,pos_km1,pos_km2;
	int plane,plane_comp,plane_last,plane1,plane2,plane3,plane4,plane5;

    //dma rma try2
    float data_ldm[ 1 * ((x0 + wx + x1)/4+1) * (z0 + wz*4 + z1)* 6];
    //dma rma try2



	/*unsigned int * test_float;*/



	/*if (LOG)
	if (MYID==0) time1=MPI_Wtime();*/
    slice_str = (wx+x0+x1)*(wz+z0+z1);
	strip_str = (wz+z0+z1);


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
		float vel_change[sum_izn*3];

		for (iix = 0; iix < nx; iix++) {
			ixbeg = wx * MX * iix + wx * xid;
			ixend = wx * MX * iix + wx * (xid + 1);
			ixend = ixend < dim_x ? ixend : dim_x;
			ixn = ixend - ixbeg;
            ixbeg0 = wx * MX * iix;
#ifdef RMA 
            for(l = 0;l<wy+y0+y1;l++){ 
                // eight_coop(Fs+((j+l-y0)*slice+(ixbeg0+1-x0)*strip+izbeg7+1-z0)*6,Fs+((j+l-y0)*slice+(ixbeg0+1-x0)*strip+izbeg0+1-z0)*6,
                // strs+l*(wx+x0+x1)*(wz+z0+z1)*6,strip*6,wz*6,x0,x1,z0*6,z1*6);
                readData_dma_with_rma4(strs+l*(wx+x0+x1)*(wz+z0+z1)*6,Fs+((1+l-y0)*slice+(ixbeg0+1-x0)*strip+izbeg7+1-z0)*6,
                Fs+((1+l-y0)*slice+(ixbeg0+1-x0)*strip+izbeg3+1-z0)*6, wx, wz*6, x0, x1, z0*6, z1*6, strip*6,
                    len76*6, len75*6, len74*6, len32*6, len31*6, len30*6);
            }

#else
            if(izn > 0){
                for(l = 0;l<wy+y0+y1;l++){
                    for (i = ixbeg+1; i <= ixend+x0+x1; i++) {
                        reply = 0;
                        athread_get(PE_MODE, Fs+((1+l-y0)*slice+(i-x0)*strip+izbeg+1-z0)*6, strs+l*(wx+x0+x1)*(wz+z0+z1)*6+(i-ixbeg-1)*(wz+z0+z1)*6, (izn+z0+z1)*6*sizeof(float), &reply, 0, 0, 0);
                        while(reply!=1);
                    }
                }
                
            }
#endif
            plane_comp = 0;
            plane_last = y0 + wy + y1;
			for (j = 1; j <= dim_y; j++) {
                //get next plane
#ifdef RMA
                readData_rma4_1(data_ldm,Fs+((j+wy+y1)*slice+(ixbeg0+1-x0)*strip+izbeg7+1-z0)*6,
                Fs+((j+wy+y1)*slice+(ixbeg0+1-x0)*strip+izbeg3+1-z0)*6, wx, wz*6, x0, x1, z0*6, z1*6, strip*6,
                    len74*6, len30*6);
#else
                reply_mask = 0;
                if(izn > 0 && j < dim_y && ixn > 0){
                    for (i = ixbeg+1; i <= ixend+x0+x1; i++) {
                        // reply_mask = 0;
                        athread_get(PE_MODE, Fs+((j+wy+y1)*slice+(i-x0)*strip+izbeg+1-z0)*6, strs+plane_last*(wx+x0+x1)*(wz+z0+z1)*6+(i-ixbeg-1)*(wz+z0+z1)*6, (izn+z0+z1)*6*sizeof(float), &reply_mask, 0, 0, 0);
                        // while(reply_mask!=1);
                    }
                }
#endif
                plane1 = plane_comp;
				plane2 = addn(plane_comp, 1, wy + y0 + y1 + 1);
				plane3 = addn(plane_comp, 2, wy + y0 + y1 + 1);
				plane4 = addn(plane_comp, 3, wy + y0 + y1 + 1);
				plane5 = addn(plane_comp, 4, wy + y0 + y1 + 1);
				plane =  plane3;

#ifdef RMA0
                readData_dma_with_rma4(rp,Frp+(j*slice_rp+(ixbeg0+1)*strip_rp+izbeg7+1)*3,
				Frp+(j*slice_rp+(ixbeg0+1)*strip_rp+izbeg3+1)*3, wx, wz*3, 0, 0, 0, 0, strip_rp*3,
                    len76*3, len75*3, len74*3, len32*3, len31*3, len30*3);
				readData_dma_with_rma4(vel,Fv+(j*slice+(ixbeg0+1)*strip+izbeg7+1)*3,
				Fv+(j*slice+(ixbeg0+1)*strip+izbeg3+1)*3, wx, wz*3, 0, 0, 0, 0, strip*3,
                    len76*3, len75*3, len74*3, len32*3, len31*3, len30*3);
#else
				if(izn > 0){
					for (i = ixbeg+1; i <= ixend; i++) {
						reply = 0;
						athread_get(PE_MODE, Frp+(j*slice_rp+i*strip_rp+izbeg+1)*3, rp+(i-ixbeg-1)*wz*3, izn*3*sizeof(float), &reply, 0, 0, 0);
						athread_get(PE_MODE, Fv+(j*slice+i*strip+izbeg+1)*3, vel+(i-ixbeg-1)*wz*3, izn*3*sizeof(float), &reply, 0, 0, 0);
						while(reply!=2);
					}
				}
#endif


				

                for (i = ixbeg+1; i <= ixend; i++) {
					for (k = izbeg+1; k <= izend; k++) {
					/// each fused array maintains its own index
                        
						idx = plane3 *slice_str + (i-ixbeg-1+x0)*strip_str + (k-izbeg-1+z0);
						idx_rp = (i-ixbeg-1)*wz+k-izbeg-1;
                        pos_jm2 = plane1 *slice_str + (i-ixbeg-1+x0)*strip_str + (k-izbeg-1+z0);
						pos_jm1 = plane2 *slice_str + (i-ixbeg-1+x0)*strip_str + (k-izbeg-1+z0);
						pos_jp2 = plane5 *slice_str + (i-ixbeg-1+x0)*strip_str + (k-izbeg-1+z0);
						pos_jp1 = plane4 *slice_str + (i-ixbeg-1+x0)*strip_str + (k-izbeg-1+z0);
						pos_im2 = idx - 2 * strip_str;
						pos_im1 = idx - strip_str;
						pos_ip2 = idx + 2 * strip_str;
						pos_ip1 = idx + strip_str;
						pos_km2 = idx - 2;
						pos_km1 = idx - 1;
						pos_kp2 = idx + 2;
						pos_kp1 = idx + 1;
						// if(Fs[ (j*slice+i*strip+k)*6+0] != strs[idx*6+0]){
						// 	printf("s 1 dma error\n");
						// }
						// if(Fs[((j-1)*slice+i*strip+k)*6+1] != strs[pos_jm1*6+1]){
						// 	printf("s 2 dma error\n");
						// }
						// if(Fs[(j*slice+(i+1)*strip+k)*6+5] != strs[pos_ip1*6+5]){
						// 	printf("s 3 dma error\n");
						// }
						// if(Frp[(j*slice_rp+i*strip_rp+k)*3+0] != rp[idx_rp*3+0]){
						// 	printf("rp  dma error\n");
						// }
						// if(Fv[(j*slice+i*strip+k)*3+2] != rp[idx_rp*3+2]){
						// 	printf("v dma error\n");
						// }
						
						// int idx = j * slice + i * strip + k;
						// int idx_rp = j * slice_rp + i * strip_rp + k;

						// sxx_x = dx * (b1 * (Fs[(idx + strip) * 6 + 0] - Fs[idx * 6 + 0]) + b2 * (Fs[(idx + 2 * strip) * 6 + 0] - Fs[(idx - strip) * 6 + 0]));
						// sxy_y = dy * (b1 * (Fs[idx * 6 + 3] - Fs[(idx - slice) * 6 + 3]) + b2 * (Fs[(idx + slice) * 6 + 3] - Fs[(idx - 2 * slice) * 6 + 3]));
						// sxz_z = dz * (b1 * (Fs[idx * 6 + 5] - Fs[(idx - 1) * 6 + 5]) + b2 * (Fs[(idx + 1) * 6 + 5] - Fs[(idx - 2) * 6 + 5]));

						// // updating components of particle velocities
						// Fv[idx * 3 + 0] += (sxx_x + sxy_y +sxz_z) / Frp[idx_rp * 3 + 0];

						// syy_y = dy * (b1 * (Fs[(idx + slice) * 6 + 1] - Fs[idx * 6 + 1]) + b2 * (Fs[(idx + 2 * slice) * 6 + 1] - Fs[(idx - slice) * 6 + 1]));
						// sxy_x = dx * (b1 * (Fs[idx * 6 + 3] - Fs[(idx - strip) * 6 + 3]) + b2 * (Fs[(idx + strip) * 6 + 3] - Fs[(idx - 2 * strip) * 6 + 3]));
						// syz_z = dz * (b1 * (Fs[idx * 6 + 4] - Fs[(idx - 1) * 6 + 4]) + b2 * (Fs[(idx + 1) * 6 + 4] - Fs[(idx - 2) * 6 + 4]));

						// Fv[idx * 3 + 1] += (syy_y + sxy_x + syz_z) / Frp[idx_rp * 3 + 1];

						// szz_z = dz * (b1 * (Fs[(idx + 1) * 6 + 2] - Fs[idx * 6 + 2]) + b2 * (Fs[(idx + 2) * 6 + 2] - Fs[(idx - 1) * 6 + 2]));
						// sxz_x = dx * (b1 * (Fs[idx * 6 + 5] - Fs[(idx - strip) * 6 + 5]) + b2 * (Fs[(idx + strip) * 6 + 5] - Fs[(idx - 2 * strip) * 6 + 5]));
						// syz_y = dy * (b1 * (Fs[idx * 6 + 4] - Fs[(idx - slice) * 6 + 4]) + b2 * (Fs[(idx + slice) * 6 + 4] - Fs[(idx - 2 * slice) * 6 + 4]));

						// Fv[idx * 3 + 2] += (szz_z + sxz_x + syz_y) / Frp[idx_rp * 3 + 2];
						sxx_x = dx * (b1 * (strs[pos_ip1 * 6 + 0] - strs[idx * 6 + 0]) + b2 * (strs[pos_ip2 * 6 + 0] - strs[pos_im1 * 6 + 0]));
						sxy_y = dy * (b1 * (strs[idx * 6 + 3] - strs[pos_jm1 * 6 + 3]) + b2 * (strs[pos_jp1 * 6 + 3] - strs[pos_jm2 * 6 + 3]));
						sxz_z = dz * (b1 * (strs[idx * 6 + 5] - strs[pos_km1 * 6 + 5]) + b2 * (strs[pos_kp1 * 6 + 5] - strs[pos_km2 * 6 + 5]));

						// updating components of particle velocities
						vel[idx_rp * 3 + 0] += (sxx_x + sxy_y +sxz_z) / rp[idx_rp * 3 + 0];

						syy_y = dy * (b1 * (strs[pos_jp1 * 6 + 1] - strs[idx * 6 + 1]) + b2 * (strs[pos_jp2 * 6 + 1] - strs[pos_jm1 * 6 + 1]));
						sxy_x = dx * (b1 * (strs[idx * 6 + 3] - strs[pos_im1 * 6 + 3]) + b2 * (strs[pos_ip1 * 6 + 3] - strs[pos_im2 * 6 + 3]));
						syz_z = dz * (b1 * (strs[idx * 6 + 4] - strs[pos_km1 * 6 + 4]) + b2 * (strs[pos_kp1 * 6 + 4] - strs[pos_km2 * 6 + 4]));

						vel[idx_rp * 3 + 1] += (syy_y + sxy_x + syz_z) / rp[idx_rp * 3 + 1];

						szz_z = dz * (b1 * (strs[pos_kp1 * 6 + 2] - strs[idx * 6 + 2]) + b2 * (strs[pos_kp2 * 6 + 2] - strs[pos_km1 * 6 + 2]));
						sxz_x = dx * (b1 * (strs[idx * 6 + 5] - strs[pos_im1 * 6 + 5]) + b2 * (strs[pos_ip1 * 6 + 5] - strs[pos_im2 * 6 + 5]));
						syz_y = dy * (b1 * (strs[idx * 6 + 4] - strs[pos_jm1 * 6 + 4]) + b2 * (strs[pos_jp1 * 6 + 4] - strs[pos_jm2 * 6 + 4]));

						vel[idx_rp * 3 + 2] += (szz_z + sxz_x + syz_y) / rp[idx_rp * 3 + 2];
					}//k
#ifdef RMAP 
                    athread_ssync (ROW_SCOPE, 0xff); 
					if(izn > 0 ){
						l_reply = 0;
                        r_reply = 0;
						athread_rma_iput(vel+(i-ixbeg-1)*wz*3,&l_reply,sizeof(float)*izn*3,_ROW*8+7,vel_change+(sum_izn-_COL*wz-izn)*3,&r_reply);
					    while(l_reply!=1);
					}                   
					athread_ssync (ROW_SCOPE, 0xff);
                    if(_COL == 7 && izn > 0){
						reply = 0;
						athread_put(PE_MODE, vel_change, Fv+(j*slice+i*strip+izbeg+1)*3, sizeof(float)*sum_izn*3, &reply, 0, 0);
					    while(reply!=1);
					}

#endif
				}//i
#ifdef RMAP
#else
				if(izn > 0){
					for (i = ixbeg+1; i <= ixend; i++) {
						reply = 0;
						athread_put(PE_MODE, vel+(i-ixbeg-1)*wz*3, Fv+(j*slice+i*strip+izbeg+1)*3, izn*3*sizeof(float), &reply, 0, 0);
						while(reply!=1);
					}
				}
#endif
#ifdef RMA
                readData_rma4_2(data_ldm, strs+plane_last*(wx+x0+x1)*(wz+z0+z1)*6, wx, wz*6, x0, x1, z0*6, z1*6, strip*6,
                    len76*6, len75*6, len74*6, len32*6, len31*6, len30*6);
#else
                if(izn > 0 && j < dim_y && ixn > 0){
                    while(reply_mask!=ixn+x1+x0);
                }
#endif
                plane_last = addn(plane_last, 1, wy + y0 + y1 + 1); 
                plane_comp = addn(plane_comp, 1, wy + y0 + y1 + 1); 
			}//j
		}//iix
	 }//iiz






}

