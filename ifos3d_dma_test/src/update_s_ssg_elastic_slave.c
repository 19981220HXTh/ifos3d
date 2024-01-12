#include <slave.h>


#define wx 10
#define wy 1
#define MX 8
#define MZ 8
//#define DEBUG
#define RMA
//#define RMA0
#define RMAP

typedef struct Param_str {
	int dim_x;
	int dim_y;
	int dim_z;
	int slice;
	int strip;
	int slice_up;
	int strip_up;
	int slice_pi;
	int strip_pi;
	float DT;
	float DX;
	float DY;
	float DZ;
	int FDCOEFF;
	// float *vel;
	// float *stress;
	// float *up;
	// float *pi;
	// float *u;
	float *Fv;
	float *Fs;
	float *Fup;
	float *Fpi;
	float *Fu;
} Param_str;


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


void update_s_kernel_fusion_slave(Param_str *param) {

	int id = athread_get_id(-1);

	int dim_x = param->dim_x;
	int dim_y = param->dim_y;
	int dim_z = param->dim_z;
	int slice = param->slice;
	int strip = param->strip;
	int slice_up = param->slice_up;
	int strip_up = param->strip_up;
	int slice_pi = param->slice_pi;
	int strip_pi = param->strip_pi;
	float DT = param->DT;
	float DX = param->DX;
	float DY = param->DY;
	float DZ = param->DZ;
	int FDCOEFF = param->FDCOEFF;

	float *Fv = param->Fv;
	float *Fs = param->Fs;
	float *Fup = param->Fup;
	float *Fpi = param->Fpi;
	float *Fu = param->Fu;

	int wz = 5;
	wz = (dim_z <= wz) ? dim_z : wz; 


	//stencil operator
	const int x0 = 2;
	const int x1 = 2;
	const int y0 = 2;
	const int y1 = 2;
	const int z0 = 2;
	const int z1 = 2;

	int ix, iy, iz, iix, iiz;
	int ixbeg, ixend, ixn;
	int izbeg, izend, izn;
	int i,j,k,l;
	int izbeg0,izbeg1,izbeg2,izbeg3,izbeg4,izbeg5,izbeg6,izbeg7,ixbeg0;
	int len76,len75,len74,len32,len31,len30;
	int izend0,sum_izn;//put change



	// int xstep_6 = xstep * 6;
	// int ystep_6 = ystep * 6;
	// int xstep_3 = xstep * 3;
	// int ystep_3 = ystep * 3;
	// int xstep_l = wz + z0 + z1;
	// int ystep_l = (wx + x0 + x1) * (wz + z0 + z1);
	// int xstep_l2 = xstep_l * 2;
	// int ystep_l3 = ystep_l * 3;
	// int xstep_l3 = xstep_l * 3;
	// int ystep_l6 = ystep_l * 6;
	// int xstep_l6 = xstep_l * 6;



	// int NX = (dim_x + wx * MX - 1) / (wx * MX);
	// int NZ = (dim_z + wz - 1) / wz;
	int yid = id / (MX * MZ);
    int xid = (id - yid * (MX * MZ)) / MZ;
    int zid = (id - yid * (MX * MZ)) % MZ;
	int nx = (dim_x + wx * MX  - 1) / (wx * MX);
	int nz = (dim_z + wz * MZ - 1) / (wz * MZ);

	float vel[(wy + y0 + y1 + 1) * (wx + x0 + x1) * (wz + z0 + z1) * 3];
	float strs[wx * wz * 6];
	float up[wx * wz * 3];
	float pi[wx * wz];
	float u[wx * wz ];

	// float *vel_s1, *strs_s1, *up_s1, *pi_s1, *u_s1;
    // float *vel_temp, *vel_s1_temp, *strs_temp, *strs_s1_temp,*up_temp,*up_s1_temp,*pi_temp,*pi_s1_temp,*u_temp,*u_s1_temp;
	// float *vel_c = vel + x0 * xstep_l3 + z0 * 3;

	float vxx, vxy, vxz, vyx, vyy, vyz, vzx, vzy, vzz;
	float vxyyx, vyzzy, vxzzx, vxxyyzz, vyyzz, vxxzz, vxxyy;
	float g, f, fipjp, fjpkp, fipkp;
	float b1, b2;

    

    int idx,idx_p,idx_n,idx_test,slice_v,strip_v;
	int pos_jp1,pos_jp2,pos_jm1,pos_jm2,pos_ip1,pos_ip2,pos_im1,pos_im2,pos_kp1,pos_kp2,pos_km1,pos_km2;
    int plane,plane_comp,plane_last,plane1,plane2,plane3,plane4,plane5;

	//dma rma try2
    float data_ldm[ 1 * ((x0 + wx + x1)/4+1) * (z0 + wz*4 + z1)* 3];
    //dma rma try2


	b1 = 9.0 / 8.0; b2 = -1.0 / 24.0; /* Taylor coefficients*/
	if (FDCOEFF == 2) {
		b1 = 1.1382; b2 = -0.046414;
	} /* Holberg coefficients E=0.1 %*/

    slice_v = (wx+x0+x1)*(wz+z0+z1);
	strip_v = (wz+z0+z1);

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
		float strs_change[sum_izn*6];

		for (iix = 0; iix < nx; iix++) {
			ixbeg = wx * MX * iix + wx * xid;
			ixend = wx * MX * iix + wx * (xid + 1);
			ixend = ixend < dim_x ? ixend : dim_x;
			ixn = ixend - ixbeg;
			ixbeg0 = wx * MX * iix;
#ifdef RMA 
			for(l = 0;l<wy+y0+y1;l++){ 
				// eight_coop(Fv+((j+l-y0)*slice+(ixbeg0+1-x0)*strip+izbeg7+1-z0)*3,Fv+((j+l-y0)*slice+(ixbeg0+1-x0)*strip+izbeg0+1-z0)*3,
				// vel+l*(wx+x0+x1)*(wz+z0+z1)*3,strip*3,wz*3,x0,x1,z0*3,z1*3);
				readData_dma_with_rma4(vel+l*(wx+x0+x1)*(wz+z0+z1)*3,Fv+((1+l-y0)*slice+(ixbeg0+1-x0)*strip+izbeg7+1-z0)*3,
				Fv+((1+l-y0)*slice+(ixbeg0+1-x0)*strip+izbeg3+1-z0)*3, wx, wz*3, x0, x1, z0*3, z1*3, strip*3,
					len76*3, len75*3, len74*3, len32*3, len31*3, len30*3);
			}

#else
			if(izn > 0){
				for(l = 0;l<wy+y0+y1;l++){
					for (i = ixbeg+1; i <= ixend+x0+x1; i++) {
						reply = 0;
						athread_get(PE_MODE, Fv+((1+l-y0)*slice+(i-x0)*strip+izbeg+1-z0)*3, vel+l*(wx+x0+x1)*(wz+z0+z1)*3+(i-ixbeg-1)*(wz+z0+z1)*3, (izn+z0+z1)*3*sizeof(float), &reply, 0, 0, 0);
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
				readData_rma4_1(data_ldm,Fv+((j+wy+y1)*slice+(ixbeg0+1-x0)*strip+izbeg7+1-z0)*3,
				Fv+((j+wy+y1)*slice+(ixbeg0+1-x0)*strip+izbeg3+1-z0)*3, wx, wz*3, x0, x1, z0*3, z1*3, strip*3,
					len74*3, len30*3);
#else
                reply_mask = 0;
                if(izn > 0 && j < dim_y && ixn > 0){
                    for (i = ixbeg+1; i <= ixend+x0+x1; i++) {
                        // reply_mask = 0;
                        athread_get(PE_MODE,  Fv+((j+wy+y1)*slice+(i-x0)*strip+izbeg+1-z0)*3, vel+plane_last*(wx+x0+x1)*(wz+z0+z1)*3+(i-ixbeg-1)*(wz+z0+z1)*3, (izn+z0+z1)*3*sizeof(float), &reply_mask, 0, 0, 0);
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
                readData_dma_with_rma4(up,Fup+(j*slice_up+(ixbeg0+1)*strip_up+izbeg7+1)*3,
				Fup+(j*slice_up+(ixbeg0+1)*strip_up+izbeg3+1)*3, wx, wz*3, 0, 0, 0, 0, strip_up*3,
                    len76*3, len75*3, len74*3, len32*3, len31*3, len30*3);
				readData_dma_with_rma4(strs,Fs+(j*slice+(ixbeg0+1)*strip+izbeg7+1)*6,
				Fs+(j*slice+(ixbeg0+1)*strip+izbeg3+1)*6, wx, wz*6, 0, 0, 0, 0, strip*6,
                    len76*6, len75*6, len74*6, len32*6, len31*6, len30*6);
				readData_dma_with_rma4(pi,Fpi+(j*slice_pi+(ixbeg0+1)*strip_pi+izbeg7+1),
				Fpi+(j*slice_pi+(ixbeg0+1)*strip_pi+izbeg3+1), wx, wz, 0, 0, 0, 0, strip_pi,
                    len76, len75, len74, len32, len31, len30);
				readData_dma_with_rma4(u,Fu+(j*slice_pi+(ixbeg0+1)*strip_pi+izbeg7+1),
				Fu+(j*slice_pi+(ixbeg0+1)*strip_pi+izbeg3+1), wx, wz, 0, 0, 0, 0, strip_pi,
                    len76, len75, len74, len32, len31, len30);
#else
				if(izn > 0){
					for (i = ixbeg+1; i <= ixend; i++) {
						reply = 0;
						athread_get(PE_MODE, Fup+(j*slice_up+i*strip_up+izbeg+1)*3, up+(i-ixbeg-1)*wz*3, izn*3*sizeof(float), &reply, 0, 0, 0);
						athread_get(PE_MODE, Fpi+(j*slice_pi+i*strip_pi+izbeg+1), pi+(i-ixbeg-1)*wz, izn*sizeof(float), &reply, 0, 0, 0);
						athread_get(PE_MODE, Fu+(j*slice_pi+i*strip_pi+izbeg+1), u+(i-ixbeg-1)*wz, izn*sizeof(float), &reply, 0, 0, 0);
						athread_get(PE_MODE, Fs+(j*slice+i*strip+izbeg+1)*6, strs+(i-ixbeg-1)*wz*6, izn*6*sizeof(float), &reply, 0, 0, 0);
						while(reply!=4);
					}
				}
#endif

				for (i = ixbeg+1; i <= ixend; i++) {
					for (k = izbeg+1; k <= izend; k++) {
						// int idx = j * slice + i * strip + k;
						// int idx_n = j * slice_up + i * strip_up + k;
						// int idx_p = j * slice_pi + i * strip_pi + k;
						// vxx = (b1 * (Fv[idx * 3 + 0] - Fv[(idx - strip) * 3 + 0]) + b2 * (Fv[(idx + strip) * 3 + 0] - Fv[(idx - 2 * strip) * 3 + 0])) / DX;
						// vxy = (b1 * (Fv[(idx + slice) * 3 + 0] - Fv[idx * 3 + 0]) + b2 * (Fv[(idx + 2 * slice) * 3 + 0] - Fv[(idx - slice) * 3 + 0])) / DY;
						// vxz = (b1 * (Fv[(idx + 1) * 3 + 0] - Fv[idx * 3 + 0]) + b2 * (Fv[(idx + 2) * 3 + 0] - Fv[(idx -1) * 3 + 0])) / DZ;
						// vyx = (b1 * (Fv[(idx + strip) * 3 + 1] - Fv[idx * 3 + 1]) + b2 * (Fv[(idx + 2 * strip) * 3 + 1] - Fv[(idx - strip) * 3 + 1])) / DX;
						// vyy = (b1 * (Fv[idx * 3 + 1] - Fv[(idx - slice) * 3 + 1]) + b2 * (Fv[(idx + slice) * 3 + 1] - Fv[(idx - 2 * slice) * 3 + 1])) / DY;
						// vyz = (b1 * (Fv[(idx + 1) * 3 + 1] - Fv[idx * 3 + 1]) + b2 * (Fv[(idx + 2) * 3 + 1] - Fv[(idx - 1) * 3 + 1])) / DZ;
						// vzx = (b1 * (Fv[(idx + strip) * 3 + 2] - Fv[idx * 3 + 2]) + b2 * (Fv[(idx + 2 * strip) * 3 + 2] - Fv[(idx - strip) * 3 + 2])) / DX;
						// vzy = (b1 * (Fv[(idx + slice) * 3 + 2] - Fv[idx * 3 + 2]) + b2 * (Fv[(idx + 2 * slice) * 3 + 2] - Fv[(idx - slice) * 3 + 2])) / DY;
						// vzz = (b1 * (Fv[idx * 3 + 2] - Fv[(idx - 1) * 3 + 2]) + b2 * (Fv[(idx + 1) * 3 + 2] - Fv[(idx - 2) * 3 + 2])) / DZ;


						// // updating components of the stress tensor, partially
						// fipjp=Fup[idx_n * 3 + 0] * DT;
						// fjpkp=Fup[idx_n * 3 + 1] * DT;
						// fipkp=Fup[idx_n * 3 + 2] * DT;
						// g = Fpi[idx_p];
						// f = 2.0*Fu[idx_p];


						// vxyyx=vxy+vyx;
						// vyzzy=vyz+vzy;
						// vxzzx=vxz+vzx;
						// vxxyyzz=vxx+vyy+vzz;
						// vyyzz=vyy+vzz;
						// vxxzz=vxx+vzz;
						// vxxyy=vxx+vyy;


						// Fs[6 * idx + 0] += DT*((g*vxxyyzz)-(f*vyyzz));
						// Fs[6 * idx + 1] += DT*((g*vxxyyzz)-(f*vxxzz));
						// Fs[6 * idx + 2] += DT*((g*vxxyyzz)-(f*vxxyy));
						// Fs[6 * idx + 3] += (fipjp*vxyyx);
						// Fs[6 * idx + 4] += (fjpkp*vyzzy);
						// Fs[6 * idx + 5] += (fipkp*vxzzx);


						idx = plane3 *slice_v + (i-ixbeg-1+x0)*strip_v + (k-izbeg-1+z0);
						idx_n = (i-ixbeg-1) * wz + k-izbeg-1;
						// idx_test = j * slice + i * strip + k;
						pos_jm2 = plane1 *slice_v + (i-ixbeg-1+x0)*strip_v + (k-izbeg-1+z0);
						pos_jm1 = plane2 *slice_v + (i-ixbeg-1+x0)*strip_v + (k-izbeg-1+z0);
						pos_jp2 = plane5 *slice_v + (i-ixbeg-1+x0)*strip_v + (k-izbeg-1+z0);
						pos_jp1 = plane4 *slice_v + (i-ixbeg-1+x0)*strip_v + (k-izbeg-1+z0);
						pos_im2 = idx - 2 * strip_v;
						pos_im1 = idx - strip_v;
						pos_ip2 = idx + 2 * strip_v;
						pos_ip1 = idx + strip_v;
						pos_km2 = idx - 2;
						pos_km1 = idx - 1;
						pos_kp2 = idx + 2;
						pos_kp1 = idx + 1;

						// if(Fs[6*(j*slice+i*strip+k)+0] != strs[6*((i-ixbeg-1) * wz + k-izbeg-1) + 0]){
						// 	printf("s 0 dma error Fs:%f  strs:%f \n",Fs[6*(j*slice+i*strip+k)+0],strs[6*idx_test + 0]);
						// }
						// if(Fs[6*(j*slice+i*strip+k)+4] != strs[6*((i-ixbeg-1) * wz + k-izbeg-1)+ 4]){
						// 	printf("s 4 dma error\n");
						// }
						vxx = (b1 * (vel[idx * 3 + 0] - vel[pos_im1 * 3 + 0]) + b2 * (vel[pos_ip1 * 3 + 0] - vel[pos_im2 * 3 + 0])) / DX;
						vxy = (b1 * (vel[pos_jp1 * 3 + 0] - vel[idx * 3 + 0]) + b2 * (vel[pos_jp2 * 3 + 0] - vel[pos_jm1 * 3 + 0])) / DY;
						vxz = (b1 * (vel[pos_kp1 * 3 + 0] - vel[idx * 3 + 0]) + b2 * (vel[pos_kp2 * 3 + 0] - vel[pos_km1 * 3 + 0])) / DZ;
						vyx = (b1 * (vel[pos_ip1 * 3 + 1] - vel[idx * 3 + 1]) + b2 * (vel[pos_ip2 * 3 + 1] - vel[pos_im1 * 3 + 1])) / DX;
						vyy = (b1 * (vel[idx * 3 + 1] - vel[pos_jm1 * 3 + 1]) + b2 * (vel[pos_jp1 * 3 + 1] - vel[pos_jm2 * 3 + 1])) / DY;
						vyz = (b1 * (vel[pos_kp1 * 3 + 1] - vel[idx * 3 + 1]) + b2 * (vel[pos_kp2 * 3 + 1] - vel[pos_km1 * 3 + 1])) / DZ;
						vzx = (b1 * (vel[pos_ip1 * 3 + 2] - vel[idx * 3 + 2]) + b2 * (vel[pos_ip2 * 3 + 2] - vel[pos_im1* 3 + 2])) / DX;
						vzy = (b1 * (vel[pos_jp1 * 3 + 2] - vel[idx * 3 + 2]) + b2 * (vel[pos_jp2 * 3 + 2] - vel[pos_jm1 * 3 + 2])) / DY;
						vzz = (b1 * (vel[idx * 3 + 2] - vel[pos_km1 * 3 + 2]) + b2 * (vel[pos_kp1 * 3 + 2] - vel[pos_km2 * 3 + 2])) / DZ;


						// updating components of the stress tensor, partially
						fipjp=up[idx_n * 3 + 0] * DT;
						fjpkp=up[idx_n * 3 + 1] * DT;
						fipkp=up[idx_n * 3 + 2] * DT;
						g = pi[idx_n];
						f = 2.0*u[idx_n];


						vxyyx=vxy+vyx;
						vyzzy=vyz+vzy;
						vxzzx=vxz+vzx;
						vxxyyzz=vxx+vyy+vzz;
						vyyzz=vyy+vzz;
						vxxzz=vxx+vzz;
						vxxyy=vxx+vyy;


						strs[6 * idx_n + 0] += DT*((g*vxxyyzz)-(f*vyyzz));
						strs[6 * idx_n + 1] += DT*((g*vxxyyzz)-(f*vxxzz));
						strs[6 * idx_n + 2] += DT*((g*vxxyyzz)-(f*vxxyy));
						strs[6 * idx_n + 3] += (fipjp*vxyyx);
						strs[6 * idx_n + 4] += (fjpkp*vyzzy);
						strs[6 * idx_n + 5] += (fipkp*vxzzx);


					    // if(Fup[(j*slice_up+i*strip_up+k)*3+0] != up[idx_n*3+0]){
						// 	printf("1 dma error");
						// }
						// if(Fup[(j*slice_up+i*strip_up+k)*3+1] != up[idx_n*3+1]){
						// 	printf("2 dma error");
						// }
						// if(Fup[(j*slice_up+i*strip_up+k)*3+2] != up[idx_n*3+2]){
						// 	printf("3 dma error");
						// }
						// if(Fpi[j*slice_pi+i*strip_pi+k] != pi[idx_p]){
						// 	printf("pi dma error");
						// }
						// if(Fu[j*slice_pi+i*strip_pi+k] != u[idx_p]){
						// 	printf("u dma error");
						// }
						// if(Fv[ (j*slice+i*strip+k)*3+0] != vel[idx*3+0]){
						// 	printf("v 1 dma error\n");
						// }
						// if(Fv[((j-1)*slice+i*strip+k)*3+1] != vel[pos_jm1*3+1]){
						// 	printf("v 2 dma error\n");
						// }
						// if(Fv[(j*slice+(i+1)*strip+k)*3+5] != vel[(idx+strip_v)*3+5]){
						// 	printf("v 3 dma error\n");
						// }
					}//k
#ifdef RMAP 
                    athread_ssync (ROW_SCOPE, 0xff); 
					if(izn > 0){
						l_reply = 0;
						r_reply = 0;
						athread_rma_iput(strs+(i-ixbeg-1)*wz*6,&l_reply,sizeof(float)*izn*6,_ROW*8+7,strs_change+(sum_izn-_COL*wz-izn)*6,&r_reply);
					    while(l_reply!=1);
					}                   	
					athread_ssync (ROW_SCOPE, 0xff);
                    if(_COL == 7){
						reply = 0;
					    athread_put(PE_MODE, strs_change, Fs+(j*slice+i*strip+izbeg+1)*6, sizeof(float)*sum_izn*6, &reply, 0, 0);
						while(reply!=1);
					}
#endif
				}//i
#ifdef RMAP
#else
				if(izn > 0){
					for (i = ixbeg+1; i <= ixend; i++) {
						reply = 0;
						athread_put(PE_MODE, strs+(i-ixbeg-1)*wz*6, Fs+(j*slice+i*strip+izbeg+1)*6, izn*6*sizeof(float), &reply, 0, 0);
						while(reply!=1);
					}
				}
#endif

#ifdef RMA
            readData_rma4_2(data_ldm, vel+plane_last*(wx+x0+x1)*(wz+z0+z1)*3, wx, wz*3, x0, x1, z0*3, z1*3, strip*3,
					len76*3, len75*3, len74*3, len32*3, len31*3, len30*3);
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
