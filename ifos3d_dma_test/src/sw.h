/*
 * macdrp_slave.h
 *
 *  Created on: Jun 12, 2018
 *      Author: bingo
 */


#define get_slv_id(tid) asm volatile ("rcsr %0, 0" : "=r"(tid))
#define get_row_id(rid) asm volatile ("rcsr %0, 1" : "=r"(rid))
#define get_col_id(cid) asm volatile ("rcsr %0, 2" : "=r"(cid))
#define REG_PUTR(var, dst) asm volatile ("putr %0,%1"::"r"(var),"r"(dst))
#define REG_PUTC(var, dst) asm volatile ("putc %0,%1"::"r"(var),"r"(dst))
#define REG_GETR(var) asm volatile ("getr %0":"=r"(var))
#define REG_GETC(var) asm volatile ("getc %0":"=r"(var))


// blk is the number of floats which are transfered for regsiter
// eg blk = 24, will transfer (wz + z0 + z1)*6 floats
#define simdcpy_z(src, dst, blk) {\
    int i, j; \
	floatv4 *src_p = (floatv4*)(src + cid*wz*6);\
	floatv4 *dst_p = (floatv4*)(dst + cid*(wz + z0 + z1)*6);\
	register floatv4 d0;\
	for (i = 0; i < blk; i++) {\
		simd_load(d0, src_p);\
		simd_store(d0, dst_p);\
		src_p ++;\
		dst_p ++;\
	}\
}

#define regcomm_put_z(RC, src, dst_id, blk) {\
    int i, j; \
	floatv4 *src_p = (floatv4*)(src);\
	register floatv4 d0;\
	for(i = 0; i < blk; i++) {\
		simd_load(d0, src_p);\
		REG_PUT##RC(d0, dst_id);\
		src_p ++;\
	}\
}

#define regcomm_get_z(RC, dst, blk) {\
    int i, j; \
	floatv4 *dst_p = (floatv4*)(dst);\
	register floatv4 d0;\
	for(i = 0; i < blk; i++) {\
		REG_GET##RC(d0);\
		simd_store(d0, dst_p);\
		dst_p ++;\
	}\
}

// eg blk = 24, will transfer (wz + z0 + z1)*6 floats
#define regshffle_strs(src, dst, blk) {\
    int i, j, k; \
	int target_id;\
	floatv4 *src_p, *dst_p;\
	register floatv4 d0;\
	simdcpy_z(src, dst, blk);\
	for(k = 1; k < 8; k++) {\
		target_id = cid ^ k;\
		src_p = (floatv4*)(src + target_id*wz*6);\
		dst_p = (floatv4*)(dst + target_id*(wz + z0 + z1)*6);\
	    for(i = 0; i < blk; i++) {\
		    if(cid < target_id) {\
		        simd_load(d0, src_p);\
		        REG_PUTR(d0, target_id);\
		        src_p++;\
		    } else {\
		        REG_GETR(d0);\
		        simd_store(d0, dst_p);\
		        dst_p++;\
		    }\
	    }\
		athread_syn(ROW_SCOPE, 0xff);\
	}\
	for(k = 1; k < 8; k++) {\
		target_id = cid ^ k;\
		src_p = (floatv4*)(src + target_id*wz*6);\
		dst_p = (floatv4*)(dst + target_id*(wz + z0 + z1)*6);\
	    for(i = 0; i < blk; i++) {\
		    if(cid > target_id) {\
		        simd_load(d0, src_p);\
		        REG_PUTR(d0, target_id);\
		        src_p++;\
		    } else {\
		        REG_GETR(d0);\
		        simd_store(d0, dst_p);\
		        dst_p++;\
		    }\
	    }\
		athread_syn(ROW_SCOPE, 0xff);\
	}\
}

#define simdcpy_z_halo(src, dst, blk) {\
    int i, j; \
	floatv4 *src_p = (floatv4*)(src + cid*wz*3);\
	floatv4 *dst_p = (floatv4*)(dst + cid*wz*3);\
	register floatv4 d0;\
	for (i = 0; i < blk; i++) {\
		simd_load(d0, src_p);\
		simd_store(d0, dst_p);\
		src_p ++;\
		dst_p ++;\
	}\
}

#define regshffle_vel(src, dst, blk) {\
    int i, j, k; \
	int target_id;\
	floatv4 *src_p, *dst_p;\
	register floatv4 d0;\
	simdcpy_z_halo(src, dst, blk);\
	for(k = 1; k < 8; k++) {\
		target_id = cid ^ k;\
		src_p = (floatv4*)(src + target_id*wz*3);\
		dst_p = (floatv4*)(dst + target_id*wz*3);\
	    for(i = 0; i < blk; i++) {\
		    if(cid < target_id) {\
		        simd_load(d0, src_p);\
		        REG_PUTR(d0, target_id);\
		        src_p++;\
		    } else {\
		        REG_GETR(d0);\
		        simd_store(d0, dst_p);\
		        dst_p++;\
		    }\
	    }\
		athread_syn(ROW_SCOPE, 0xff);\
	}\
	for(k = 1; k < 8; k++) {\
		target_id = cid ^ k;\
		src_p = (floatv4*)(src + target_id*wz*3);\
		dst_p = (floatv4*)(dst + target_id*wz*3);\
	    for(i = 0; i < blk; i++) {\
		    if(cid > target_id) {\
		        simd_load(d0, src_p);\
		        REG_PUTR(d0, target_id);\
		        src_p++;\
		    } else {\
		        REG_GETR(d0);\
		        simd_store(d0, dst_p);\
		        dst_p++;\
		    }\
	    }\
		athread_syn(ROW_SCOPE, 0xff);\
	}\
}

// eg blk = 24, will transfer wz*3 floats
#define regshffle_vel_bak(src, dst, blk) {\
    int k; \
	int target_id;\
	floatv4 *regshff_src_p, *regshff_dst_p;\
	simdcpy_z_halo(src, dst, blk);\
	for(k = 1; k < 8; k++) {\
		target_id = cid ^ k;\
		regshff_src_p = (floatv4*)(src + target_id*wz*3);\
		regshff_dst_p = (floatv4*)(dst + target_id*wz*3);\
		if(cid < target_id) {\
			regcomm_put_z(R, regshff_src_p, target_id, blk);\
		} else {\
			regcomm_get_z(R, regshff_dst_p, blk);\
		}\
		athread_syn(ROW_SCOPE, 0xff);\
	}\
	for(k = 1; k < 8; k++) {\
		target_id = cid ^ k;\
		regshff_src_p = (floatv4*)(src + target_id*wz*3);\
		regshff_dst_p = (floatv4*)(dst + target_id*wz*3);\
		if(cid > target_id) {\
			regcomm_put_z(R, regshff_src_p, target_id, blk);\
		} else {\
			regcomm_get_z(R, regshff_dst_p, blk);\
		}\
		athread_syn(ROW_SCOPE, 0xff);\
	}\
}


/*
#define simdcpy_z_halo(src, dst) {\
	int i, j;\
	floatv4 *src_p = (floatv4*)(src + cid * 12 * 9);\
	floatv4 *dst_p = (floatv4*)(dst + (MY + (cid >> 1)) * 16 * 9);\
	register floatv4 d0;\
	for (i = 0; i < 4; i++) {\
		for (j = 0; j < 9; j++) {\
			simd_load(d0, src_p);\
			simd_store(d0, dst_p);\
			src_p ++;\
			dst_p ++;\
		}\
	}\
}

#define regshuffle_halo(src, halo_side) {\
	int k;\
	floatv4 *regshff_src_p, *regshff_dst_p;\
	float *src_buf, *dst_buf, *tmp_buf;\
	src_buf = buffer_w_line; dst_buf = buffer_w_sq; tmp_buf = buffer_w_accept;\
	if (~cid & 1) {\
		dma(get_desc, (long)(src + (MY - (cid >> 1)) * ystep_w), (long)(tmp_buf));\
		dma_wait(&get_reply, 1);\
		get_reply = 0;\
		simdcpy_z_cid7(src_buf, dst_buf);\
	}\
	for (k = 1; k < 8; k++) {\
		int target_id = cid ^ k;\
			regshff_src_p = (floatv4*)(src_buf + target_id * 12 * wsize);\
			regshff_dst_p = (floatv4*)(dst_buf + (MY + (target_id >> 1)) * 16 * wsize);\
			if (target_id < cid) {\
				if (~cid & 1) regcomm_put_u16(R, regshff_src_p, target_id);\
			} else {\
				if (~target_id & 1) regcomm_get_u16(R, regshff_dst_p);\
			}\
			athread_syn(COL_SCOPE, 0xff);\
	}\
	for (k = 1; k < 8; k++) {\
		int target_id = cid ^ k;\
		regshff_src_p = (floatv4*)(src_buf + target_id * 12 * wsize);\
		regshff_dst_p = (floatv4*)(dst_buf + (MY + (target_id >> 1)) * 16 * wsize);\
		if (target_id > cid) {\
			if (~cid & 1) regcomm_put_u16(R, regshff_src_p, target_id);\
		} else {\
			if (~target_id & 1) regcomm_get_u16(R, regshff_dst_p);\
		}\
		athread_syn(COL_SCOPE, 0xff);\
	}\
}
*/


