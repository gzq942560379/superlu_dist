#include <crts.h>
#include <stdio.h>
#include "sw/slave_param.h"

#ifndef MAX_BLOCK_SIZE
#define MAX_BLOCK_SIZE 256
#endif

__thread_local int indirect[MAX_BLOCK_SIZE];
__thread_local int indirect2[MAX_BLOCK_SIZE];
__thread_local int segment_ptr[MAX_BLOCK_SIZE+1];
__thread_local int segment_offset[MAX_BLOCK_SIZE];
__thread_local int segment_count;
__thread_local int non_zero_col_num;
__thread_local doublecomplex* nzvals[MAX_BLOCK_SIZE];
__thread_local doublecomplex * non_zero_ucols[MAX_BLOCK_SIZE];
__thread_local int fnzs[MAX_BLOCK_SIZE];

__thread_local double SRC[MAX_BLOCK_SIZE * 2];
__thread_local doublecomplex DEST[MAX_BLOCK_SIZE];

void  remain_scatter_2_index_compress_dma(remain_scatter_param_t* param){
    int l_remain_block_num = param->l_remain_block_num;
    int u_block_num = param->u_block_num;
    int u_start = param->u_start;
    int ldt = param->ldt;
    Ublock_info_t* Ublock_info = param->Ublock_info;
    Remain_info_t* Remain_info = param->Remain_info;
    doublecomplex* bigV = param->bigV;
    int* xsup = param->xsup;
    int klst = param->klst;
    int nbrow = param->nbrow;
    int_t* usub = param->usub;
    int_t* lsub = param->lsub;
    int_t** Ufstnz_br_ptr = param->Ufstnz_br_ptr;
    doublecomplex** Unzval_br_ptr = param->Unzval_br_ptr;
    int_t** Lrowind_bc_ptr = param->Lrowind_bc_ptr;
    doublecomplex** Lnzval_bc_ptr = param->Lnzval_bc_ptr;
    gridinfo_t* grid = param->grid;

    int local_j_start = (u_block_num * CRTS_cid)/8 + u_start;
    int local_j_end = (u_block_num * (CRTS_cid+1))/8 + u_start;
    int local_i_start = (l_remain_block_num * CRTS_rid)/8;
    int local_i_end = (l_remain_block_num * (CRTS_rid+1))/8;

    for(int j = local_j_start;j<local_j_end;j++){
        for(int lb = local_i_start;lb < local_i_end;lb++){
            // Getting U block information
            int_t rukp =  Ublock_info[j].rukp;
			int_t iukp =  Ublock_info[j].iukp;
			int jb   =  Ublock_info[j].jb;
			int nsupc = SuperSize(jb);
			int ljb = LBj (jb, grid);
			int st_col = j == u_start ? 0 : Ublock_info[j-1].full_u_cols;

			/* Getting L block L(i,k) information */
			int_t lptr = Remain_info[lb].lptr;
			int ib   = Remain_info[lb].ib;
			int temp_nbrow = lsub[lptr+1];
			lptr += LB_DESCRIPTOR;
			int cum_nrow = lb == 0 ? 0 : Remain_info[lb-1].FullRow;
			doublecomplex* tempv = bigV + (st_col * nbrow + cum_nrow); /* Sherry */
            // scatter_u
			if ( ib < jb ) {
                int_t ilst = FstBlockC (ib + 1);
                int_t lib = LBi (ib, grid);
                int_t *index = Ufstnz_br_ptr[lib];

                /* Reinitilize the pointers to the beginning of the k-th column/row of
                * L/U factors.
                * usub[] - index array for panel U(k,:)
                */
                int_t iuip_lib = BR_HEADER;
                int_t ruip_lib = 0;

                int_t ijb = index[iuip_lib];
                while (ijb < jb) {   /* Search for destination block. */
                    ruip_lib += index[iuip_lib + 1];
                    iuip_lib += UB_DESCRIPTOR + SuperSize (ijb);
                    ijb = index[iuip_lib];
                }
                iuip_lib += UB_DESCRIPTOR;

                int_t* indirect_index = &lsub[lptr];
                int len = temp_nbrow;
                // segment compress
                segment_ptr[0] = 0;
                segment_offset[0] = indirect_index[0];
                int segment_prev = indirect_index[0]; 
                int segment_count = 1;
                for(int i = 1; i < len; i++){
                    int segment_cur = indirect_index[i] - i;
                    if(segment_prev != segment_cur){
                        segment_ptr[segment_count] = i;
                        segment_offset[segment_count] = segment_cur;
                        segment_count += 1;
                        segment_prev = segment_cur;
                    }
                }
                segment_ptr[segment_count] = len;

                non_zero_col_num = 0;
                for (int jj = 0; jj < nsupc; ++jj) {
                    int_t segsize = klst - usub[iukp + jj];
                    int fnz = index[iuip_lib+jj];
                    if(segsize){
                        fnzs[non_zero_col_num] = fnz;
                        non_zero_ucols[non_zero_col_num] = &Unzval_br_ptr[lib][ruip_lib];
                        non_zero_col_num += 1;
                    }
                    ruip_lib += ilst - fnz;
                }

                for (int jj_ptr = 0; jj_ptr < non_zero_col_num; ++jj_ptr) {
                    int fnz = fnzs[jj_ptr];
                    doublecomplex *ucol_cur = non_zero_ucols[jj_ptr];
                    double *tempv_cur = (double*)(tempv + nbrow * jj_ptr);
                    int dest_nbrow = ilst - fnz;
                    CRTS_dma_get(SRC, tempv_cur, temp_nbrow * sizeof(doublecomplex));
                    CRTS_dma_get(DEST, ucol_cur, dest_nbrow * sizeof(doublecomplex));
                    for(int ptr = 0; ptr < segment_count; ++ptr){
                        int i_start = segment_ptr[ptr];
                        int i_end = segment_ptr[ptr+1];
                        int offset = segment_offset[ptr];
                        int rel = offset - fnz;
                        double *UCOL = (double*)&ucol_cur[rel];
                        #pragma omp simd
                        for(int i = i_start*2; i < i_end*2; i++){
                            UCOL[i] -= SRC[i];
                        }
                    }
                    CRTS_dma_put(ucol_cur, DEST, dest_nbrow * sizeof(doublecomplex));
                }

            // scatter_l
            }else{
                int_t *index = Lrowind_bc_ptr[ljb];
                int_t ldv = index[1];       /* LDA of the destination lusup. */
                int_t lptrj = BC_HEADER;
                int_t luptrj = 0;
                int_t ijb = index[lptrj];

                while (ijb != ib)  /* Search for destination block L(i,j) */
                {
                    luptrj += index[lptrj + 1];
                    lptrj += LB_DESCRIPTOR + index[lptrj + 1];
                    ijb = index[lptrj];
                }

                int_t fnz = FstBlockC (ib);
                lptrj += LB_DESCRIPTOR;
                int_t dest_nbrow=index[lptrj - 1];

                for (int i = 0; i < dest_nbrow; ++i) {
                    int_t rel = index[lptrj + i] - fnz;
                    indirect[rel] = i;
                }

                /* can be precalculated? */
                for (int i = 0; i < temp_nbrow; ++i) { /* Source index is a subset of dest. */
                    int_t rel = lsub[lptr + i] - fnz;
                    indirect2[i] = indirect[rel];
                }
                
                int* indirect_index = indirect2;
                int len = temp_nbrow;

                segment_ptr[0] = 0;
                segment_offset[0] = indirect_index[0];
                int segment_prev = indirect_index[0]; 
                int segment_count = 1;
                for(int i = 1; i < len; i++){
                    int segment_cur = indirect_index[i] - i;
                    if(segment_prev != segment_cur){
                        segment_ptr[segment_count] = i;
                        segment_offset[segment_count] = segment_cur;
                        segment_count += 1;
                        segment_prev = segment_cur;
                    }
                }
                segment_ptr[segment_count] = len;

                doublecomplex* nzval = Lnzval_bc_ptr[ljb] + luptrj; /* Destination block L(i,j) */

                non_zero_col_num = 0;
                for (int jj = 0; jj < nsupc; ++jj) {
                    int_t segsize = klst - usub[iukp + jj];
                    if(segsize){
                        nzvals[non_zero_col_num] = nzval + jj * ldv;
                        non_zero_col_num += 1;
                    }
                }
                for (int jj_ptr = 0; jj_ptr < non_zero_col_num; ++jj_ptr) {
                    doublecomplex *nzval_cur = nzvals[jj_ptr];
                    double *tempv_cur = (double*)(tempv + jj_ptr * nbrow);
                    CRTS_dma_get(SRC, tempv_cur, temp_nbrow * sizeof(doublecomplex));
                    CRTS_dma_get(DEST, nzval_cur, dest_nbrow * sizeof(doublecomplex));
                    for(int ptr = 0; ptr < segment_count; ++ptr){
                        int i_start = segment_ptr[ptr];
                        int i_end = segment_ptr[ptr+1];
                        int offset = segment_offset[ptr];
                        double *NZVAL = (double*)(DEST + offset);
                        for(int i = i_start*2; i < i_end*2; i++){
                            NZVAL[i] -= SRC[i];
                        } 
                    }
                    CRTS_dma_put(nzval_cur, DEST, dest_nbrow * sizeof(doublecomplex));
                }
            }
        }   
    }
}