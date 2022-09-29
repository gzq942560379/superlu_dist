#include <slave.h>
#include <stdio.h>
#include "sw/slave_param.h"

#ifndef MAX_BLOCK_SIZE 
#define MAX_BLOCK_SIZE 256
#endif

__thread_local int segment_ptr[MAX_BLOCK_SIZE + 1];
__thread_local int segment_offset[MAX_BLOCK_SIZE];

__thread_local int indirect[MAX_BLOCK_SIZE];
__thread_local int indirect2[MAX_BLOCK_SIZE];

__thread_local doublecomplex * non_zero_ucols[MAX_BLOCK_SIZE];
__thread_local int fnzs[MAX_BLOCK_SIZE];
__thread_local doublecomplex* nzvals[MAX_BLOCK_SIZE];

__thread_local double SRC[MAX_BLOCK_SIZE * 2];
__thread_local double TARGET[MAX_BLOCK_SIZE * 2];

void  remain_scatter_1_ldm(remain_scatter_param_t* param){
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

    int local_j_start = (u_block_num * _COL)/8 + u_start;
    int local_j_end = (u_block_num * (_COL+1))/8 + u_start;
    int local_i_start = (l_remain_block_num * _ROW)/8;
    int local_i_end = (l_remain_block_num * (_ROW+1))/8;

    for(int j = local_j_start;j<local_j_end;j++){
        for(int i = local_i_start;i < local_i_end;i++){
            int_t rukp =  Ublock_info[j].rukp;
			int_t iukp =  Ublock_info[j].iukp;
			int jb   =  Ublock_info[j].jb;
			int nsupc = SuperSize(jb);
			int ljb = LBj (jb, grid);
			int st_col;
			int ncols;
			if ( j>u_start ) {
				ncols = Ublock_info[j].full_u_cols - Ublock_info[j-1].full_u_cols;
				st_col = Ublock_info[j-1].full_u_cols;
			} else {
				ncols = Ublock_info[j].full_u_cols;
				st_col = 0;
			}
			/* Getting L block L(i,k) information */
			int_t lptr = Remain_info[i].lptr;
			int ib   = Remain_info[i].ib;
			int temp_nbrow = lsub[lptr+1];
			lptr += LB_DESCRIPTOR;
			int cum_nrow = (i==0 ? 0 : Remain_info[i-1].FullRow);

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

                assert(temp_nbrow > 0), "Assert temp_nbrow > 0 in ", __FILE__, ":", __LINE__;

                int segment_count = 1;
                int segment_prev = lsub[lptr];
                for(int i = 1; i < temp_nbrow; i++){
                    int segment_cur = lsub[lptr + i] - i;
                    if(segment_prev != segment_cur){
                        segment_count += 1;
                        segment_prev = segment_cur;
                    }
                }
                segment_ptr[0] = 0;
                segment_offset[0] = lsub[lptr];
                segment_prev = lsub[lptr];
                int segment_index = 1;
                for(int i = 1; i < temp_nbrow; i++){
                    int segment_cur = lsub[lptr + i] - i;
                    if(segment_prev != segment_cur){
                        segment_ptr[segment_index] = i;
                        segment_offset[segment_index] = segment_cur;
                        segment_index += 1;
                        segment_prev = segment_cur;
                    }
                }
                segment_ptr[segment_index] = temp_nbrow;
                assert(segment_index == segment_count);

                int non_zero_col_num = 0;

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
                    athread_dma_get(SRC, tempv + nbrow * jj_ptr, 2 * temp_nbrow * sizeof(double));
                    for(int ptr = 0; ptr < segment_count; ++ptr){
                        int i_start = segment_ptr[ptr];
                        int i_end = segment_ptr[ptr+1];
                        int i_len = (i_end - i_start) * 2;
                        int offset = segment_offset[ptr];
                        int rel = offset - fnz;
                        // double *UCOL = (double*)&ucol_cur[rel + i_start];
                        double *SRC_ = &SRC[i_start * 2];
                        athread_dma_get(TARGET, &ucol_cur[rel + i_start], i_len * sizeof(double));
                        for(int i = 0; i < i_len; i++){
                            TARGET[i] -= SRC_[i];
                        }
                        athread_dma_put(&ucol_cur[rel + i_start], TARGET, i_len * sizeof(double));
                    }
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

                /*
                * Build indirect table. This is needed because the indices are not sorted
                * in the L blocks.
                */
                int_t fnz = FstBlockC (ib);
                lptrj += LB_DESCRIPTOR;
                int_t dest_nbrow=index[lptrj - 1];

                assert(dest_nbrow >= temp_nbrow), "Assert : dest_nbrow >= temp_nbrow in ",__FILE__, ":",__LINE__;
                for (int i = 0; i < dest_nbrow; ++i) {
                    int rel = index[lptrj + i] - fnz;
                    indirect[rel] = i;
                }
                /* can be precalculated? */
                for (int i = 0; i < temp_nbrow; ++i) { /* Source index is a subset of dest. */
                    int rel = lsub[lptr + i] - fnz;
                    indirect2[i] = indirect[rel];
                }

                // compress indirect2 to segment
                int segment_count = 1;
                int segment_prev = indirect2[0];
                for(int i = 1; i < temp_nbrow; i++){
                    int segment_cur = indirect2[i] - i;
                    if(segment_prev != segment_cur){
                        segment_count += 1;
                        segment_prev = segment_cur;
                    }
                }
                segment_ptr[0] = 0;
                segment_offset[0] = indirect2[0];
                segment_prev = indirect2[0]; 
                int segment_index = 1;
                for(int i = 1; i < temp_nbrow; i++){
                    int segment_cur = indirect2[i] - i;
                    if(segment_prev != segment_cur){
                        segment_ptr[segment_index] = i;
                        segment_offset[segment_index] = segment_cur;
                        segment_index += 1;
                        segment_prev = segment_cur;
                    }
                }
                segment_ptr[segment_index] = temp_nbrow;
                assert(segment_index == segment_count);

                doublecomplex *nzval = Lnzval_bc_ptr[ljb] + luptrj; /* Destination block L(i,j) */
                
                int non_zero_col_num = 0;
                for (int jj = 0; jj < nsupc; ++jj) {
                    int_t segsize = klst - usub[iukp + jj];
                    if(segsize){
                        nzvals[non_zero_col_num] = nzval + jj * ldv;
                        non_zero_col_num += 1;
                    }
                }

                for (int jj_ptr = 0; jj_ptr < non_zero_col_num; ++jj_ptr) {
                    doublecomplex *nzval_cur = nzvals[jj_ptr];
                    // doublecomplex *tempv_cur = tempv + jj_ptr * nbrow;
                    athread_dma_get(SRC, tempv + jj_ptr * nbrow, 2 * temp_nbrow * sizeof(double));
                    for(int ptr = 0; ptr < segment_count; ++ptr){
                        int i_start = segment_ptr[ptr];
                        int i_end = segment_ptr[ptr+1];
                        int i_len = (i_end - i_start) * 2;
                        int offset = segment_offset[ptr];
                        double *NZVAL = (double*)(nzval_cur + offset + i_start);
                        double* SRC_ = &SRC[i_start*2];
                        athread_dma_get(TARGET, NZVAL, i_len * sizeof(double));
                        for(int i = 0; i < i_len; i++){
                            TARGET[i] -= SRC_[i];
                        } 
                        athread_dma_put(NZVAL, TARGET, i_len * sizeof(double));
                    }
                }
            }
        }   
    }
}
