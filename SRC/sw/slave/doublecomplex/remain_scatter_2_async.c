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

__thread_local double SRC[MAX_BLOCK_SIZE * 4];
__thread_local double TARGET[MAX_BLOCK_SIZE * 4];


volatile __thread_local athread_rply_t put_reply = 0;
__thread_local int put_count = 0;
volatile __thread_local athread_rply_t get_target_reply[2];
volatile __thread_local athread_rply_t get_source_reply[2];

void  remain_scatter_2_async(remain_scatter_param_t* param){
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

    assert(ldt <= MAX_BLOCK_SIZE);

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

                // get next SRC
                int jj_ptr_next = 0;
                int fnz_next = fnzs[jj_ptr_next];
                doublecomplex *ucol_next = non_zero_ucols[jj_ptr_next];
                int source_flag_next = jj_ptr_next & 1;
                double* SRC_next = &SRC[MAX_BLOCK_SIZE * 2 * source_flag_next];
                get_source_reply[source_flag_next] = 0;
                athread_dma_iget(SRC_next, tempv + nbrow * jj_ptr_next, 2 * temp_nbrow * sizeof(double), &get_source_reply[source_flag_next]);
                
                jj_ptr_next += 1;
                for (; jj_ptr_next <= non_zero_col_num; ++jj_ptr_next) {
                    // prepare current SRC
                    int jj_ptr_cur = jj_ptr_next - 1;
                    int fnz_cur = fnz_next;
                    doublecomplex *ucol_cur = ucol_next;
                    int source_flag_cur = jj_ptr_cur & 1;
                    double* SRC_cur = &SRC[MAX_BLOCK_SIZE * 2 * source_flag_cur];
                    
                    // get next SRC
                    if(jj_ptr_next < non_zero_col_num){
                        fnz_next = fnzs[jj_ptr_next];
                        ucol_next = non_zero_ucols[jj_ptr_next];
                        source_flag_next = jj_ptr_next & 1;
                        SRC_next = &SRC[MAX_BLOCK_SIZE * 2 * source_flag_next];
                        get_source_reply[source_flag_next] = 0;
                        athread_dma_iget(SRC_next, tempv + nbrow * jj_ptr_next, 2 * temp_nbrow * sizeof(double), &get_source_reply[source_flag_next]);
                    }
                    // wait currnet SRC
                    athread_dma_wait_value (&get_source_reply[source_flag_cur], 1);

                    // process current SRC
                    {
                        int ptr_next = 0;
                        // get next TARGET
                        int i_start_next = segment_ptr[ptr_next];
                        int i_end_next = segment_ptr[ptr_next+1];
                        int i_len_next = (i_end_next - i_start_next) * 2;
                        int offset_next = segment_offset[ptr_next];
                        int rel_next = offset_next - fnz_cur;
                        int target_flag_next = ptr_next & 1;
                        double* TARGET_next = &TARGET[MAX_BLOCK_SIZE * 2 * target_flag_next];
                        get_target_reply[target_flag_next] = 0;
                        athread_dma_iget(TARGET_next, &ucol_cur[rel_next + i_start_next], i_len_next * sizeof(double), &get_target_reply[target_flag_next]);
                        ptr_next += 1;  

                        for(; ptr_next <= segment_count; ++ptr_next){
                            // prepare current TARGET
                            int ptr_cur = ptr_next - 1;
                            int i_len_cur = i_len_next;
                            int i_start_cur = i_start_next;
                            int rel_cur = rel_next;
                            int target_flag_cur = ptr_cur & 1;
                            double* TARGET_cur  = &TARGET[MAX_BLOCK_SIZE * 2 * target_flag_cur];
                            double *SRC_offset = &SRC_cur[i_start_cur * 2];

                            // get next TARGET
                            if(ptr_next < segment_count){
                                target_flag_next = ptr_next & 1;
                                i_start_next = segment_ptr[ptr_next];
                                i_end_next = segment_ptr[ptr_next+1];
                                i_len_next = (i_end_next - i_start_next) * 2;
                                offset_next = segment_offset[ptr_next];
                                rel_next = offset_next - fnz_cur;
                                TARGET_next = &TARGET[MAX_BLOCK_SIZE * 2 * target_flag_next];
                                get_target_reply[target_flag_next] = 0;
                                athread_dma_iget(TARGET_next, &ucol_cur[rel_next + i_start_next], i_len_next * sizeof(double), &get_target_reply[target_flag_next]);
                            }

                            // wait current TARGET
                            athread_dma_wait_value (&get_target_reply[target_flag_cur], 1);

                            // process current TARGET
                            for(int i = 0; i < i_len_cur; i++){
                                TARGET_cur[i] -= SRC_offset[i];
                            }
                            
                            // write current TARGET
                            // athread_dma_iput(&ucol_cur[rel_cur + i_start_cur], TARGET_cur, i_len_cur * sizeof(double), &put_reply);
                            // put_count += 1;
                            athread_dma_put(&ucol_cur[rel_cur + i_start_cur], TARGET_cur, i_len_cur * sizeof(double));
                        }
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
                        // athread_dma_iput(NZVAL, TARGET, i_len * sizeof(double), &put_reply);
                        // put_count += 1;
                        athread_dma_put(NZVAL, TARGET, i_len * sizeof(double));
                    }
                }
            }
        }   
    }
    // athread_dma_wait_value (&put_reply, put_count);
}
