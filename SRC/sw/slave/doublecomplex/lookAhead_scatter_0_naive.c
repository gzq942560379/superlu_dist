#include <slave.h>
#include <stdio.h>
#include "sw/slave_param.h"

void  lookAhead_scatter_0_naive(lookAhead_scatter_param_t* param){
    int l_lookAhead_block_num = param->l_lookAhead_block_num;
    int u_block_num = param->u_block_num;
    int u_start = param->u_start;
    int ldt = param->ldt;
    Ublock_info_t* Ublock_info = param->Ublock_info;
    int* lookAhead_lptr = param->lookAhead_lptr;
    int* lookAhead_ib = param->lookAhead_ib;
    int* lookAheadFullRow = param->lookAheadFullRow;
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

    int* indirect = ldm_malloc(ldt*sizeof(int));
    int* indirect2 = ldm_malloc(ldt*sizeof(int));

    int local_ij_start = (l_lookAhead_block_num * u_block_num * _MYID)/64;
    int local_ij_end = (l_lookAhead_block_num * u_block_num * (_MYID+1))/64;

    for(int ij = local_ij_start; ij<local_ij_end; ij++){
        int j = ij/l_lookAhead_block_num + u_start;
        int i = ij%l_lookAhead_block_num;

        // Getting U block information
        int_t rukp =  Ublock_info[j].rukp;
        int_t iukp =  Ublock_info[j].iukp;
        int jb   =  Ublock_info[j].jb;
        int nsupc = SuperSize(jb);
        int ljb = LBj (jb, grid);
        int st_col = j>u_start ? Ublock_info[j-1].full_u_cols : 0;

        /* Getting L block L(i,k) information */
        int_t lptr = lookAhead_lptr[i];
        int ib   = lookAhead_ib[i];
        int temp_nbrow = lsub[lptr+1];
        lptr += LB_DESCRIPTOR;
        int cum_nrow = (i==0 ? 0 : lookAheadFullRow[i-1]);

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
            int_t* segment_ptr = ldm_malloc((segment_count+1) * sizeof(int_t));
            int_t* segment_offset = ldm_malloc(segment_count * sizeof(int_t));
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
            for (int jj = 0; jj < nsupc; ++jj) {
                int segsize = klst - usub[iukp + jj];
                int fnz = index[iuip_lib+jj];
                if (segsize) {          /* Nonzero segment in U(k,j). */
                    doublecomplex *ucol = &Unzval_br_ptr[lib][ruip_lib];
                    // for (int i = 0; i < temp_nbrow; ++i) {
                    //     int rel = lsub[lptr + i] - fnz;
                    //     z_sub(&ucol[rel], &ucol[rel], &tempv[i]);
                    // } /* for i = 0:temp_nbropw */
                    for(int ptr = 0; ptr < segment_count; ++ptr){
                        int i_start = segment_ptr[ptr];
                        int i_end = segment_ptr[ptr+1];
                        int offset = segment_offset[ptr];
                        int rel = offset - fnz;
                        double *UCOL = (double*)&ucol[rel];
                        double *TEMPV = (double*)tempv;
                        for(int i = i_start * 2; i < i_end * 2; i++){
                            // z_sub(&UCOL[i], &UCOL[i], &tempv[i]);
                            UCOL[i] -= TEMPV[i];
                        }
                    }
                    tempv += nbrow; /* Jump LDA to next column */
                }  /* if segsize */
                ruip_lib += ilst - fnz;
                // printf("ilst, fnz, ilst - fnz, ruip_lib : %d %d %d %d\n",ilst,fnz,ilst - fnz,ruip_lib);
            }  /* for jj = 0:nsupc */

            ldm_free(segment_offset,(segment_count+1) * sizeof(int_t));
            ldm_free(segment_ptr,(segment_count+1) * sizeof(int_t));

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
                // observe
                // printf("i -> indirect2[i] indirect2[i]-i : %d, %d, %d\n",i ,indirect2[i],indirect2[i]-i);
                // assert(i == indirect2[i]), "Assert i == indirect2[i] in " ,__FILE__, ":",__LINE__;
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
            int_t* segment_ptr = ldm_malloc((segment_count+1) * sizeof(int_t));
            int_t* segment_offset = ldm_malloc(segment_count * sizeof(int_t));
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
            for (int jj = 0; jj < nsupc; ++jj) {
                int_t segsize = klst - usub[iukp + jj];
                if (segsize) {
                    // for (int i = 0; i < temp_nbrow; ++i) {
                    //     z_sub(&nzval[indirect2[i]], &nzval[indirect2[i]], &tempv[i]);
                    // }
                    for(int ptr = 0; ptr < segment_count; ++ptr){
                        int i_start = segment_ptr[ptr];
                        int i_end = segment_ptr[ptr+1];
                        int offset = segment_offset[ptr];
                        double *NZVAL = (double*)(nzval + offset);
                        double *TEMPV = (double*)(tempv);
                        for(int i = i_start*2; i < i_end*2; i++){
                            // z_sub(&NZVAL[i], &NZVAL[i], &tempv[i]);
                            NZVAL[i] -= TEMPV[i];
                        } 
                    }
                    tempv += nbrow;
                }
                nzval += ldv;
            }
            ldm_free(segment_offset,(segment_count+1) * sizeof(int_t));
            ldm_free(segment_ptr,segment_count * sizeof(int_t));
        }
    }
    ldm_free(indirect,ldt*sizeof(int));
    ldm_free(indirect2,ldt*sizeof(int));

}