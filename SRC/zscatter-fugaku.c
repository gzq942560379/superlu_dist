/*! \file
Copyright (c) 2003, The Regents of the University of California, through
Lawrence Berkeley National Laboratory (subject to receipt of any required
approvals from U.S. Dept. of Energy)

All rights reserved.

The source code is distributed under BSD license, see the file License.txt
at the top-level directory.
*/

/*! @file
 * \brief Scatter the computed blocks into LU destination.
 *
 * <pre>
 * -- Distributed SuperLU routine (version 6.1.1) --
 * Lawrence Berkeley National Lab, Univ. of California Berkeley.
 * October 1, 2014
 *
 * Modified:
 *   September 18, 2017, enable SIMD vectorized scatter operation.
 *
 */
#include <math.h>
#include "superlu_zdefs.h"
#include "utils.h"

void
zscatter_l_fugaku (
            int ib,         /* row block number of source block L(i,k) */
            int ljb,        /* local column block number of dest. block L(i,j) */
            int nsupc,      /* number of columns in destination supernode */
            int_t iukp,     /* point to destination supernode's index[] */
            int_t* xsup,
            int klst,
            int nbrow,      /* LDA of the block in tempv[] */
            int_t lptr,     /* Input, point to index[] location of block L(i,k) */
            int temp_nbrow, /* number of rows of source block L(i,k) */
            int_t* usub,
            int_t* lsub,
            doublecomplex *tempv,
            int* indirect_thread,int* indirect2,
            int_t ** Lrowind_bc_ptr, doublecomplex **Lnzval_bc_ptr,
            gridinfo_t * grid)
{
    int max_supernode_size = sp_ienv_dist(3);
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

    assert(dest_nbrow >= temp_nbrow);
    // printf("(dest_nbrow ,temp_nbrow) : %d,%d\n",dest_nbrow,temp_nbrow);

    for (int i = 0; i < dest_nbrow; ++i) {
        int rel = index[lptrj + i] - fnz;
        indirect_thread[rel] = i;
    }
    /* can be precalculated? */
    for (int i = 0; i < temp_nbrow; ++i) { /* Source index is a subset of dest. */
        int rel = lsub[lptr + i] - fnz;
        indirect2[i] = indirect_thread[rel];
        // observe
        // printf("i -> indirect2[i] indirect2[i]-i : %d, %d, %d\n",i ,indirect2[i],indirect2[i]-i);
        // assert(i == indirect2[i]), "Assert i == indirect2[i] in " ,__FILE__, ":",__LINE__;
    }

    // compress indirect2 to segment
    // int segment_count = 1;
    // int segment_prev = indirect2[0];
    // for(int i = 1; i < temp_nbrow; i++){
    //     int segment_cur = indirect2[i] - i;
    //     if(segment_prev != segment_cur){
    //         segment_count += 1;
    //         segment_prev = segment_cur;
    //     }
    // }
    // int_t* segment_ptr = SUPERLU_MALLOC((max_supernode_size+1) * sizeof(int_t));
    // int_t* segment_offset = SUPERLU_MALLOC(max_supernode_size * sizeof(int_t));
    // segment_ptr[0] = 0;
    // segment_offset[0] = indirect2[0];
    // int segment_prev = indirect2[0]; 
    // int segment_count = 1;
    // for(int i = 1; i < temp_nbrow; i++){
    //     int segment_cur = indirect2[i] - i;
    //     if(segment_prev != segment_cur){
    //         segment_ptr[segment_count] = i;
    //         segment_offset[segment_count] = segment_cur;
    //         segment_count += 1;
    //         segment_prev = segment_cur;
    //     }
    // }
    // segment_ptr[segment_count] = temp_nbrow;

    indirect_index_segment_compress_t segment_compress;
    indirect_index_segment_compress_init(&segment_compress, indirect2, temp_nbrow);

    // assert(segment_index == segment_count);

    // check 
    // printf("segment_count temp_nbrow : %d,%d\n",segment_count,temp_nbrow);
    // for(int ptr = 0; ptr < segment_count; ++ptr){
    //     int i_start = segment_ptr[ptr];
    //     int i_end = segment_ptr[ptr+1];
    //     int offset = segment_offset[ptr];
    //     for(int i = i_start; i < i_end; i++){
    //         int true_offset = indirect2[i] - i;
    //         printf("i : true_offset, offset %d %d %d\n",i, true_offset, offset);
    //         assert(true_offset == offset), "Assert true_offset == offset in ", __FILE__,":",__LINE__; 
    //     }
    // }

    doublecomplex *nzval = Lnzval_bc_ptr[ljb] + luptrj; /* Destination block L(i,j) */

    // for (int jj = 0; jj < nsupc; ++jj) {
    //     int_t segsize = klst - usub[iukp + jj];
    //     if (segsize) {
    //         // for (int i = 0; i < temp_nbrow; ++i) {
    //         //     z_sub(&nzval[indirect2[i]], &nzval[indirect2[i]], &tempv[i]);
    //         // }
    //         for(int ptr = 0; ptr < segment_count; ++ptr){
    //             int i_start = segment_ptr[ptr];
    //             int i_end = segment_ptr[ptr+1];
    //             int offset = segment_offset[ptr];
    //             doublecomplex *NZVAL = nzval + offset;
    //             for(int i = i_start; i < i_end; i++){
    //                 z_sub(&NZVAL[i], &NZVAL[i], &tempv[i]);
    //             } 
    //         }
    //         tempv += nbrow;
    //     }
    //     nzval += ldv;
    // }
    // printf("segment_count : %d\n",segment_compress.segment_count);

    int non_zero_col_num = 0;
    doublecomplex* nzvals[max_supernode_size];
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
        for(int ptr = 0; ptr < segment_compress.segment_count; ++ptr){
            int i_start = segment_compress.segment_ptr[ptr];
            int i_end = segment_compress.segment_ptr[ptr+1];
            int offset = segment_compress.segment_offset[ptr];
            double *NZVAL = (double*)(nzval_cur + offset);
            // for(int i = i_start; i < i_end; i++){
            //     z_sub(&NZVAL[i], &NZVAL[i], &tempv_cur[i]);
            // } 
            #pragma omp simd
            for(int i = i_start*2; i < i_end*2; i++){
                NZVAL[i] -= tempv_cur[i];
            } 
        }
    }

    indirect_index_segment_compress_destroy(&segment_compress);
} /* zscatter_l */


void
zscatter_u_fugaku ( 
            int ib,
            int jb,
            int nsupc,
            int_t iukp,
            int_t * xsup,
            int klst,
            int nbrow,      /* LDA of the block in tempv[] */
            int_t lptr,     /* point to index location of block L(i,k) */
            int temp_nbrow, /* number of rows of source block L(i,k) */
            int_t* lsub,
            int_t* usub,
            doublecomplex* tempv,
            int_t ** Ufstnz_br_ptr, doublecomplex **Unzval_br_ptr,
            gridinfo_t * grid)
{
    int max_supernode_size = sp_ienv_dist(3);
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
    /* Skip descriptor. Now point to fstnz index of block U(i,j). */
    iuip_lib += UB_DESCRIPTOR;

    assert(temp_nbrow > 0);

    // observe
    // for(int i = 0;i < temp_nbrow;i++){
    //     printf("i lsub[lptr + i] lsub[lptr + i]-i : %d %d %d\n",i,lsub[lptr + i],lsub[lptr + i]-i);
    // }

    // offset compress to segment
    // int segment_count = 1;
    // int segment_prev = lsub[lptr];
    // for(int i = 1; i < temp_nbrow; i++){
    //     int segment_cur = lsub[lptr + i] - i;
    //     if(segment_prev != segment_cur){
    //         segment_count += 1;
    //         segment_prev = segment_cur;
    //     }
    // }
    // int_t* segment_ptr = SUPERLU_MALLOC((max_supernode_size+1) * sizeof(int_t));
    // int_t* segment_offset = SUPERLU_MALLOC(max_supernode_size * sizeof(int_t));
    // segment_ptr[0] = 0;
    // segment_offset[0] = lsub[lptr];
    // int segment_prev = lsub[lptr];
    // int segment_count = 1;
    // for(int i = 1; i < temp_nbrow; i++){
    //     int segment_cur = lsub[lptr + i] - i;
    //     if(segment_prev != segment_cur){
    //         segment_ptr[segment_count] = i;
    //         segment_offset[segment_count] = segment_cur;
    //         segment_count += 1;
    //         segment_prev = segment_cur;
    //     }
    // }
    // segment_ptr[segment_count] = temp_nbrow;
    // assert(segment_count == segment_count);

    indirect_index_segment_compress_t segment_compress;
    indirect_index_segment_compress_init(&segment_compress, &lsub[lptr], temp_nbrow);
    
    // if(segment_count > 3)
    // printf("segment_count temp_nbrow : %d %d\n",segment_count,temp_nbrow);
    // assert(segment_count < 4), "Assert segment_count < 4 in ", __FILE__, ":", __LINE__;

    // check 
    // printf("segment_count : %d\n",segment_count);
    // for(int ptr = 0; ptr < segment_count; ++ptr){
    //     int i_start = segment_ptr[ptr];
    //     int i_end = segment_ptr[ptr+1];
    //     int offset = segment_offset[ptr];
    //     for(int i = i_start; i < i_end; i++){
    //         int true_offset = lsub[lptr + i] - i;
    //         printf("i : true_offset, offset %d %d %d\n",i, true_offset, offset);
    //         assert(true_offset == offset), "Assert true_offset == offset in ", __FILE__,":",__LINE__; 
    //     }
    // }

    // for (int jj = 0; jj < nsupc; ++jj) {
    //     int segsize = klst - usub[iukp + jj];
    //     int fnz = index[iuip_lib+jj];
    //     if (segsize) {          /* Nonzero segment in U(k,j). */
    //         doublecomplex *ucol = &Unzval_br_ptr[lib][ruip_lib];
    //         // for (int i = 0; i < temp_nbrow; ++i) {
    //         //     int rel = lsub[lptr + i] - fnz;
    //         //     z_sub(&ucol[rel], &ucol[rel], &tempv[i]);
    //         // } /* for i = 0:temp_nbropw */
    //         for(int ptr = 0; ptr < segment_count; ++ptr){
    //             int i_start = segment_ptr[ptr];
    //             int i_end = segment_ptr[ptr+1];
    //             int offset = segment_offset[ptr];
    //             int rel = offset - fnz;
    //             doublecomplex *UCOL = &ucol[rel];
    //             for(int i = i_start; i < i_end; i++){
    //                 z_sub(&UCOL[i], &UCOL[i], &tempv[i]);
    //             }
    //         }
    //         tempv += nbrow; /* Jump LDA to next column */
    //     }  /* if segsize */
    //     ruip_lib += ilst - fnz;
    //     // printf("ilst, fnz, ilst - fnz, ruip_lib : %d %d %d %d\n",ilst,fnz,ilst - fnz,ruip_lib);
    // }  /* for jj = 0:nsupc */
    // printf("segment_count : %d\n",segment_compress.segment_count);

    int non_zero_col_num = 0;
    doublecomplex * non_zero_ucols[max_supernode_size];
    int fnzs[max_supernode_size];
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
        for(int ptr = 0; ptr < segment_compress.segment_count; ++ptr){
            int i_start = segment_compress.segment_ptr[ptr];
            int i_end = segment_compress.segment_ptr[ptr+1];
            int offset = segment_compress.segment_offset[ptr];
            int rel = offset - fnz;
            double *UCOL = (double*)&ucol_cur[rel];
            // for(int i = i_start; i < i_end; i++){
            //     z_sub(&UCOL[i], &UCOL[i], &tempv_cur[i]);
            // }
            #pragma omp simd
            for(int i = i_start*2; i < i_end*2; i++){
                UCOL[i] -= tempv_cur[i];
            }
        }
    }

    if(segment_compress.segment_count == 1){
        for (int jj_ptr = 0; jj_ptr < non_zero_col_num; ++jj_ptr) {
        int fnz = fnzs[jj_ptr];
        doublecomplex *ucol_cur = non_zero_ucols[jj_ptr];
        double *tempv_cur = (double*)(tempv + nbrow * jj_ptr);
            for(int ptr = 0; ptr < segment_compress.segment_count; ++ptr){
                int i_start = segment_compress.segment_ptr[ptr];
                int i_end = segment_compress.segment_ptr[ptr+1];
                int offset = segment_compress.segment_offset[ptr];
                int rel = offset - fnz;
                double *UCOL = (double*)&ucol_cur[rel];
                // for(int i = i_start; i < i_end; i++){
                //     z_sub(&UCOL[i], &UCOL[i], &tempv_cur[i]);
                // }
                #pragma omp simd
                for(int i = i_start*2; i < i_end*2; i++){
                    UCOL[i] -= tempv_cur[i];
                }
            }
        }
    }

    indirect_index_segment_compress_destroy(&segment_compress);
} /* zscatter_u */
