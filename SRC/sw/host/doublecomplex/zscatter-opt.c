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
#include "sw/superlu_zdefs_ex.h"
#include "sw/utils.h"

#ifndef MAX_BLOCK_SIZE
#define MAX_BLOCK_SIZE 256
#endif

void
zscatter_l_opt (
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
    int iam = grid->iam;
    
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

    int_t fnz = FstBlockC (ib);
    lptrj += LB_DESCRIPTOR;
    int_t dest_nbrow=index[lptrj - 1];

    assert(dest_nbrow >= temp_nbrow);

    for (int i = 0; i < dest_nbrow; ++i) {
        int rel = index[lptrj + i] - fnz;
        indirect_thread[rel] = i;
    }
    /* can be precalculated? */
    for (int i = 0; i < temp_nbrow; ++i) { /* Source index is a subset of dest. */
        int rel = lsub[lptr + i] - fnz;
        indirect2[i] = indirect_thread[rel];
    }

    indirect_index_segment_compress_t segment_compress;
    indirect_index_segment_compress_init(&segment_compress, indirect2, temp_nbrow);

    doublecomplex *nzval = Lnzval_bc_ptr[ljb] + luptrj; /* Destination block L(i,j) */

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
            #pragma omp simd
            for(int i = i_start*2; i < i_end*2; i++){
                NZVAL[i] -= tempv_cur[i];
            } 
        }
    }

    indirect_index_segment_compress_destroy(&segment_compress);
} /* zscatter_l */



void zscatter_u_opt ( 
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
    int iam = grid->iam;
    int_t ilst = FstBlockC (ib + 1);
    int_t lib = LBi (ib, grid);
    int_t *index = Ufstnz_br_ptr[lib];

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

    indirect_index_segment_compress_t segment_compress;
    indirect_index_segment_compress_init(&segment_compress, &lsub[lptr], temp_nbrow);
    
    int non_zero_col_num = 0;
    doublecomplex * non_zero_ucols[MAX_BLOCK_SIZE];
    int fnzs[MAX_BLOCK_SIZE];
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

    indirect_index_segment_compress_destroy(&segment_compress);
} /* zscatter_u */
