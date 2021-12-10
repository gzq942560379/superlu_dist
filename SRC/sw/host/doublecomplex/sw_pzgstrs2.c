/*! \file
Copyright (c) 2003, The Regents of the University of California, through
Lawrence Berkeley National Laboratory (subject to receipt of any required
approvals from U.S. Dept. of Energy)

All rights reserved.

The source code is distributed under BSD license, see the file License.txt
at the top-level directory.
*/

/*! @file
 * \brief Performs panel LU factorization.
 *
 * <pre>
 * -- Distributed SuperLU routine (version 7.0) --
 * Lawrence Berkeley National Lab, Univ. of California Berkeley.
 * August 15, 2014
 *
 * Modified:
 *   September 30, 2017
 *   May 10, 2019 version 7.0.0
 *
 * <pre>
 * Purpose
 * =======
 *   Panel factorization -- block column k
 *
 *   Factor diagonal and subdiagonal blocks and test for exact singularity.
 *   Only the column processes that own block column *k* participate
 *   in the work.
 *
 * Arguments
 * =========
 * options (input) superlu_dist_options_t* (global)
 *         The structure defines the input parameters to control
 *         how the LU decomposition will be performed.
 *
 * k0     (input) int (global)
 *        Counter of the next supernode to be factorized.
 *
 * k      (input) int (global)
 *        The column number of the block column to be factorized.
 *
 * thresh (input) double (global)
 *        The threshold value = s_eps * anorm.
 *
 * Glu_persist (input) Glu_persist_t*
 *        Global data structures (xsup, supno) replicated on all processes.
 *
 * grid   (input) gridinfo_t*
 *        The 2D process mesh.
 *
 * Llu    (input/output) zLocalLU_t*
 *        Local data structures to store distributed L and U matrices.
 *
 * U_diag_blk_send_req (input/output) MPI_Request*
 *        List of send requests to send down the diagonal block of U.
 *
 * tag_ub (input) int
 *        Upper bound of MPI tag values.
 *
 * stat   (output) SuperLUStat_t*
 *        Record the statistics about the factorization.
 *        See SuperLUStat_t structure defined in util.h.
 *
 * info   (output) int*
 *        = 0: successful exit
 *        < 0: if info = -i, the i-th argument had an illegal value
 *        > 0: if info = i, U(i,i) is exactly zero. The factorization has
 *             been completed, but the factor U is exactly singular,
 *             and division by zero will occur if it is used to solve a
 *             system of equations.
 * </pre>
 */

#include <math.h>
#include <stdbool.h>
#include "superlu_zdefs.h"

/*****************************************************************************
 * The following pzgstrf2_trsm is in version 6 and earlier.
 *****************************************************************************/
/*! \brief
 *
 * <pre>
 * Purpose
 * =======
 *   Panel factorization -- block column k
 *
 *   Factor diagonal and subdiagonal blocks and test for exact singularity.
 *   Only the column processes that own block column *k* participate
 *   in the work.
 *
 * Arguments
 * =========
 * options (input) superlu_dist_options_t* (global)
 *         The structure defines the input parameters to control
 *         how the LU decomposition will be performed.
 *
 * k0     (input) int (global)
 *        Counter of the next supernode to be factorized.
 *
 * k      (input) int (global)
 *        The column number of the block column to be factorized.
 *
 * thresh (input) double (global)
 *        The threshold value = s_eps * anorm.
 *
 * Glu_persist (input) Glu_persist_t*
 *        Global data structures (xsup, supno) replicated on all processes.
 *
 * grid   (input) gridinfo_t*
 *        The 2D process mesh.
 *
 * Llu    (input/output) zLocalLU_t*
 *        Local data structures to store distributed L and U matrices.
 *
 * U_diag_blk_send_req (input/output) MPI_Request*
 *        List of send requests to send down the diagonal block of U.
 *
 * tag_ub (input) int
 *        Upper bound of MPI tag values.
 *
 * stat   (output) SuperLUStat_t*
 *        Record the statistics about the factorization.
 *        See SuperLUStat_t structure defined in util.h.
 *
 * info   (output) int*
 *        = 0: successful exit
 *        < 0: if info = -i, the i-th argument had an illegal value
 *        > 0: if info = i, U(i,i) is exactly zero. The factorization has
 *             been completed, but the factor U is exactly singular,
 *             and division by zero will occur if it is used to solve a
 *             system of equations.
 * </pre>
 */


/*****************************************************************************
 * The following pdgstrf2_sw is improved for SW by guozhuoqiang
 *****************************************************************************/
void pzgstrs2_sw
(int_t k0, int_t k, Glu_persist_t * Glu_persist, gridinfo_t * grid,
 zLocalLU_t * Llu, Ublock_info_t *Ublock_info, SuperLUStat_t * stat)
{
//  assert define USE_VENDOR_BLAS
#ifndef USE_VENDOR_BLAS
    printf("pzgstrs2_omp only support USE_VENDOR_BLAS !!!");
    exit(-1);
#endif 
    int iam, pkk;
    int incx = 1;
    int nsupr;                /* number of rows in the block L(:,k) (LDA) */
    int segsize;
    int nsupc;                /* number of columns in the block */
    int_t luptr, iukp, rukp;
    int_t b, gb, j, klst, knsupc, lk, nb;
    int_t *xsup = Glu_persist->xsup;
    int_t *usub;
    doublecomplex *lusup, *uval;

    /* Quick return. */
    lk = LBi (k, grid);         /* Local block number */
    if (!Llu->Unzval_br_ptr[lk]) return;

    /* Initialization. */
    iam = grid->iam;
    pkk = PNUM (PROW (k, grid), PCOL (k, grid), grid);
    //int k_row_cycle = k / grid->nprow;  /* for which cycle k exist (to assign rowwise thread blocking) */
    //int gb_col_cycle;  /* cycle through block columns  */
    klst = FstBlockC (k + 1);
    knsupc = SuperSize (k);
    usub = Llu->Ufstnz_br_ptr[lk];  /* index[] of block row U(k,:) */
    uval = Llu->Unzval_br_ptr[lk];
    if (iam == pkk) {
        lk = LBj (k, grid);
        nsupr = Llu->Lrowind_bc_ptr[lk][1]; /* LDA of lusup[] */
        lusup = Llu->Lnzval_bc_ptr[lk];
    } else {
        nsupr = Llu->Lsub_buf_2[k0 % (1 + stat->num_look_aheads)][1];   /* LDA of lusup[] */
        lusup = Llu->Lval_buf_2[k0 % (1 + stat->num_look_aheads)];
    }
    /* Master thread: set up pointers to each block in the row */
    nb = usub[0];
    iukp = BR_HEADER;
    rukp = 0;
    int* blocks_index_pointers = SUPERLU_MALLOC (3 * nb * sizeof(int));
    int* blocks_value_pointers = blocks_index_pointers + nb;
    int* nsupc_temp = blocks_value_pointers + nb;
    for (b = 0; b < nb; b++) { /* set up pointers to each block */
        blocks_index_pointers[b] = iukp + UB_DESCRIPTOR;
        blocks_value_pointers[b] = rukp;
        gb = usub[iukp];
        rukp += usub[iukp+1];
        nsupc = SuperSize( gb );
        nsupc_temp[b] = nsupc;
        iukp += (UB_DESCRIPTOR + nsupc);  /* move to the next block */
    }
    /* Loop through all the blocks in the row. */
    // count trsv tasks
    int task_count = 0;
    for (b = 0; b < nb; ++b) {
        iukp = blocks_index_pointers[b];
        rukp = blocks_value_pointers[b];
        /* Loop through all the segments in the block. */
        for (j = 0; j < nsupc_temp[b]; j++) {
            segsize = klst - usub[iukp + j];
            if (segsize) {
                task_count++;
		        rukp += segsize;
            } /* end if segsize > 0 */
        } /* end for j in parallel ... */
    } /* end for b ... */

    // alloc task memory and record
    int task_index = 0;
    int* A_indexes = SUPERLU_MALLOC(task_count * sizeof(int));
    int* segsizes = SUPERLU_MALLOC(task_count * sizeof(int));
    int* nsuprs = SUPERLU_MALLOC(task_count * sizeof(int));
    int* rukps = SUPERLU_MALLOC(task_count * sizeof(int));
    for (b = 0; b < nb; ++b) {
        iukp = blocks_index_pointers[b];
        rukp = blocks_value_pointers[b];
        /* Loop through all the segments in the block. */
        for (j = 0; j < nsupc_temp[b]; j++) {
            segsize = klst - usub[iukp + j];
            if (segsize) {
                int_t luptr = (knsupc - segsize) * (nsupr + 1);
                A_indexes[task_index] = luptr;
                segsizes[task_index] = segsize;
                nsuprs[task_index] = nsupr;
                rukps[task_index] = rukp;
		        rukp += segsize;
                task_index += 1;
            } /* end if segsize > 0 */
        } /* end for j in parallel ... */
    } /* end for b ... */

    assert(task_index == task_count);

    // check task
    assert(rukps[0] == 0);
    for(task_index = 1; task_index < task_count; ++task_index){
        assert(A_indexes[task_index] == 0);
        assert(segsizes[task_index] == segsizes[task_index - 1]);
        assert(nsuprs[task_index] == nsuprs[task_index - 1]);
        assert(rukps[task_index] == (rukps[task_index - 1] + segsizes[0]));
    }
//     // use trsv
//     for(task_index = 0; task_index < task_count; ++task_index){
// #ifdef USE_VENDOR_BLAS
//         ztrsv_ ("L", "N", "U", &segsizes[task_index], &lusup[A_indexes[task_index]], &nsuprs[task_index],
//                 &uval[rukps[task_index]], &incx, 1, 1, 1);
// #else
//         ztrsv_ ("L", "N", "U", &segsizes[task_index], &lusup[A_indexes[task_index]], &nsuprs[task_index],
//                 &uval[rukps[task_index]], &incx);
// #endif
//     }

    // AX = B, X = A^{-1}B
    int m = segsizes[0];
    int n = task_count;
    doublecomplex alpha = {1.0, 0.0};
    doublecomplex* A = lusup; // m * m matrix 
    int lda = nsuprs[0];
    doublecomplex* BX = uval;  // matrix B X , m * n   
    int ldb = segsizes[0];
    // use trsm
#ifdef USE_VENDOR_BLAS
    ztrsm_("L", "L", "N", "U", &m, &n, &alpha, A, &lda, BX, &ldb, 1, 1, 1, 1);
#else
    ztrsm_("L", "L", "N", "U", &m, &n, &alpha, A, &lda, BX, &ldb);
#endif

    SUPERLU_FREE(A_indexes);
    SUPERLU_FREE(segsizes);
    SUPERLU_FREE(nsuprs);
    SUPERLU_FREE(rukps);
} /* pzgstrs2_omp */

