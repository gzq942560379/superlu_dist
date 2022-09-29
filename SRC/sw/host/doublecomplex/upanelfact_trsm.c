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
 * -- Distributed SuperLU routine (version 7.2) --
 * Lawrence Berkeley National Lab, Univ. of California Berkeley.
 * August 15, 2014
 *
 * Modified:
 *   September 30, 2017
 *   May 10, 2019  v7.0.0
 *   December 12, 2021  v7.2.0
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
#include "sw/superlu_zdefs_ex.h"

void upanelfact_trsm(int_t k0, int_t k, Glu_persist_t * Glu_persist, gridinfo_t* grid, zLocalLU_t* Llu, SuperLUStat_t* stat)
{
	/* Initialization. */
    int iam = grid->iam;
	int_t krow = PROW( k, grid );
	int_t kcol = PCOL( k, grid );	
	int_t myrow = MYROW( iam, grid );
    int_t mycol = MYCOL( iam, grid );
	if(myrow == krow){
		/* Quick return. */
		int_t lki = LBi(k, grid ); /* Local block number */
		if ( !Llu->Unzval_br_ptr[lki] ) return;

    	int_t *xsup = Glu_persist->xsup;
		int_t klst = FstBlockC( k+1 );
		int_t knsupc= SuperSize( k );
		int_t* usub = Llu->Ufstnz_br_ptr[lki]; /* index[] of block row U(k,:) */
		doublecomplex* uval = Llu->Unzval_br_ptr[lki];

        int_t ldl;
        doublecomplex* lval;

        if (mycol == kcol) {
            int_t lkj = LBj (k, grid);
            ldl = Llu->Lrowind_bc_ptr[lkj][1]; /* LDA of lusup[] */
            lval = Llu->Lnzval_bc_ptr[lkj];
        } else {
            ldl = Llu->Lsub_buf_2[k0 % (1 + stat->num_look_aheads)][1];   /* LDA of lusup[] */
            lval = Llu->Lval_buf_2[k0 % (1 + stat->num_look_aheads)];
        }

		int_t index = BR_HEADER;
		int_t nub = usub[0];
		int ucols = 0;
		for (int_t b = 0; b < nub; ++b) {
			int_t gbc = usub[index];
			int_t nsupc = SuperSize(gbc);
            int_t block_nzval_len = usub[index+1];
            int_t block_ucols = 0;
			index += UB_DESCRIPTOR;
			/* Loop through all the segments in the block. */
			for (int_t j = 0; j < nsupc; ++j) {
				int_t segsize = klst - usub[index++];
				if ( segsize ) { /* Nonzero segment. */
                    assert(segsize == knsupc);
                    block_ucols += 1;
				}
			}
            assert(block_nzval_len == block_ucols * knsupc);
            ucols += block_ucols;
		} /* for b ... */
		zMatrix_t mLkk, mAUk_;
		zMatrix_init(&mLkk, knsupc, knsupc, ldl, lval);
		zMatrix_init(&mAUk_, knsupc, ucols, knsupc, uval);
		doublecomplex ONE = {1.0, 0.0};
        
    	ztrsm_("L", "L", "N", "U", &mLkk.row, &mAUk_.col, &ONE, mLkk.val, &mLkk.ld, mAUk_.val, &mAUk_.ld, 1, 1, 1, 1);

		stat->ops[FACT] += 4 * mLkk.row * (mLkk.row + 1) * mAUk_.col;
	}
} /* upanelfact_trsm */
