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
#include "superlu_zdefs.h"
//#include "cblas.h"

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
/* This pdgstrf2 is based on TRSM function */
void
sw_pzgstrf2_trsm
    (superlu_dist_options_t * options, int_t k0, int_t k, double thresh,
     Glu_persist_t * Glu_persist, gridinfo_t * grid, zLocalLU_t * Llu,
     MPI_Request * U_diag_blk_send_req, int tag_ub,
     SuperLUStat_t * stat, int *info)
{
    /* printf("entering pzgstrf2 %d \n", grid->iam); */
    int cols_left, iam, l, pkk, pr;
    int incx = 1, incy = 1;

    int nsupr;            /* number of rows in the block (LDA) */
    int nsupc;            /* number of columns in the block */
    int luptr;
    int_t i, myrow, krow, j, jfst, jlst, u_diag_cnt;
    int_t *xsup = Glu_persist->xsup;
    doublecomplex *lusup, temp;
    doublecomplex *ujrow, *ublk_ptr;   /* pointer to the U block */
    doublecomplex alpha = {-1.0, 0.0}, zero = {0.0, 0.0}, one = {1.0, 0.0};
    int_t Pr;
    MPI_Status status;
    MPI_Comm comm = (grid->cscp).comm;
    double t1, t2;

    /* Initialization. */
    iam = grid->iam;
    Pr = grid->nprow;
    myrow = MYROW (iam, grid);
    krow = PROW (k, grid);
    pkk = PNUM (PROW (k, grid), PCOL (k, grid), grid);
    j = LBj (k, grid);          /* Local block number */
    jfst = FstBlockC (k);
    jlst = FstBlockC (k + 1);
    lusup = Llu->Lnzval_bc_ptr[j];
    nsupc = SuperSize (k);
    if (Llu->Lrowind_bc_ptr[j])
        nsupr = Llu->Lrowind_bc_ptr[j][1];
    else
        nsupr = 0;
    ublk_ptr = ujrow = Llu->ujrow;

    luptr = 0;                  /* Point to the diagonal entries. */
    cols_left = nsupc;          /* supernode size */
    int ld_ujrow = nsupc;       /* leading dimension of ujrow */
    u_diag_cnt = 0;
    incy = ld_ujrow;

    if ( U_diag_blk_send_req &&
	 U_diag_blk_send_req[myrow] != MPI_REQUEST_NULL ) {
        /* There are pending sends - wait for all Isend to complete */

        for (pr = 0; pr < Pr; ++pr) {
            if (pr != myrow) {
                MPI_Wait (U_diag_blk_send_req + pr, &status);
            }
	}
	/* flag no more outstanding send request. */
	U_diag_blk_send_req[myrow] = MPI_REQUEST_NULL;
    }

    if (iam == pkk) {            /* diagonal process */
	/* ++++ First step compute diagonal block ++++++++++ */
        for (j = 0; j < jlst - jfst; ++j) {  /* for each column in panel */
            /* Diagonal pivot */
            i = luptr;
           if ( options->ReplaceTinyPivot == YES ) {
                if ( slud_z_abs1(&lusup[i]) < thresh &&
		     lusup[i].r != 0.0 && lusup[i].i != 0.0 ) { /* Diagonal */
#if ( PRNTlevel>=2 )
                    printf ("(%d) .. col %d, tiny pivot %e  ",
                            iam, jfst + j, lusup[i]);
#endif
                    /* Keep the new diagonal entry with the same sign. */
                    if ( lusup[i].r < 0 ) lusup[i].r = -thresh;
                    else lusup[i].r = thresh;
                    lusup[i].i = 0.0;
                    ++(stat->TinyPivots);
                }
            }

            /* storing U in full form  */
            int st;
            for (l = 0; l < cols_left; ++l, i += nsupr, ++u_diag_cnt) {
                st = j * ld_ujrow + j;
                ublk_ptr[st + l * ld_ujrow] = lusup[i]; /* copy one row of U */
            }

            /* Test for singularity. */
            if ( ujrow[0].r == 0.0 && ujrow[0].i == 0.0 ) {
                *info = j + jfst + 1;
            } else {              /* Scale the j-th column within diag. block. */
                slud_z_div(&temp, &one, &ujrow[0]);
                for (i = luptr + 1; i < luptr - j + nsupc; ++i)
                    zz_mult(&lusup[i], &lusup[i], &temp);
                stat->ops[FACT] += 6*(nsupc-j-1) + 10;
            }

            /* Rank-1 update of the trailing submatrix within diag. block. */
            if (--cols_left) {
                /* l = nsupr - j - 1;  */
                l = nsupc - j - 1;  /* Piyush */
                zgeru_(&l, &cols_left, &alpha, &lusup[luptr+1], &incx,
                       &ujrow[ld_ujrow], &incy, &lusup[luptr + nsupr + 1],
                       &nsupr);
                stat->ops[FACT] += 8 * l * cols_left;
            }

            /* ujrow = ublk_ptr + u_diag_cnt;  */
            ujrow = ujrow + ld_ujrow + 1; /* move to next row of U */
            luptr += nsupr + 1; /* move to next column */

        }                       /* for column j ...  first loop */

	    /* ++++ Second step compute off-diagonal block with communication  ++*/

        ublk_ptr = ujrow = Llu->ujrow;

        if (U_diag_blk_send_req && iam == pkk)  { /* Send the U block downward */
            /** ALWAYS SEND TO ALL OTHERS - TO FIX **/
            for (pr = 0; pr < Pr; ++pr) {
                if (pr != krow) {
                    /* tag = ((k0<<2)+2) % tag_ub;        */
                    /* tag = (4*(nsupers+k0)+2) % tag_ub; */
                    MPI_Isend (ublk_ptr, nsupc * nsupc, SuperLU_MPI_DOUBLE_COMPLEX, pr,
                               SLU_MPI_TAG (4, k0) /* tag */ ,
                               comm, U_diag_blk_send_req + pr);

                }
            }


	    /* flag outstanding Isend */
            U_diag_blk_send_req[krow] = (MPI_Request) TRUE; /* Sherry */
        }

        /* pragma below would be changed by an MKL call */

        l = nsupr - nsupc;
        // n = nsupc;
	    doublecomplex alpha = {1.0, 0.0};
#if defined (USE_VENDOR_BLAS)
        ztrsm_ ("R", "U", "N", "N", &l, &nsupc,
                &alpha, ublk_ptr, &ld_ujrow, &lusup[nsupc], &nsupr,
		1, 1, 1, 1);
#else
        ztrsm_ ("R", "U", "N", "N", &l, &nsupc,
                &alpha, ublk_ptr, &ld_ujrow, &lusup[nsupc], &nsupr);
#endif
	    stat->ops[FACT] += 4.0 * ((flops_t) nsupc * (nsupc+1) * l);
    } else {  /* non-diagonal process */
        /* ================================================================== *
         * Receive the diagonal block of U for panel factorization of L(:,k). *
         * Note: we block for panel factorization of L(:,k), but panel        *
	     * factorization of U(:,k) do not block                               *
         * ================================================================== */

        /* tag = ((k0<<2)+2) % tag_ub;        */
        /* tag = (4*(nsupers+k0)+2) % tag_ub; */
        // printf("hello message receiving%d %d\n",(nsupc*(nsupc+1))>>1,SLU_MPI_TAG(4,k0));
        MPI_Recv (ublk_ptr, (nsupc * nsupc), SuperLU_MPI_DOUBLE_COMPLEX, krow,
                  SLU_MPI_TAG (4, k0) /* tag */ ,
                  comm, &status);
        if (nsupr > 0) {
            doublecomplex alpha = {1.0, 0.0};

#if defined (USE_VENDOR_BLAS)
            ztrsm_ ("R", "U", "N", "N", &nsupr, &nsupc,
                    &alpha, ublk_ptr, &ld_ujrow, lusup, &nsupr, 1, 1, 1, 1);
#else
            ztrsm_ ("R", "U", "N", "N", &nsupr, &nsupc,
                    &alpha, ublk_ptr, &ld_ujrow, lusup, &nsupr);
#endif
	    stat->ops[FACT] += 4.0 * ((flops_t) nsupc * (nsupc+1) * nsupr);
        }

    } /* end if pkk ... */


    /* printf("exiting pzgstrf2 %d \n", grid->iam);  */

}  /* PZGSTRF2_trsm */
