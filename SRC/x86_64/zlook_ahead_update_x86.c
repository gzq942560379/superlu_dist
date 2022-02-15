/*! \file
Copyright (c) 2003, The Regents of the University of California, through
Lawrence Berkeley National Laboratory (subject to receipt of any required
approvals from U.S. Dept. of Energy)

All rights reserved.

The source code is distributed under BSD license, see the file License.txt
at the top-level directory.
*/

/************************************************************************/
/*! @file
 * \brief Look-ahead update of the Schur complement.
 *
 * <pre>
 * -- Distributed SuperLU routine (version 5.4) --
 * Lawrence Berkeley National Lab, Univ. of California Berkeley.
 * October 1, 2014
 *
 * Modified:
 *  September 18, 2017
 *  June 1, 2018  add parallel AWPM pivoting; add back arrive_at_ublock()
 *
 */

#include <assert.h>  /* assertion doesn't work if NDEBUG is defined */

iukp = iukp0; /* point to the first block in index[] */
rukp = rukp0; /* point to the start of nzval[] */
j = jj0 = 0;  /* After the j-loop, jj0 points to the first block in U
                 outside look-ahead window. */

#ifndef USE_VENDOR_BLAS
assert(false);
#endif

#ifndef ISORT
assert(false);
#endif

while (j < nub && iperm_u[j] <= k0 + num_look_aheads)
{
#if defined(USE_X86) && defined(OPT_gather_avoid)
    zMatrix_t U;
    zMatrix_t L_all;
    zMatrix_t *L_list = SUPERLU_MALLOC(nlb * sizeof(zMatrix_t));
#endif
    doublecomplex zero = {0.0, 0.0};

    /* Search is needed because a permutation perm_u is involved for j  */
    /* Search along the row for the pointers {iukp, rukp} pointing to
     * block U(k,j).
     * j    -- current block in look-ahead window, initialized to 0 on entry
     * iukp -- point to the start of index[] metadata
     * rukp -- point to the start of nzval[] array
     * jb   -- block number of block U(k,j), update destination column
     */
    arrive_at_ublock(
		     j, &iukp, &rukp, &jb, &ljb, &nsupc,
         	    iukp0, rukp0, usub, perm_u, xsup, grid
		    );

    j++;
    jj0++;
    jj = iukp;

    while (usub[jj] == klst) ++jj; /* Skip zero segments */

    ldu = klst - usub[jj++];
    ncols = 1;

    /* This loop computes ldu. */
    for (; jj < iukp + nsupc; ++jj) { /* for each column jj in block U(k,j) */
        segsize = klst - usub[jj];
        if (segsize) {
            ++ncols;
            if (segsize > ldu)  ldu = segsize;
        }
    }
#if ( DEBUGlevel>=3 )
    ++num_update;
#endif

#if ( DEBUGlevel>=3 )
    printf ("(%d) k=%d,jb=%d,ldu=%d,ncols=%d,nsupc=%d\n",
	    iam, k, jb, ldu, ncols, nsupc);
    ++num_copy;
#endif

#if defined(USE_X86) && defined(OPT_gather_avoid)
    zMatrix_init(&U,ldu,ncols,ldu,&uval[rukp]);
#else
    /* Now copy one block U(k,j) to bigU for GEMM, padding zeros up to ldu. */
    tempu = bigU; /* Copy one block U(k,j) to bigU for GEMM */
    for (jj = iukp; jj < iukp + nsupc; ++jj) {
        segsize = klst - usub[jj];
        if (segsize) {
            lead_zero = ldu - segsize;
            for (i = 0; i < lead_zero; ++i) tempu[i] = zero;
            tempu += lead_zero;
            for (i = 0; i < segsize; ++i) {
                tempu[i] = uval[rukp + i];
            }
            rukp += segsize;
            tempu += segsize;
        }
    }
#endif
    tempu = bigU; /* set back to the beginning of the buffer */

    nbrow = lsub[1]; /* number of row subscripts in L(:,k) */
    if (myrow == krow) nbrow = lsub[1] - lsub[3]; /* skip diagonal block for those rows. */
    // double ttx =SuperLU_timer_();

    lptr = lptr0; /* point to the start of index[] in supernode L(:,k) */
    luptr = luptr0;

    int all_row = 0;
    LBlock_info_t* LBlock_info = SUPERLU_MALLOC(nlb * sizeof(LBlock_info_t));


    for (lb = 0; lb < nlb; lb++) { /* Loop through each block in L(:,k) */

        ib = lsub[lptr];        /* block number of L(i,k) */
        int temp_nbrow = lsub[lptr + 1];    /* Number of full rows. */

#if defined(USE_X86) && defined(OPT_gather_avoid)
        zMatrix_init(&L_list[lb], temp_nbrow, ldu, nsupr, &lusup[luptr + (knsupc - ldu) * nsupr]);
#endif

        LBlock_info[lb].ib = ib;
        LBlock_info[lb].lptr = lptr;
        LBlock_info[lb].FullRow = ((lb == 0) ?  temp_nbrow : (LBlock_info[lb-1].FullRow + temp_nbrow));
        
	    all_row += temp_nbrow;
        lptr += LB_DESCRIPTOR;  /* Skip descriptor. */
        lptr += temp_nbrow;
        luptr += temp_nbrow;

    } /* end parallel for lb = 0, nlb ... all blocks in L(:,k) */
#if defined(USE_X86) && defined(OPT_gather_avoid)
    zMatrix_init(&L_all, all_row, ldu, nsupr, &lusup[luptr0 + (knsupc - ldu) * nsupr]);
#endif

    if(all_row > 0 && ldu > 0 && ncols > 0){
        /* calling gemm */
        double flops = 8.0 * (flops_t)all_row * ldu * ncols;
        stat->ops[FACT] += flops;
        lookaheadupdateflops += flops;

#if defined(USE_X86) && defined(OPT_gather_avoid)
        zMatrix_t* L = &L_list[lb];
        assert(L_all.col == U.row);
        zgemm_("N", "N", &L_all.row, &U.col, &L_all.col, &alpha,
            L_all.val, &L_all.ld,
            U.val, &U.ld, &beta, bigV, &all_row, 1, 1);        
#else
        zgemm_("N", "N", &all_row, &ncols, &ldu, &alpha,
            &lusup[luptr0 + (knsupc - ldu) * nsupr], &nsupr,
            tempu, &ldu, &beta, bigV, &all_row, 1, 1);
#endif

#if defined(USE_X86) && defined(parallel_scatter) && defined(_OPENMP)
#pragma omp parallel default(shared)
#endif
        {

#if defined(USE_X86) && defined(parallel_scatter) && defined(_OPENMP)
            int thread_id = omp_get_thread_num();
            int* indirect_thread = indirect + (ldt + CACHELINE/sizeof(int)) * thread_id;
            int* indirect2_thread = indirect2 + (ldt + CACHELINE/sizeof(int)) * thread_id;
#pragma omp for private (lb,lptr,ib) schedule(dynamic)
#else /* not use _OPENMP */
            int thread_id = 0;
            int* indirect_thread = indirect;
            int* indirect2_thread = indirect2;
#endif   
            for (lb = 0; lb < nlb; lb++) { /* Loop through each block in L(:,k) */
                ib = LBlock_info[lb].ib;
                // ib = LBlock_info[lb].ib;        /* block number of L(i,k) */
                lptr = LBlock_info[lb].lptr;
                int temp_nbrow = lsub[lptr + 1];    /* Number of full rows. */
                lptr += LB_DESCRIPTOR;  /* Skip descriptor. */
                int cum_nrow = (lb==0 ? 0 : LBlock_info[lb-1].FullRow);
                
                doublecomplex * tempv = bigV + cum_nrow;

                /* Now scattering the output. */
                if (ib < jb) {    /* A(i,j) is in U. */
#if defined(USE_X86) && defined(OPT_scatter_index_compress)
                    zscatter_u_x86 (ib, jb,
                            nsupc, iukp, xsup,
                            klst, all_row,
                            lptr, temp_nbrow, lsub,
                            usub, tempv, Ufstnz_br_ptr, Unzval_br_ptr, grid);
#else
                    zscatter_u (ib, jb,
                            nsupc, iukp, xsup,
                            klst, all_row,
                            lptr, temp_nbrow, lsub,
                            usub, tempv, Ufstnz_br_ptr, Unzval_br_ptr, grid);
#endif
                } else {          /* A(i,j) is in L. */
#if defined(USE_X86) && defined(OPT_scatter_index_compress)
                    zscatter_l_x86 (ib, ljb, nsupc, iukp, xsup, klst, all_row, lptr,
                            temp_nbrow, usub, lsub, tempv,
                            indirect_thread, indirect2_thread,
                            Lrowind_bc_ptr, Lnzval_bc_ptr, grid);
#else
                    zscatter_l (ib, ljb, nsupc, iukp, xsup, klst, all_row, lptr,
                            temp_nbrow, usub, lsub, tempv,
                            indirect_thread, indirect2_thread,
                            Lrowind_bc_ptr, Lnzval_bc_ptr, grid);
#endif
                }
            } /* end parallel for lb = 0, nlb ... all blocks in L(:,k) */
        }
    }
    SUPERLU_FREE(LBlock_info);

    iukp += nsupc; /* Mov to block U(k,j+1) */

    /* =========================================== *
     * == factorize L(:,j) and send if possible == *
     * =========================================== */
    kk = jb; /* destination column that is just updated */
    kcol = PCOL (kk, grid);
    kk0 = iperm_u[j - 1];
    look_id = kk0 % (1 + num_look_aheads);

    if (look_ahead[kk] == k0 && kcol == mycol) {
        /* current column is the last dependency */
        look_id = kk0 % (1 + num_look_aheads);

        /* Factor diagonal and subdiagonal blocks and test for exact
           singularity.  */
        factored[kk] = 0;

        double tt1 = SuperLU_timer_();

        PZGSTRF2(options, kk0, kk, thresh, Glu_persist, grid, Llu,
                  U_diag_blk_send_req, tag_ub, stat, info);

        pdgstrf2_timer += SuperLU_timer_() - tt1;

        /* stat->time7 += SuperLU_timer_() - ttt1; */

        /* Multicasts numeric values of L(:,kk) to process rows. */
        send_req = send_reqs[look_id];
        msgcnt = msgcnts[look_id];

        lk = LBj (kk, grid);    /* Local block number. */
        lsub1 = Lrowind_bc_ptr[lk];
        lusup1 = Lnzval_bc_ptr[lk];
        if (lsub1) {
            msgcnt[0] = lsub1[1] + BC_HEADER + lsub1[0] * LB_DESCRIPTOR;
            msgcnt[1] = lsub1[1] * SuperSize (kk);
        } else {
            msgcnt[0] = 0;
            msgcnt[1] = 0;
        }

        scp = &grid->rscp;      /* The scope of process row. */
        for (pj = 0; pj < Pc; ++pj) {
            if (ToSendR[lk][pj] != EMPTY) {

#if ( PROFlevel>=1 )
                TIC (t1);
#endif
                MPI_Isend (lsub1, msgcnt[0], mpi_int_t, pj,
                           SLU_MPI_TAG (0, kk0) /* (4*kk0)%tag_ub */ ,
                           scp->comm, &send_req[pj]);
                MPI_Isend (lusup1, msgcnt[1], SuperLU_MPI_DOUBLE_COMPLEX, pj,
                           SLU_MPI_TAG (1, kk0) /* (4*kk0+1)%tag_ub */ ,
                           scp->comm, &send_req[pj + Pc]);
#if ( PROFlevel>=1 )
                TOC (t2, t1);
                stat->utime[COMM] += t2;
                msg_cnt += 2;
                msg_vol += msgcnt[0] * iword + msgcnt[1] * dword;
#endif

#if ( DEBUGlevel>=2 )
                printf ("[%d] -2- Send L(:,%4d): #lsub %4d, #lusup %4d to Pj %2d, tags %d:%d \n",
                        iam, kk, msgcnt[0], msgcnt[1], pj,
			    SLU_MPI_TAG(0,kk0), SLU_MPI_TAG(1,kk0));
#endif
            }  /* end if ( ToSendR[lk][pj] != EMPTY ) */
        } /* end for pj ... */
    } /* end if( look_ahead[kk] == k0 && kcol == mycol ) */

#if defined(USE_X86) && defined(OPT_gather_avoid)
    SUPERLU_FREE(L_list);
#endif

} /* end while j < nub and perm_u[j] <k0+NUM_LOOK_AHEAD */
