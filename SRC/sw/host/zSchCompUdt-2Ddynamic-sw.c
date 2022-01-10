/*! \file
Copyright (c) 2003, The Regents of the University of California, through
Lawrence Berkeley National Laboratory (subject to receipt of any required
approvals from U.S. Dept. of Energy)

All rights reserved.

The source code is distributed under BSD license, see the file License.txt
at the top-level directory.
*/

/*! @file
 * \brief This file contains the main loop of pdgstrf which involves rank k
 *        update of the Schur complement.
 *        Uses 2D partitioning for the scatter phase.
 *
 * <pre>
 * -- Distributed SuperLU routine (version 5.4) --
 * Lawrence Berkeley National Lab, Univ. of California Berkeley.
 * October 1, 2014
 *
 * Modified:
 *   September 14, 2017
 *   - First gather U-panel, then depending on "ldu" (excluding leading zeros),
 *     gather only trailing columns of the L-panel corresponding to the nonzero
 *     of U-rows.
 *   - Padding zeros for nice dimensions of GEMM.
 *
 *  June 1, 2018  add parallel AWPM pivoting; add back arrive_at_ublock()
 */

/*
 * Buffers:
 *     [ lookAhead_L_buff | Remain_L_buff ] : stores the gathered L-panel
 *                                            (A matrix in C := A*B )
 *     bigU : stores the U-panel (B matrix in C := A*B)
 *     bigV : stores the block GEMM result (C matrix in C := A*B)
 */

#ifndef USE_VENDOR_BLAS
	assert(false);
#endif
if ( msg0 && msg2 ) { /* L(:,k) and U(k,:) are not empty. */
    int cum_nrow = 0; /* cumulative number of nonzero rows in L(:,k) */
    int temp_nbrow;   /* nonzero rows in current block L(i,k) */
    lptr  = lptr0;
    luptr = luptr0;
    int Lnbrow, Rnbrow; /* number of nonzero rows in look-ahead window,
			   and remaining part.  */

    /*******************************************************************
     * Separating L blocks into the top part within look-ahead window
     * and the remaining ones.
     *******************************************************************/

	int lookAheadBlk=0, RemainBlk=0;

	tt_start = SuperLU_timer_();

	/* Sherry -- can this loop be threaded?? */
	/* Loop through all blocks in L(:,k) to set up pointers to the start
	* of each block in the data arrays.
	*   - lookAheadFullRow[i] := number of nonzero rows from block 0 to i
	*   - lookAheadStRow[i] := number of nonzero rows before block i
	*   - lookAhead_lptr[i] := point to the start of block i in L's index[]
	*   - (ditto Remain_Info[i])
	*/
	
	int* lookAhead_block_ids = SUPERLU_MALLOC(nlb * sizeof(int));
	int* Remain_block_ids = SUPERLU_MALLOC(nlb * sizeof(int));
    for (int i = 0; i < nlb; ++i) {
		ib = lsub[lptr];            /* Block number of L(i,k). */
		temp_nbrow = lsub[lptr+1];  /* Number of full rows. */

		int look_up_flag = 1; /* assume ib is outside look-up window */
		for (int j = k0+1; j < SUPERLU_MIN (k0 + num_look_aheads + 2, nsupers ); ++j) {
			if ( ib == perm_c_supno[j] ) {
				look_up_flag = 0; /* flag ib within look-up window */
				break;            /* Sherry -- can exit the loop?? */
			}
		}

		if ( look_up_flag == 0 ) { /* ib is within look-up window */
			if (lookAheadBlk==0) {
				lookAheadFullRow[lookAheadBlk] = temp_nbrow;
			} else {
				lookAheadFullRow[lookAheadBlk] = temp_nbrow + lookAheadFullRow[lookAheadBlk-1];
			}
			lookAheadStRow[lookAheadBlk] = cum_nrow;
			lookAhead_lptr[lookAheadBlk] = lptr;
			lookAhead_ib[lookAheadBlk] = ib;

			lookAhead_block_ids[lookAheadBlk] = i;

			lookAheadBlk++;
		} else { /* ib is not in look-up window */
			if ( RemainBlk==0 ) {
				Remain_info[RemainBlk].FullRow = temp_nbrow;
			} else {
				Remain_info[RemainBlk].FullRow =
					temp_nbrow + Remain_info[RemainBlk-1].FullRow;
			}
				RemainStRow[RemainBlk] = cum_nrow;
				// Remain_lptr[RemainBlk] = lptr;
			Remain_info[RemainBlk].lptr = lptr;
			// Remain_ib[RemainBlk] = ib;
			Remain_info[RemainBlk].ib = ib;

			Remain_block_ids[RemainBlk] = i;

			RemainBlk++;
		}

		cum_nrow += temp_nbrow;
		lptr += LB_DESCRIPTOR;  /* Skip descriptor. */
		lptr += temp_nbrow;     /* Move to next block */
		luptr += temp_nbrow;
    }  /* for i ... set up pointers for all blocks in L(:,k) */

	for(int i = 1;i<lookAheadBlk;i++){
		assert(lookAhead_block_ids[i] = lookAhead_block_ids[i-1] + 1);
	}
	for(int i = 1;i<RemainBlk;i++){
		assert(Remain_block_ids[i] = Remain_block_ids[i-1] + 1);
	}	

	SUPERLU_FREE(lookAhead_block_ids);
	SUPERLU_FREE(Remain_block_ids);

    lptr = lptr0;
    luptr = luptr0;
	/* leading dimension of L look-ahead buffer, same as Lnbrow */
	//int LDlookAhead_LBuff = lookAheadBlk==0 ? 0 :lookAheadFullRow[lookAheadBlk-1];
	Lnbrow = lookAheadBlk==0 ? 0 : lookAheadFullRow[lookAheadBlk-1];
	/* leading dimension of L remaining buffer, same as Rnbrow */
	//int LDRemain_LBuff = RemainBlk==0 ? 0 : Remain_info[RemainBlk-1].FullRow;
	Rnbrow = RemainBlk==0 ? 0 : Remain_info[RemainBlk-1].FullRow;
	/* assert( cum_nrow == (LDlookAhead_LBuff + LDRemain_LBuff) );*/
	/* Piyush fix */
	//int LDlookAhead_LBuff = lookAheadBlk==0? 0 : lookAheadFullRow[lookAheadBlk-1];

	nbrow = Lnbrow + Rnbrow; /* total number of rows in L */
	LookAheadRowSepMOP += 2*knsupc*(nbrow);

	/***********************************************
	* Gather U blocks (AFTER LOOK-AHEAD WINDOW)   *
	***********************************************/
	tt_start = SuperLU_timer_();

	int u_block_number = nub - jj0;
	
#if defined(USE_SW) && defined(OPT_gather_avoid)
	zMatrix_t U_all;
	zMatrix_t L_lookAhead_all;
	zMatrix_t L_Remain_all;
	zMatrix_t* U_list = SUPERLU_MALLOC((u_block_number)*sizeof(zMatrix_t));
	zMatrix_t* L_lookAhead_list = SUPERLU_MALLOC(lookAheadBlk*sizeof(zMatrix_t));
	zMatrix_t* L_Remain_list = SUPERLU_MALLOC(RemainBlk*sizeof(zMatrix_t));
#endif
	
	if ( nbrow > 0 ) { /* L(:,k) is not empty */
	/*
	* Counting U blocks
	*/
		ldu = 0; /* Calculate ldu for U(k,:) after look-ahead window. */
		ncols = 0; /* Total number of nonzero columns in U(k,:) */
		int temp_ncols = 0;

	/* jj0 contains the look-ahead window that was updated in
	dlook_ahead_update.c. Now the search can continue from that point,
	not to start from block 0. */


	/* if ( iam==0 ) printf("--- k0 %d, k %d, jj0 %d, nub %d\n", k0, k, jj0, nub);*/

	/*
	 * Loop through all blocks in U(k,:) to set up pointers to the start
	 * of each block in the data arrays, store them in Ublock_info[j]
	 * for block U(k,j).
	*/

		for (j = jj0; j < nub; ++j) { /* jj0 starts after look-ahead window. */
			temp_ncols = 0;
			/* Cannot remove following call, since perm_u != Identity  */
			arrive_at_ublock(
					j, &iukp, &rukp, &jb, &ljb, &nsupc,
					iukp0, rukp0, usub, perm_u, xsup, grid
					);
			Ublock_info[j].iukp = iukp;
			Ublock_info[j].rukp = rukp;
			Ublock_info[j].jb = jb;
			/* if ( iam==0 )
			printf("j %d: Ublock_info[j].iukp %d, Ublock_info[j].rukp %d,"
			"Ublock_info[j].jb %d, nsupc %d\n",
			j, Ublock_info[j].iukp, Ublock_info[j].rukp,
			Ublock_info[j].jb, nsupc); */

			/* Prepare to call GEMM. */
			jj = iukp;
			for (; jj < iukp+nsupc; ++jj) {
				segsize = klst - usub[jj];
				assert(segsize == 0 || ldu == 0 || ldu == segsize); // related to OPT_gather_avoid
				if ( segsize ) {
					++temp_ncols;
					if ( segsize > ldu ) ldu = segsize;
				}
			}
			Ublock_info[j].full_u_cols = temp_ncols;
			ncols += temp_ncols;
			
#if defined(USE_SW) && defined(OPT_gather_avoid)
			zMatrix_init(&U_list[j-jj0], ldu,temp_ncols, ldu, &uval[rukp]);
			// printf("(iukp, rukp, jb, full_u_cols), (%d %d %d %d)\n",Ublock_info[j].iukp, Ublock_info[j].rukp, Ublock_info[j].jb, Ublock_info[j].full_u_cols);
#endif
		} /* end for j ... compute ldu & ncols */

		// printf("(ldu, ncols, klst, nsupc), (%d, %d, %d, %d)\n",ldu, ncols, klst, nsupc);

		/* Now doing prefix sum on full_u_cols.
		* After this, full_u_cols is the number of nonzero columns
		* from block 0 to block j.
		*/
		for ( j = jj0+1; j < nub; ++j) {
			Ublock_info[j].full_u_cols += Ublock_info[j-1].full_u_cols;
		}

	 	/* Padding zeros to make {m,n,k} multiple of vector length. */
		jj = 8; //n;
		// printf("(gemm_padding, Rnbrow, ncols, ldu, jj), (%d, %d, %d, %d, %d)\n",gemm_padding, Rnbrow, ncols, ldu, jj);
		
		assert(gemm_padding == 0);
		if (gemm_padding > 0 && Rnbrow > jj && ncols > jj && ldu > jj) {

			gemm_n_pad = ncols + (ncols % GEMM_PADLEN);
			//gemm_n_pad = ncols;
			//gemm_k_pad = ldu + (ldu % GEMM_PADLEN);
			gemm_k_pad = ldu;

			for (i = Rnbrow; i < gemm_m_pad; ++i)  // padding A matrix
				for (j = 0; j < gemm_k_pad; ++j)
					Remain_L_buff[i + j*gemm_m_pad] = zero;
			for (i = 0; i < Rnbrow; ++i)
				for (j = ldu; j < gemm_k_pad; ++j)
					Remain_L_buff[i + j*gemm_m_pad] = zero;
			for (i = ldu; i < gemm_k_pad; ++i)     // padding B matrix
				for (j = 0; j < gemm_n_pad; ++j)
					bigU[i + j*gemm_k_pad] = zero;
			for (i = 0; i < ldu; ++i)
				for (j = ncols; j < gemm_n_pad; ++j)
					bigU[i + j*gemm_k_pad] = zero;
		} else {
			gemm_m_pad = Rnbrow;
			gemm_n_pad = ncols;
			gemm_k_pad = ldu;
		}
		// printf("(gemm_m_pad, gemm_n_pad, gemm_k_pad), (%d, %d, %d)\n",gemm_m_pad, gemm_n_pad, gemm_k_pad);

		tempu = bigU; /* buffer the entire row block U(k,:) */

         /* Gather U(k,:) into buffer bigU[] to prepare for GEMM */

		// check whether need gather?
#if defined(USE_SW) && defined(OPT_gather_avoid)
		// no need to gather
		zMatrix_init(&U_all,ldu,ncols,ldu,&uval[Ublock_info[jj0].rukp]);
#else

		if(ldu > 0){
			for (j = jj0; j < nub; ++j) { /* jj0 starts after look-ahead window. */
				if (j==jj0) 
					tempu = bigU;
				//else tempu = bigU + ldu * Ublock_info[j-1].full_u_cols;
				else 
					tempu = bigU + gemm_k_pad * Ublock_info[j-1].full_u_cols;
				/* == processing each of the rema ining columns in parallel == */
				iukp = Ublock_info[j].iukp;
				rukp = Ublock_info[j].rukp;
				jb = Ublock_info[j].jb;
				nsupc = SuperSize (jb );
				/* Copy from U(k,j) to tempu[], padding zeros.  */
				
				for (jj = iukp; jj < iukp+nsupc; ++jj) {
					segsize = klst - usub[jj];
					if ( segsize ) {
						lead_zero = ldu - segsize;
						assert(lead_zero == 0);
						for (i = 0; i < lead_zero; ++i) 
							tempu[i] = zero;
						//tempu += lead_zero;
						for (i = 0; i < segsize; ++i)
							tempu[i+lead_zero] = uval[rukp+i];
						rukp += segsize;
						tempu += gemm_k_pad;
					}
				}

			}   /* parallel for j = jj0 .. nub */
		}
#endif
        GatherMOP += 2*ldu*ncols;

    }  /* end if (nbrow>0), end gather U blocks */

    GatherUTimer += SuperLU_timer_() - tt_start;
    int jj_cpu = nub;       /* limit between CPU and GPU */
    int thread_id;
    /*tempv = bigV;*/

    /**********************
     * Gather L blocks    *
     **********************/
    tt_start = SuperLU_timer_();

    /* Loop through the look-ahead blocks to copy Lval into the buffer */
	if(ldu > 0){
#if defined(USE_SW) && defined(OPT_gather_avoid)
		zMatrix_init(&L_lookAhead_all, Lnbrow, ldu, nsupr, &lusup[luptr + lookAheadStRow[0] + (knsupc-ldu) * nsupr]);
#endif
		for (i = 0; i < lookAheadBlk; ++i) {
			int StRowDest, temp_nbrow;
			if ( i==0 ) {
				StRowDest = 0;
				temp_nbrow = lookAheadFullRow[0];
			} else {
				StRowDest   = lookAheadFullRow[i-1];
				temp_nbrow  = lookAheadFullRow[i]-lookAheadFullRow[i-1];
			}
			int StRowSource = lookAheadStRow[i];

#if defined(USE_SW) && defined(OPT_gather_avoid)
			// L lookAhead no copy
			zMatrix_init(&L_lookAhead_list[i],temp_nbrow,ldu,nsupr,&lusup[luptr + StRowSource + (knsupc-ldu) * nsupr]);
#else
			/* Now copying one block into L lookahead buffer */
			/* #pragma omp parallel for (gives slow down) */
			// printf("knsupc, ldu, knsupc-ldu : %d, %d, %d\n",knsupc,ldu,knsupc-ldu);
			for (j = knsupc-ldu; j < knsupc; ++j) { 
				tempu = &lookAhead_L_buff[StRowDest + (j - (knsupc-ldu)) * Lnbrow];
				tempv = &lusup[luptr + StRowSource + j * nsupr ];
				for (jj = 0; jj < temp_nbrow; ++jj) tempu[jj] = tempv[jj];
			} /* end for j ... */
#endif
			// printf("temp_nbrow, nsupr : %d, %d\n",temp_nbrow,nsupr);
		} /* parallel for i ... gather Lval blocks from lookahead window */
#if defined(USE_SW) && defined(OPT_gather_avoid)
		zMatrix_init(&L_Remain_all, Rnbrow, ldu, nsupr, &lusup[luptr + RemainStRow[0] + (knsupc-ldu) * nsupr]);
#endif
// 		/* Loop through the remaining blocks to copy Lval into the buffer */
		for (int i = 0; i < RemainBlk; ++i) {
			int StRowDest, temp_nbrow;
			if ( i==0 )  {
				StRowDest  = 0;
				temp_nbrow = Remain_info[0].FullRow;
			} else  {
				StRowDest   = Remain_info[i-1].FullRow;
				temp_nbrow  = Remain_info[i].FullRow - Remain_info[i-1].FullRow;
			}
			int StRowSource = RemainStRow[i];

#if defined(USE_SW) && defined(OPT_gather_avoid)
			zMatrix_init(&L_Remain_list[i],temp_nbrow,ldu,nsupr, &lusup[luptr + StRowSource + (knsupc-ldu) * nsupr]);
#else
			/* Now copying a block into L remaining buffer */
			// #pragma omp parallel for (gives slow down)
			for (int j = knsupc-ldu; j < knsupc; ++j) {
				tempu = &Remain_L_buff[StRowDest + (j - (knsupc-ldu)) * gemm_m_pad];
				tempv = &lusup[luptr + j*nsupr + StRowSource];
				for (jj = 0; jj < temp_nbrow; ++jj) tempu[jj] = tempv[jj];
			} /* end for j ... */
#endif
		} /* parallel for i ... copy Lval into the remaining buffer */
	}

	tt_end = SuperLU_timer_();
	GatherLTimer += tt_end - tt_start;


	/*************************************************************************
	* Perform GEMM (look-ahead L part, and remain L part) followed by Scatter
	*************************************************************************/
	tempu = bigU;  /* setting to the start of padded U(k,:) */

	if ( Lnbrow>0 && ldu>0 && ncols>0 ) { /* Both L(:,k) and U(k,:) nonempty */
		/***************************************************************
		* Updating blocks in look-ahead window of the LU(look-ahead-rows,:)
		***************************************************************/

		/* Count flops for total GEMM calls */
		ncols = Ublock_info[nub-1].full_u_cols;
		flops_t flps = 8.0 * (flops_t)Lnbrow * ldu * ncols;
		LookAheadScatterMOP += 3 * Lnbrow * ncols; /* scatter-add */
		schur_flop_counter += flps;
		stat->ops[FACT]    += flps;
		LookAheadGEMMFlOp  += flps;

#if ( PRNTlevel>=1 )
		gemm_max_m = SUPERLU_MAX(gemm_max_m, Lnbrow);
		gemm_max_n = SUPERLU_MAX(gemm_max_n, ncols);
		gemm_max_k = SUPERLU_MAX(gemm_max_k, ldu);
		tt_start = SuperLU_timer_();
#endif


#if defined(USE_SW) && defined(OPT_gather_avoid)
		assert(L_lookAhead_all.col == U_all.row);
		// perform gemm 
		zgemm_("N", "N", &L_lookAhead_all.row, &U_all.col, &L_lookAhead_all.col, &alpha,
				L_lookAhead_all.val, &L_lookAhead_all.ld,
				U_all.val, &U_all.ld, &beta, bigV, &Lnbrow, 1, 1);
#else
		zgemm_("N", "N", &Lnbrow, &ncols, &ldu, &alpha,
				lookAhead_L_buff, &Lnbrow,
				bigU, &ldu, &beta, bigV, &Lnbrow, 1, 1);
#endif

#if ( PRNTlevel>=1 )
		tt_end = SuperLU_timer_();
		LookAheadGEMMTimer += tt_end - tt_start;
#if ( PROFlevel>=1 )
		gemm_stats[gemm_count].m = Lnbrow;
		gemm_stats[gemm_count].n = ncols;
		gemm_stats[gemm_count].k = ldu;
		gemm_stats[gemm_count++].microseconds = (tt_end - tt_start) * 1e6;
#endif
		tt_start = SuperLU_timer_();
#endif

#if defined(USE_SW) && defined(OPT_lookAhead_scatter_slave)
		lookAhead_scatter_param_t lookAhead_scatter_param;
		lookAhead_scatter_param.l_lookAhead_block_num = lookAheadBlk;
		lookAhead_scatter_param.u_block_num = nub-jj0;
		lookAhead_scatter_param.u_start = jj0;
		lookAhead_scatter_param.ldt = ldt;
		lookAhead_scatter_param.Ublock_info = Ublock_info;
		lookAhead_scatter_param.lookAhead_lptr = lookAhead_lptr;
		lookAhead_scatter_param.lookAhead_ib = lookAhead_ib;
		lookAhead_scatter_param.lookAheadFullRow = lookAheadFullRow;
		lookAhead_scatter_param.bigV = bigV;
		lookAhead_scatter_param.xsup = xsup;
		lookAhead_scatter_param.klst = klst;
		lookAhead_scatter_param.nbrow = Lnbrow;
		lookAhead_scatter_param.usub = usub;
		lookAhead_scatter_param.lsub = lsub;
		lookAhead_scatter_param.Ufstnz_br_ptr = Ufstnz_br_ptr;
		lookAhead_scatter_param.Unzval_br_ptr = Unzval_br_ptr;
		lookAhead_scatter_param.Lrowind_bc_ptr = Lrowind_bc_ptr;
		lookAhead_scatter_param.Lnzval_bc_ptr = Lnzval_bc_ptr;
		lookAhead_scatter_param.grid = grid;
    	athread_spawn(lookAhead_scatter_0_naive, &lookAhead_scatter_param);
		athread_join();
#else
	   	/* Each thread is assigned one loop index ij, responsible for
	      block update L(lb,k) * U(k,j) -> tempv[]. */
		// printf("\nIn look ahead :\n");
	   	for (int ij = 0; ij < lookAheadBlk*(nub-jj0); ++ij) {
	       	/* jj0 starts after look-ahead window. */
            int j   = ij/lookAheadBlk + jj0;
            int lb  = ij%lookAheadBlk;

            /* Getting U block U(k,j) information */
            /* unsigned long long ut_start, ut_end; */
            int_t rukp =  Ublock_info[j].rukp;
            int_t iukp =  Ublock_info[j].iukp;
            int jb   =  Ublock_info[j].jb;
            int nsupc = SuperSize(jb);
            int ljb = LBj (jb, grid);  /* destination column block */
            int st_col;
            int ncols;  /* Local variable counts only columns in the block */
            if ( j > jj0 ) { /* jj0 starts after look-ahead window. */
                ncols  = Ublock_info[j].full_u_cols-Ublock_info[j-1].full_u_cols;
                st_col = Ublock_info[j-1].full_u_cols;
            } else {
                ncols  = Ublock_info[j].full_u_cols;
                st_col = 0;
            }

            /* Getting L block L(i,k) information */
            int_t lptr = lookAhead_lptr[lb];
            int ib   = lookAhead_ib[lb];
            int temp_nbrow = lsub[lptr+1];
            lptr += LB_DESCRIPTOR;
            int cum_nrow = (lb==0 ? 0 : lookAheadFullRow[lb-1]);

	    	/* Block-by-block GEMM in look-ahead window */
	    	doublecomplex* tempv1 = bigV + st_col * Lnbrow + cum_nrow;

            if ( ib < jb ) {
#if defined(USE_SW) && defined(OPT_scatter_index_compress)
                zscatter_u_sw (
					ib, jb,
					nsupc, iukp, xsup,
					klst, Lnbrow,
					lptr, temp_nbrow, lsub,
					usub, tempv1,
					Ufstnz_br_ptr, Unzval_br_ptr,
					grid
			    );
#else
                zscatter_u (
					ib, jb,
					nsupc, iukp, xsup,
					klst, Lnbrow,
					lptr, temp_nbrow, lsub,
					usub, tempv1,
					Ufstnz_br_ptr, Unzval_br_ptr,
					grid
			    );
#endif

            } else {
#if defined(USE_SW) && defined(OPT_scatter_index_compress)
                zscatter_l_sw (
					ib, ljb,
					nsupc, iukp, xsup,
					klst, Lnbrow,
					lptr, temp_nbrow,
					usub, lsub, tempv1,
					indirect, indirect2,
					Lrowind_bc_ptr, Lnzval_bc_ptr,
					grid
				);
#else
                zscatter_l (
					ib, ljb,
					nsupc, iukp, xsup,
					klst, Lnbrow,
					lptr, temp_nbrow,
					usub, lsub, tempv1,
					indirect, indirect2,
					Lrowind_bc_ptr, Lnzval_bc_ptr,
					grid
				);
#endif
			}
	   } /* end omp for ij = ... */
	
#endif

#if ( PRNTlevel>=1 )
		LookAheadScatterTimer += SuperLU_timer_() - tt_start;
#endif
	} /* end if Lnbrow>0 ... look-ahead GEMM and scatter */


    /***************************************************************
     * Updating remaining rows and columns on CPU.
     ***************************************************************/
    ncols = jj_cpu==0 ? 0 : Ublock_info[jj_cpu-1].full_u_cols;

    if ( Rnbrow>0 && ldu>0 ) { /* There are still blocks remaining ... */
		double flps = 8.0 * (double)Rnbrow * ldu * ncols;
		schur_flop_counter  += flps;
		stat->ops[FACT]     += flps;

#if ( PRNTlevel>=1 )
		RemainGEMM_flops += flps;
		gemm_max_m = SUPERLU_MAX(gemm_max_m, Rnbrow);
		gemm_max_n = SUPERLU_MAX(gemm_max_n, ncols);
		gemm_max_k = SUPERLU_MAX(gemm_max_k, ldu);
		tt_start = SuperLU_timer_();
		/* printf("[%d] .. k0 %d, before large GEMM: %d-%d-%d, RemainBlk %d\n",
			iam, k0,Rnbrow,ldu,ncols,RemainBlk);  fflush(stdout);
		assert( Rnbrow*ncols < bigv_size ); */
#endif

#if defined(USE_SW) && defined(OPT_gather_avoid) 
		assert(L_Remain_all.col = U_all.row);
		zgemm_("N", "N", &L_Remain_all.row, &U_all.col, &L_Remain_all.col, &alpha,
				L_Remain_all.val, &L_Remain_all.ld,
				U_all.val, &U_all.ld, &beta, bigV, &Rnbrow, 1, 1);
#else
		zgemm_("N", "N", &gemm_m_pad, &gemm_n_pad, &gemm_k_pad, &alpha,
				Remain_L_buff, &gemm_m_pad,
				bigU, &gemm_k_pad, &beta, bigV, &Rnbrow, 1, 1);
#endif

#if ( PRNTlevel>=1 )
		tt_end = SuperLU_timer_();
		RemainGEMMTimer += tt_end - tt_start;
#if ( PROFlevel>=1 )
		//fprintf(fgemm, "%8d%8d%8d %16.8e\n", Rnbrow, ncols, ldu,
		// (tt_end - tt_start)*1e6); // time in microsecond
		//fflush(fgemm);
		gemm_stats[gemm_count].m = Rnbrow;
		gemm_stats[gemm_count].n = ncols;
		gemm_stats[gemm_count].k = ldu;
		gemm_stats[gemm_count++].microseconds = (tt_end - tt_start) * 1e6;
#endif
		tt_start = SuperLU_timer_();
#endif

		/* Scatter into destination block-by-block. */

		/* Ideally, should organize the loop as:
				for (j = 0; j < jj_cpu; ++j) {
				for (lb = 0; lb < RemainBlk; ++lb) {
					L(lb,k) X U(k,j) -> tempv[]
					}
				}
			But now, we use collapsed loop to achieve more parallelism.
			Total number of block updates is:
			(# of RemainBlk in L(:,k)) X (# of blocks in U(k,:))
		*/


		/* Each thread is assigned one loop index ij, responsible for
			block update L(lb,k) * U(k,j) -> tempv[]. */

#if defined(USE_SW) && defined(OPT_remain_scatter_slave)
		remain_scatter_param_t remain_scatter_param;
		remain_scatter_param.l_remain_block_num = RemainBlk;
		remain_scatter_param.u_block_num = jj_cpu-jj0;
		remain_scatter_param.u_start = jj0;
		remain_scatter_param.ldt = ldt;
		remain_scatter_param.Ublock_info = Ublock_info;
		remain_scatter_param.Remain_info = Remain_info;
		remain_scatter_param.bigV = bigV;
		remain_scatter_param.xsup = xsup;
		remain_scatter_param.klst = klst;
		remain_scatter_param.nbrow = Rnbrow;
		remain_scatter_param.usub = usub;
		remain_scatter_param.lsub = lsub;
		remain_scatter_param.Ufstnz_br_ptr = Ufstnz_br_ptr;
		remain_scatter_param.Unzval_br_ptr = Unzval_br_ptr;
		remain_scatter_param.Lrowind_bc_ptr = Lrowind_bc_ptr;
		remain_scatter_param.Lnzval_bc_ptr = Lnzval_bc_ptr;
		remain_scatter_param.grid = grid;
    	// athread_spawn(remain_scatter_0_naive, &remain_scatter_param);
    	// athread_spawn(remain_scatter_1_ldm, &remain_scatter_param);
    	athread_spawn(remain_scatter_2_async, &remain_scatter_param);
		athread_join();
#else
		// printf("(ldt):%d\n",ldt);
		// printf("\nIn Remain :\n");
		// printf("(RemainBlk, jj_cpu-jj0):%d,%d\n",RemainBlk,jj_cpu-jj0);
		for (int ij = 0; ij < RemainBlk*(jj_cpu-jj0); ++ij) {
			/* jj_cpu := nub, jj0 starts after look-ahead window. */
			int j   = ij / RemainBlk + jj0; /* j-th block in U panel */
			int lb  = ij % RemainBlk;       /* lb-th block in L panel */

			/* Getting U block U(k,j) information */
			/* unsigned long long ut_start, ut_end; */
			int_t rukp =  Ublock_info[j].rukp;
			int_t iukp =  Ublock_info[j].iukp;
			int jb   =  Ublock_info[j].jb;
			int nsupc = SuperSize(jb);
			int ljb = LBj (jb, grid);
			int st_col;
			int ncols;
			if ( j>jj0 ) {
				ncols = Ublock_info[j].full_u_cols - Ublock_info[j-1].full_u_cols;
				st_col = Ublock_info[j-1].full_u_cols;
			} else {
				ncols = Ublock_info[j].full_u_cols;
				st_col = 0;
			}

			/* Getting L block L(i,k) information */
			int_t lptr = Remain_info[lb].lptr;
			int ib   = Remain_info[lb].ib;
			int temp_nbrow = lsub[lptr+1];
			lptr += LB_DESCRIPTOR;
			int cum_nrow = (lb==0 ? 0 : Remain_info[lb-1].FullRow);


			// const zMatrix_t* L = &L_Remain_list[lb];
			// const zMatrix_t* U = &U_list[j-jj0];
			// assert(L->col == U->row);
			// printf("(m,n,k):(%d,%d,%d)\n",L->row,U->col,U->row);
			// printf("%d,%d,%d\n",L->row,U->col,U->row);

			/* tempv1 points to block(i,j) in bigV : LDA == Rnbrow */
			doublecomplex* tempv1 = bigV + (st_col * Rnbrow + cum_nrow); /* Sherry */

			// printf("[%d] .. before scatter: ib %d, jb %d, temp_nbrow %d, Rnbrow %d\n", iam, ib, jb, temp_nbrow, Rnbrow); fflush(stdout);
			/* Now scattering the block */

			if ( ib < jb ) {
#if defined(USE_SW) && defined(OPT_scatter_index_compress)
				zscatter_u_sw (
					ib, jb,
					nsupc, iukp, xsup,
					//klst, Rnbrow, /*** klst, temp_nbrow, Sherry */
					klst, Rnbrow, /*** klst, temp_nbrow, Sherry */
					lptr, temp_nbrow, /* row dimension of the block */
					lsub, usub, tempv1,
					Ufstnz_br_ptr, Unzval_br_ptr,
					grid
					);
#else
				zscatter_u (
					ib, jb,
					nsupc, iukp, xsup,
					//klst, Rnbrow, /*** klst, temp_nbrow, Sherry */
					klst, Rnbrow, /*** klst, temp_nbrow, Sherry */
					lptr, temp_nbrow, /* row dimension of the block */
					lsub, usub, tempv1,
					Ufstnz_br_ptr, Unzval_br_ptr,
					grid
					);
#endif
			} else {
#if defined(USE_SW) && defined(OPT_scatter_index_compress)
				zscatter_l_sw(
					ib, ljb,
					nsupc, iukp, xsup,
					//klst, temp_nbrow, Sherry
					klst, Rnbrow, /*** temp_nbrow, Sherry */
					lptr, temp_nbrow, /* row dimension of the block */
					usub, lsub, tempv1,
					indirect, indirect2,
					Lrowind_bc_ptr,Lnzval_bc_ptr,
					grid
					);
#else
				zscatter_l(
					ib, ljb,
					nsupc, iukp, xsup,
					//klst, temp_nbrow, Sherry
					klst, Rnbrow, /*** temp_nbrow, Sherry */
					lptr, temp_nbrow, /* row dimension of the block */
					usub, lsub, tempv1,
					indirect, indirect2,
					Lrowind_bc_ptr,Lnzval_bc_ptr,
					grid
					);
#endif
			}

		} /* end omp for (int ij =...) */

#endif

#if ( PRNTlevel>=1 )
		RemainScatterTimer += SuperLU_timer_() - tt_start;
#endif
    } /* end if Rnbrow>0 ... update remaining block */


#if defined(USE_SW) && defined(OPT_gather_avoid)
	SUPERLU_FREE(U_list);
	SUPERLU_FREE(L_lookAhead_list);
	SUPERLU_FREE(L_Remain_list);
#endif
}  /* end if L(:,k) and U(k,:) are not empty */
