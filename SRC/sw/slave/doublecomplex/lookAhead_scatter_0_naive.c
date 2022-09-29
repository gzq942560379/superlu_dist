#include <crts.h>
#include <stdio.h>
#include "sw/slave_param.h"

#ifndef MAX_BLOCK_SIZE
#define MAX_BLOCK_SIZE 256
#endif

__thread_local int indirect[MAX_BLOCK_SIZE];
__thread_local int indirect2[MAX_BLOCK_SIZE];

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

    int local_ij_start = (l_lookAhead_block_num * u_block_num * CRTS_tid)/64;
    int local_ij_end = (l_lookAhead_block_num * u_block_num * (CRTS_tid+1))/64;

    for(int ij = local_ij_start; ij<local_ij_end; ij++){
        int j = ij/l_lookAhead_block_num + u_start;
        int lb = ij%l_lookAhead_block_num;

        // Getting U block information
        int_t rukp =  Ublock_info[j].rukp;
        int_t iukp =  Ublock_info[j].iukp;
        int jb   =  Ublock_info[j].jb;
        int nsupc = SuperSize(jb);
        int ljb = LBj (jb, grid);
        int st_col = j > u_start ? Ublock_info[j-1].full_u_cols : 0;

        /* Getting L block L(i,k) information */
        int_t lptr = lookAhead_lptr[lb];
        int ib   = lookAhead_ib[lb];
        int temp_nbrow = lsub[lptr+1];
        lptr += LB_DESCRIPTOR;
        int cum_nrow = (lb == 0 ? 0 : lookAheadFullRow[lb-1]);

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

            // tempv = bigV + (cum_nrow + cum_ncol*nbrow);
            for (int_t jj = 0; jj < nsupc; ++jj) {
                int_t segsize = klst - usub[iukp + jj];
                int_t fnz = index[iuip_lib++];
                if (segsize) {          /* Nonzero segment in U(k,j). */
                    doublecomplex* ucol = &Unzval_br_ptr[lib][ruip_lib];
                    for (int i = 0; i < temp_nbrow; ++i) {
                        int_t rel = lsub[lptr + i] - fnz;
                        z_sub(&ucol[rel], &ucol[rel], &tempv[i]);
                    } /* for i = 0:temp_nbropw */
                    tempv += nbrow; /* Jump LDA to next column */
                }  /* if segsize */
                ruip_lib += ilst - fnz;
            }  /* for jj = 0:nsupc */

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

            int_t fnz = FstBlockC (ib);
            lptrj += LB_DESCRIPTOR;
            int_t dest_nbrow=index[lptrj - 1];

            for (int i = 0; i < dest_nbrow; ++i) {
                int_t rel = index[lptrj + i] - fnz;
                indirect[rel] = i;
            }

            /* can be precalculated? */
            for (int i = 0; i < temp_nbrow; ++i) { /* Source index is a subset of dest. */
                int_t rel = lsub[lptr + i] - fnz;
                indirect2[i] = indirect[rel];
            }
            
            doublecomplex* nzval = Lnzval_bc_ptr[ljb] + luptrj; /* Destination block L(i,j) */

            for (int_t jj = 0; jj < nsupc; ++jj) {
                int_t segsize = klst - usub[iukp + jj];
                if (segsize) {
                    for (int i = 0; i < temp_nbrow; ++i) {
                        z_sub(&nzval[indirect2[i]], &nzval[indirect2[i]], &tempv[i]);
                    }
                    tempv += nbrow;
                }
                nzval += ldv;
            }
        }
    }
    ldm_free(indirect,ldt*sizeof(int));
    ldm_free(indirect2,ldt*sizeof(int));

}