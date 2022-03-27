#pragma once

// colunm
typedef struct {
        doublecomplex* val;
        int row;
        int col;
        int ld;
} zMatrix_t;

static void zMatrix_init(zMatrix_t* matirx, int row, int col, int ld, doublecomplex* val){
        matirx->row = row;
        matirx->col = col;
        matirx->ld = ld;
        matirx->val = val;
}

static void zMatrix_set_val(zMatrix_t* matirx, doublecomplex* val){
        matirx->val = val;
}

void pzgstrs2_fugaku
(int_t k0, int_t k, Glu_persist_t * Glu_persist, gridinfo_t * grid,
 zLocalLU_t * Llu, Ublock_info_t *Ublock_info, SuperLUStat_t * stat);

 void 
zscatter_l_fugaku (
            int ib,    /* row block number of source block L(i,k) */
            int ljb,   /* local column block number of dest. block L(i,j) */
            int nsupc, /* number of columns in destination supernode */
            int_t iukp, /* point to destination supernode's index[] */
            int_t* xsup,
            int klst,
            int nbrow,  /* LDA of the block in tempv[] */
            int_t lptr, /* Input, point to index[] location of block L(i,k) */
            int temp_nbrow, /* number of rows of source block L(i,k) */
            int_t* usub,
            int_t* lsub,
            doublecomplex *tempv,
            int* indirect_thread,int* indirect2,
            int_t ** Lrowind_bc_ptr, doublecomplex **Lnzval_bc_ptr,
            gridinfo_t * grid);

void
zscatter_u_fugaku (int ib,
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
           gridinfo_t * grid);