#pragma once

#include "dcomplex.h"
#include "superlu_defs.h"
#include "superlu_zdefs.h"

typedef struct{
    int l_lookAhead_block_num;
    int u_block_num;
    int u_start;
    int ldt;
    Ublock_info_t* Ublock_info;
    int* lookAhead_lptr;
    int* lookAhead_ib;
    int* lookAheadFullRow;
    doublecomplex* bigV;
    int* xsup;
    int klst;
    int nbrow;
    int_t* usub;
    int_t* lsub;
    int_t ** Ufstnz_br_ptr;
    doublecomplex **Unzval_br_ptr;
    int_t ** Lrowind_bc_ptr;
    doublecomplex **Lnzval_bc_ptr;
    gridinfo_t * grid;
}lookAhead_scatter_param_t;

typedef struct{
    int l_remain_block_num;
    int u_block_num;
    int u_start;
    int ldt;
    Ublock_info_t* Ublock_info;
    Remain_info_t* Remain_info;
    doublecomplex* bigV;
    int* xsup;
    int klst;
    int nbrow;
    int_t* usub;
    int_t* lsub;
    int_t ** Ufstnz_br_ptr;
    doublecomplex **Unzval_br_ptr;
    int_t ** Lrowind_bc_ptr;
    doublecomplex **Lnzval_bc_ptr;
    gridinfo_t * grid;
}remain_scatter_param_t;
