#include <stdio.h>
#include <stdlib.h>
#include "sw/superlu_defs_ex.h"

static int int_t_cmp(const void* a, const void* b)
{
    return (int)(*(const int_t*)a - *(const int_t*)b);
}

void egraph_global_build(egraph_t* graph_global, int_t n, Glu_persist_t* Glu_persist, int_t** Lrowind_bc_ptr, gridinfo_t * grid){
    // if(options->SymPattern != YES || options->RowPerm != NOROWPERM){
    //     fprintf(stderr, "egraph only support Symmetric matrix\n");
    //     fflush(stderr);
    //     abort();
    // }
    // ASSERT(options->SymPattern == YES);
    // ASSERT(options->RowPerm != NOROWPERM);

    int iam = grid->iam;
    int np = grid->nprow * grid->npcol;
    int_t nsupers = Glu_persist->supno[n-1] + 1;
    int_t ncb = CEILING(nsupers, grid->npcol);
    int_t max_ncb = (nsupers + grid->npcol - 1) / grid->npcol;
    int myrow = MYROW(iam, grid);
    int mycol = MYCOL(iam, grid);

    // build graph local
    egraph_t graph_local;
    graph_local.node = ncb;
    graph_local.edge_ptr = SUPERLU_MALLOC(sizeof(int_t) * (ncb + 1));
    int_t* edge_count_row = &graph_local.edge_ptr[1];

    // count local edge
    int_t total_local_edge = 0;
    for(int_t lbc = 0; lbc < ncb; lbc++){
        edge_count_row[lbc] = 0;
        int_t gbc = lbc * grid->npcol + mycol;
        int_t* lsub = Lrowind_bc_ptr[lbc];
        if(lsub){
            int_t nlb = lsub[0]; // number of L blocks in L(:,gb)
            int_t total_rows = lsub[1];
            int_t lsub_cur_ptr = BC_HEADER;
            for(int_t i = 0; i < nlb; i++){
                int_t gbr = lsub[lsub_cur_ptr];
                int_t rows = lsub[lsub_cur_ptr + 1];
                if(gbr != gbc){  // skip diagonal block
                    total_local_edge += 1;
                    edge_count_row[lbc] += 1;
                }
                lsub_cur_ptr += LB_DESCRIPTOR + rows;
            }
        }
    }
    // calculate edge ptr
    graph_local.edge_ptr[0] = 0;
    for(int_t lbc = 0; lbc < ncb; lbc++){
        graph_local.edge_ptr[lbc + 1] = edge_count_row[lbc] + graph_local.edge_ptr[lbc];
    }

    // alloc local graph space
    graph_local.edge = total_local_edge;
    graph_local.edges = SUPERLU_MALLOC(sizeof(int_t) * total_local_edge);

    // fill in edges
    int_t tmp_index = 0;
    for(int_t lbc = 0; lbc < ncb; lbc++){
        int_t gbc = lbc * grid->npcol + mycol;
        int_t* lsub = Lrowind_bc_ptr[lbc];
        if(lsub){
            int_t nlb = lsub[0]; // number of L blocks in L(:,gb)
            int_t total_rows = lsub[1];
            int_t lsub_cur_ptr = BC_HEADER;
            for(int_t i = 0; i < nlb; i++){
                int_t gbr = lsub[lsub_cur_ptr];
                int_t rows = lsub[lsub_cur_ptr + 1];
                if(gbr != gbc){  // skip diagonal block
                    graph_local.edges[tmp_index++] = gbr;
                }
                lsub_cur_ptr += LB_DESCRIPTOR + rows;
            }
        }
    }
    assert(tmp_index == total_local_edge);

    // check edges is sorted
    for(int_t u = 0; u < graph_local.node; ++u){
        if(graph_local.edge_ptr[u+1] - graph_local.edge_ptr[u] > 1){
            for(int_t i = graph_local.edge_ptr[u] + 1; u < graph_local.edge_ptr[u+1]; ++u){
                assert(graph_local.edges[i] > graph_local.edges[i-1]);
            }
        }
    }

    // // print
    // printf("graph local : \n");
    // for(int_t lbc = 0; lbc < ncb; lbc++){
    //     int_t gbc = lbc * grid->npcol + mycol;
    //     printf("%d : ", gbc);
    //     for(int_t ptr = graph_local.edge_ptr[lbc]; ptr < graph_local.edge_ptr[lbc + 1]; ++ptr){
    //         printf("%d ", graph_local.edges[ptr]);
    //     }
    //     printf("\n");
    //     // break;
    // }
    // fflush(stdout);

    int_t* temp_edge_ptrs = NULL;
    int_t* temp_edges = NULL;
    int* edges_counts_preproc = NULL;
    int* edges_displs_preproc = NULL;

    // Merge graph in diagonal process
    // Gatherv graph
    MPI_Comm col_comm = (grid->cscp).comm;
    egraph_t graph_col;
    graph_col.node = graph_local.node;
    MPI_Allreduce(&graph_local.edge, &graph_col.edge, 1, mpi_int_t, MPI_SUM, col_comm);
    temp_edge_ptrs = SUPERLU_MALLOC((ncb + 1) * grid->nprow * sizeof(int_t));
        
    MPI_Allgather(graph_local.edge_ptr, (ncb + 1), mpi_int_t, 
        temp_edge_ptrs, (ncb + 1), mpi_int_t, col_comm);

    edges_counts_preproc = SUPERLU_MALLOC(grid->nprow * sizeof(int));
    edges_displs_preproc = SUPERLU_MALLOC((grid->nprow + 1) * sizeof(int));
    temp_edges = SUPERLU_MALLOC(graph_col.edge * sizeof(int_t));
    edges_displs_preproc[0] = 0;
    for(int pr = 0; pr < grid->nprow; pr++){
        edges_counts_preproc[pr] = temp_edge_ptrs[pr * (ncb + 1) + ncb];
        edges_displs_preproc[pr+1] = edges_displs_preproc[pr] + edges_counts_preproc[pr];
    }
    assert(edges_displs_preproc[grid->nprow] == graph_col.edge);  

    MPI_Allgatherv(graph_local.edges, graph_local.edge, mpi_int_t, 
        temp_edges, edges_counts_preproc, edges_displs_preproc, mpi_int_t,
        col_comm);
        
    // merge edge pre proc 
    graph_col.edge_ptr = SUPERLU_MALLOC((ncb + 1) * sizeof(int_t));
    graph_col.edges = SUPERLU_MALLOC(graph_col.edge * sizeof(int_t));

    tmp_index = 0;
    graph_col.edge_ptr[0] = 0;
    for(int_t lbc = 0; lbc < ncb; lbc++){
        for(int pr = 0; pr < grid->nprow; ++pr){
            int_t* edges_proc = &temp_edges[edges_displs_preproc[pr]];
            int_t* edges_ptr_proc = &temp_edge_ptrs[pr * (ncb + 1)];
            for(int_t ptr = edges_ptr_proc[lbc]; ptr < edges_ptr_proc[lbc + 1]; ptr++){
                graph_col.edges[tmp_index++] = edges_proc[ptr];
            }        
        }
        // sort
        graph_col.edge_ptr[lbc + 1] = tmp_index;
        size_t nitems = (size_t)(graph_col.edge_ptr[lbc + 1] -  graph_col.edge_ptr[lbc]);
        qsort((void*)&graph_col.edges[graph_col.edge_ptr[lbc]], nitems, sizeof(int_t), int_t_cmp);
    }
    assert(tmp_index == graph_col.edge);

    // // print
    // printf("graph col : \n");
    // for(int_t lbc = 0; lbc < ncb; lbc++){
    //     int_t gbc = lbc * grid->npcol + mycol;
    //     printf("%d : ", gbc);
    //     for(int_t ptr = graph_col.edge_ptr[lbc]; ptr < graph_col.edge_ptr[lbc + 1]; ++ptr){
    //         printf("%d ", graph_col.edges[ptr]);
    //     }
    //     printf("\n");
    //     // break;
    // }      
    // fflush(stdout);

    // free space
    egraph_destroy(&graph_local);
    SUPERLU_FREE(temp_edge_ptrs);
    SUPERLU_FREE(edges_counts_preproc);
    SUPERLU_FREE(edges_displs_preproc);
    SUPERLU_FREE(temp_edges);

    // row comm node merge
    MPI_Comm row_comm = (grid->rscp).comm;
    
    graph_global->node = nsupers;
    MPI_Allreduce(&graph_col.edge, &graph_global->edge, 1, mpi_int_t, MPI_SUM, row_comm);
    

    temp_edge_ptrs = SUPERLU_MALLOC((ncb + 1) * grid->npcol * sizeof(int_t));
    MPI_Allgather(graph_col.edge_ptr, (ncb + 1), mpi_int_t, 
        temp_edge_ptrs, (ncb + 1), mpi_int_t, row_comm);


    edges_counts_preproc = SUPERLU_MALLOC(grid->npcol * sizeof(int));
    edges_displs_preproc = SUPERLU_MALLOC((grid->npcol + 1) * sizeof(int));    
    temp_edges = SUPERLU_MALLOC(graph_global->edge * sizeof(int_t));
    edges_displs_preproc[0] = 0;
    for(int pc = 0; pc < grid->npcol; ++pc){
        edges_counts_preproc[pc] = temp_edge_ptrs[pc * (ncb + 1) + ncb];
        edges_displs_preproc[pc+1] = edges_displs_preproc[pc] + edges_counts_preproc[pc];
    }
    assert(edges_displs_preproc[grid->npcol] == graph_global->edge);


    MPI_Allgatherv(graph_col.edges, graph_col.edge, mpi_int_t, 
        temp_edges, edges_counts_preproc, edges_displs_preproc, mpi_int_t,
        row_comm);   

    // merge node
    graph_global->edge_ptr = SUPERLU_MALLOC((nsupers + 1) * sizeof(int_t));
    graph_global->edges = SUPERLU_MALLOC(graph_global->edge * sizeof(int_t));
    tmp_index = 0;
    graph_global->edge_ptr[0] = 0;
    for(int_t k = 0; k < nsupers; ++k){
        int kcol = PCOL(k, grid);
        int lk = LBj(k, grid);
        int_t* edges_ptr_proc = &temp_edge_ptrs[kcol * (ncb + 1)];
        int_t* edges_proc = &temp_edges[edges_displs_preproc[kcol]];
        for(int_t ptr = edges_ptr_proc[lk]; ptr < edges_ptr_proc[lk + 1]; ptr++){
            graph_global->edges[tmp_index++] = edges_proc[ptr];
        }
        graph_global->edge_ptr[k + 1] = tmp_index;
    }    
    assert(tmp_index == graph_global->edge);


    
    // // print
    // if(!iam){
    //     printf("graph global : \n");
    //     for(int_t k = 0; k < nsupers; k++){
    //         printf("%d : ", k);
    //         for(int_t ptr = graph_global->edge_ptr[k]; ptr < graph_global->edge_ptr[k + 1]; ++ptr){
    //             printf("%d ", graph_global->edges[ptr]);
    //         }
    //         printf("\n");
    //     }
    //     fflush(stdout);
    // }

    // free space
    egraph_destroy(&graph_col);
    SUPERLU_FREE(temp_edge_ptrs);
    SUPERLU_FREE(edges_counts_preproc);
    SUPERLU_FREE(edges_displs_preproc);
    SUPERLU_FREE(temp_edges);
}

void egraph_global_reverse(egraph_t* graph_r, const egraph_t* graph){
    graph_r->node = graph->node;
    graph_r->edge = graph->edge;
    graph_r->edge_ptr = SUPERLU_MALLOC((graph_r->node+1) * sizeof(int_t));
    graph_r->edges = SUPERLU_MALLOC(graph_r->edge * sizeof(int_t));
    int_t* tmp_edge_index = SUPERLU_MALLOC((graph_r->node+1) * sizeof(int_t));
    int_t* tmp_edge_count = tmp_edge_index + 1;

    for(int_t i = 0; i < graph->node + 1; ++i){
        tmp_edge_index[i] = 0;
    }

    for(int_t u = 0; u < graph->node; ++u){
        for(int_t ptr = graph->edge_ptr[u]; ptr < graph->edge_ptr[u+1]; ++ptr){
            int_t v = graph->edges[ptr];
            tmp_edge_count[v]++;
        }
    }

    tmp_edge_index[0] = 0;
    graph_r->edge_ptr[0] = 0;
    for(int_t i = 0; i < graph_r->node; ++i){
        tmp_edge_index[i+1] = tmp_edge_index[i] + tmp_edge_count[i];
        graph_r->edge_ptr[i+1] = tmp_edge_index[i+1];
    }

    for(int_t u = 0; u < graph->node; ++u){
        for(int_t ptr = graph->edge_ptr[u]; ptr < graph->edge_ptr[u+1]; ++ptr){
            int_t v = graph->edges[ptr];
            graph_r->edges[tmp_edge_index[v]++] = u;
        }
    }

    // printf("print graph_global reverse : \n");
    // for(int_t u = 0; u < graph_r->node; ++u){
    //     printf("%d :", u);
    //     for(int_t ptr = graph_r->edge_ptr[u]; ptr < graph_r->edge_ptr[u+1]; ++ptr){
    //         int_t v = graph_r->edges[ptr];
    //         printf(" %d", v);
    //     }
    //     printf("\n");
    // }
    // fflush(stdout);

    SUPERLU_FREE(tmp_edge_index);
}


void egraph_destroy(egraph_t* graph){
    SUPERLU_FREE(graph->edge_ptr);
    SUPERLU_FREE(graph->edges);
}