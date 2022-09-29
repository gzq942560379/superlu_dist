#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include "sw/superlu_defs_ex.h"

void bcast_flat_tree_build(bcast_t* bcast, int* dests, int dests_size, int root, MPI_Comm comm){
    int proc_rank, proc_size;
    MPI_Comm_rank(comm, &proc_rank);
    MPI_Comm_size(comm, &proc_size);
    bcast->comm = comm;
    bcast->root = root;
    bcast->in = false;
    if(proc_rank == root){
        bcast->sendto_size = 0;
        if(dests_size > 0){
            for(int i = 0; i < dests_size; i++){
                if(dests[i] != proc_rank){
                    bcast->sendto_size+=1;
                }else{
                    bcast->in = true;
                }
            }
            bcast->sendto = SUPERLU_MALLOC(bcast->sendto_size * sizeof(int));
            int index = 0;
            for(int i = 0; i < dests_size; i++){
                if(dests[i] != proc_rank){
                    bcast->sendto[index++] = dests[i];
                }
            }
            assert(index == bcast->sendto_size);
        }
        else{
            bcast->sendto = NULL;
        }
        bcast->recvfrom = -1;
    }else{
        bcast->sendto_size = 0;
        bcast->sendto = NULL;
        bcast->recvfrom = root;
        bcast->in = false;
        for(int i = 0; i < dests_size; i++){
            if(proc_rank == dests[i])
                bcast->in = true;
        }
    }
}

void bcast_init(bcast_t* bcast){
    bcast->sendto_size = 0;
    bcast->sendto = NULL;
    bcast->recvfrom = -1;
    bcast->in = -1;
    bcast->root = -1;
}

void bcast_blocked_forward(bcast_t* bcast, void* buffer, int count, MPI_Datatype datatype, int tag, int* recv_count_p){

    int proc_rank, proc_size;
    MPI_Comm_rank(bcast->comm, &proc_rank);
    MPI_Comm_size(bcast->comm, &proc_size);

    int typesize;
    MPI_Type_size(datatype, &typesize);
    *recv_count_p = 0;

    // printf("my rank : %d\n", proc_rank);
    // printf("bcast root rank : %d\n", bcast->root);
    // fflush(stdout);

    if(proc_rank == bcast->root){
        if(bcast->in){
            *recv_count_p = count;
        }
        for(int i = 0; i < bcast->sendto_size; ++i){
            // printf("sending message to %d\n", bcast->sendto[i]);
            // fflush(stdout);
			MPI_Send(buffer, count, datatype, bcast->sendto[i], tag, bcast->comm);
            // printf("send complete to %d , count %d\n", bcast->sendto[i], count);
            // fflush(stdout);
        }
    }else{
        if(bcast->in){
            // printf("recving message from %d\n", bcast->recvfrom);
            // fflush(stdout);
            MPI_Status status;
		    MPI_Recv(buffer, count, datatype, bcast->recvfrom, tag, bcast->comm, &status);
			MPI_Get_count(&status, datatype, recv_count_p);
            // printf("recv complete from %d , count %d\n", bcast->recvfrom, *recv_count_p);
            // fflush(stdout);
        }

    }

}

void bcast_destory(bcast_t* bcast){
    if(bcast->sendto != NULL)
        SUPERLU_FREE(bcast->sendto);
}

void bcast_build_LB(bcast_t** LB_p, const egraph_t* graph, int_t nsupers, gridinfo_t* grid){
    int iam = grid->iam;
	int_t myrow = MYROW( iam, grid );
    
    bool* temp_dests_flag = SUPERLU_MALLOC(sizeof(bool) * grid->npcol);
    int* temp_dests = SUPERLU_MALLOC(sizeof(int) * grid->npcol);
    int temp_dests_size;
    // LB
    bcast_t* LB = SUPERLU_MALLOC(sizeof(bcast_t) * nsupers);
    for(int_t gbc = 0; gbc < nsupers; gbc++){
        bcast_init(&LB[gbc]);
        int_t gbc_root = PCOL(gbc, grid);
        // reset flag
        for(int pc = 0; pc < grid->npcol; ++pc){
            temp_dests_flag[pc] = false;
        }
        for(int_t ptr = graph->edge_ptr[gbc]; ptr < graph->edge_ptr[gbc + 1]; ++ptr){
            int_t gb_update = graph->edges[ptr];
            int_t dest_pc = PCOL(gb_update, grid);
            temp_dests_flag[dest_pc] = true;
        }
        temp_dests_size = 0;
        for(int pc = 0; pc < grid->npcol; ++pc){
            if(temp_dests_flag[pc] == true){
                // global rank
                temp_dests[temp_dests_size++] = PNUM(myrow, pc, grid);
            }
        }
        bcast_flat_tree_build(&LB[gbc], temp_dests, temp_dests_size, PNUM(myrow, gbc_root, grid), grid->comm);
    }
    *LB_p = LB;
    SUPERLU_FREE(temp_dests);
    SUPERLU_FREE(temp_dests_flag);
}

void bcast_build_UB(bcast_t** UB_p, const egraph_t* graph, int_t nsupers, gridinfo_t* grid){
    int iam = grid->iam;
	int_t mycol = MYCOL( iam, grid );

    bool* temp_dests_flag = SUPERLU_MALLOC(sizeof(bool) * grid->nprow);
    int* temp_dests = SUPERLU_MALLOC(sizeof(int) * grid->nprow);
    int temp_dests_size;
    // UB
    bcast_t* UB = SUPERLU_MALLOC(sizeof(bcast_t) * nsupers);
    for(int_t gbr = 0; gbr < nsupers; gbr++){
        bcast_init(&UB[gbr]);
        int_t gbr_root = PROW(gbr, grid);
        // reset flag
        for(int pr = 0; pr < grid->nprow; ++pr){
            temp_dests_flag[pr] = false;
        }
        for(int_t ptr = graph->edge_ptr[gbr]; ptr < graph->edge_ptr[gbr + 1]; ++ptr){
            int_t gb_update = graph->edges[ptr];
            int_t dest_pr = PROW(gb_update, grid);
            temp_dests_flag[dest_pr] = true;
        }
        temp_dests_size = 0;
        for(int pr = 0; pr < grid->nprow; ++pr){
            if(temp_dests_flag[pr] == true){
                // global rank
                temp_dests[temp_dests_size++] = PNUM(pr, mycol, grid);
            }
        }
        bcast_flat_tree_build(&UB[gbr], temp_dests, temp_dests_size, PNUM(gbr_root, mycol, grid), grid->comm);
    }
    *UB_p = UB;

    SUPERLU_FREE(temp_dests);
    SUPERLU_FREE(temp_dests_flag);
}

void bcast_destory_all(bcast_t* bcasts, int_t size){
    for(int_t i = 0; i < size; ++i){
        bcast_destory(&bcasts[i]);
    }
    SUPERLU_FREE(bcasts);
}