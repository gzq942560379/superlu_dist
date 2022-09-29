#include <mpi.h>
#include <stdio.h>

int get_tag_up(){
    int tag_ub;
	int mpi_flag;
    MPI_Comm_get_attr (MPI_COMM_WORLD, MPI_TAG_UB, (void*)&tag_ub, &mpi_flag);
	if(!mpi_flag){
        fprintf (stderr, "Could not get TAG_UB\n");
        return (-1);
	}
    return tag_ub;
}