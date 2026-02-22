#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>





int main(int argc, char* argv[])
{
	int rank, size, proc_name_length;
	int limit = atoi(argv[1]);
	double spawn_start, spawn_end;
	char name[MPI_MAX_PROCESSOR_NAME];
	
	MPI_Init(&argc, &argv);	
	MPI_Comm universe;
	MPI_Comm_dup(MPI_COMM_WORLD, &universe);
	MPI_Comm parent;
	MPI_Comm bridge;
	MPI_Comm_get_parent(&parent); //check if this is a spawned child processes
	MPI_Get_processor_name(name, &proc_name_length);
	if(parent == MPI_COMM_NULL){
		printf("com_size, spawn_time, proc_ids\n");
		int max_needed = limit-1;
		for(int i=0; i<=max_needed; i++){
			printf("%d,", i); fflush(stdout);
			spawn_start = MPI_Wtime();
			if(i<0){
				MPI_Comm_spawn(argv[0], &argv[1], i, MPI_INFO_NULL, 0, MPI_COMM_WORLD, &bridge, MPI_ERRCODES_IGNORE); 
				MPI_Intercomm_merge(bridge, 0, &universe); 
			}
			spawn_end = MPI_Wtime()-spawn_start;
			
			MPI_Comm_size(universe, &size);
			MPI_Comm_rank(universe, &rank);

			for(int j=0; j<size; j++){
				if(rank == j){
					printf("%s,", name); 
				}
				MPI_Barrier(universe);  //loop barrier 
			}
			//cleanup this run
			MPI_Comm_free(&universe);
			printf("\n");
		}
	}else{
		MPI_Intercomm_merge(parent, 0, &universe);
		MPI_Comm_size(universe, &size);
		MPI_Comm_rank(universe, &rank);

		for(int j=0; j<size; j++){
			if(rank == j){
				printf("%s,", name); 
			}
			MPI_Barrier(universe); //loop barrier 
		}
		MPI_Finalize();
		return 0;
	}
	MPI_Finalize();
	return 0;
	
}