#include "mpi.h" 
#include <stdio.h> 
int main(int argc, char *argv[]) 
{ 
   int parent_size; 
   int rank = -1;
   int local_size;
   int child_rank;
   MPI_Comm parent; 
   
   MPI_Init(&argc, &argv); 
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   MPI_Comm_size(MPI_COMM_WORLD, &local_size);
   
   printf("CHILD CHECKPOINT - %d \n", rank);
   
   MPI_Comm_get_parent(&parent); 
   if (parent == MPI_COMM_NULL){printf("No parent!\n");} 
   else{
		MPI_Comm_remote_size(parent, &parent_size); 
		if (parent_size != 1) printf("Something's wrong with the parent\n"); 
		MPI_Comm_rank(parent, &child_rank);
   }
   
   printf("CHILD CHECKPOINT 2 - %d \n", rank);
   
   MPI_Barrier(parent);
   
   int message;
   //MPI_Recv(&message, 1, MPI_INT, 0, rank, parent, MPI_STATUS_IGNORE);
   
   //printf("%d - message received - %d \n", rank, message);
	MPI_Comm universe;
	
	MPI_Intercomm_merge(parent, 1, &universe);
    
	//MPI_Comm dupe = parent;
	
	
   	int universe_size = -1;
	int universe_rank = -1;
	MPI_Comm_rank(universe, &universe_rank);
	MPI_Comm_size(universe, &universe_size);

	printf("CHILD %d universe %d of %d \n", rank, universe_rank, universe_size);

   MPI_Finalize(); 
   return 0; 
} 