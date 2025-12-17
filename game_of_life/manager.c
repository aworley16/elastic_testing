/* manager */ 
#include <stdio.h> 
#include "mpi.h" 
int main(int argc, char *argv[]) 
{ 
   MPI_Init(&argc, &argv); 
   int original_rank = -1;
   MPI_Comm_rank(MPI_COMM_WORLD, &original_rank);
   
   MPI_Comm everyone;           /* inter-communicator */    
   MPI_Comm universe;           /* intracommunicator containing MPI_COMM_WORLD and the spawned processes */	
 
   if(original_rank == 0){
	   MPI_Comm_spawn("worker.exe", MPI_ARGV_NULL, 1, 
                 MPI_INFO_NULL, 0, MPI_COMM_SELF, &everyone, 
                 MPI_ERRCODES_IGNORE);
				  
		int remote_size = -1;
		int parent_rank = -1;
		MPI_Comm_remote_size(everyone, &remote_size);
		MPI_Comm_rank(everyone, &parent_rank);
		printf("remote_size = %d \n", remote_size);
		printf("parent rank = %d\n", parent_rank);
		
		MPI_Barrier(everyone);
		
		//int message = 99;
		//int message2 = 44;
		//MPI_Send(&message, 1, MPI_INT, 0, 0, everyone);
		//MPI_Send(&message2, 1, MPI_INT, 1, 1, everyone);
   
		
		MPI_Intercomm_merge(everyone, 1, &universe);
		
		//see size and rank of universe;
		int universe_size = -1;
		int universe_rank = -1;
		MPI_Comm_rank(universe, &universe_rank);
		MPI_Comm_size(universe, &universe_size);
		printf("ORIGINAL_RANK 0, universe %d of %d \n", universe_rank, universe_size);
   }	
 
   MPI_Finalize(); 
   return 0; 
} 