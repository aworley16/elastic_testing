#define _DEFAULT_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <mpi.h>
#include<time.h>


int main(){

	//const char* pset_name = "mpi://WORLD"; //??? 

	int size;
	int rank;
	MPI_Session shandle = MPI_SESSION_NULL;
	MPI_Init(NULL, NULL);
	
	MPI_Comm bridge;
	MPI_Comm parent;
	MPI_Comm child;
	
	MPI_Comm_get_parent(&parent);
	if(parent != MPI_COMM_NULL) 
	{ 
		MPI_Intercomm_merge(parent, 0, &bridge); //merge with parent comm (current universe)
		MPI_Comm_size(bridge, &size);     
		MPI_Comm_rank(bridge, &rank);   
		printf("WORLD SPAWNED process %d of %d\n", rank, size); 
	}
    else //if original dup MPI_COMM_WORLD so we have a handle that we can manipulate 
	{
		MPI_Comm_spawn("./a.out", MPI_ARGV_NULL, 1, MPI_INFO_NULL, 0, MPI_COMM_WORLD, &child, MPI_ERRCODES_IGNORE);
		MPI_Intercomm_merge(child, 0, &bridge); //create new universe
		MPI_Comm_size(bridge, &size);     
		MPI_Comm_rank(bridge, &rank);   
		printf("WORLD original process %d of %d\n", rank, size); 
	}
	
	//spin up a session see if mpi://world was modded. 
	int err = MPI_Session_init(MPI_INFO_NULL, MPI_ERRORS_RETURN, &shandle);
	
	MPI_Group world_group;
	MPI_Comm new_comm;

	int i=0;
	int n_psets=0;
	int psetlen=0;

	char* pset_name = NULL;
	MPI_Session_get_num_psets(shandle, MPI_INFO_NULL, &n_psets);
	for(i=0, pset_name=NULL; i< n_psets; i++)
	{
		psetlen = 0;
		MPI_Session_get_nth_pset(shandle, MPI_INFO_NULL, i, &psetlen, NULL);
		pset_name = (char *)malloc(sizeof(char) * psetlen);
		MPI_Session_get_nth_pset(shandle, MPI_INFO_NULL, i, &psetlen, pset_name);
		printf(" pset -- ");
		for(int j=0; j<psetlen; j++){printf("%c", pset_name[j]);}
		printf(" \n");
		fflush(stdout);
		free(pset_name);
	}

	MPI_Group newgroup;
	MPI_Group_from_session_pset(shandle, "mpi://world", &newgroup);
    MPI_Comm_create_from_group(newgroup, "MyCommTag", MPI_INFO_NULL, MPI_ERRORS_RETURN, &new_comm);

    MPI_Comm_rank(new_comm, &rank);
	MPI_Comm_size(new_comm, &size);
	if(parent == MPI_COMM_NULL) {	
		printf("Hello from original rank %d of %d\n", rank, size);
	}
	else{
		printf("Hello from spawned rank %d of %d\n", rank, size);
	}
	
	
	//MPI_Group_free(&group);
	MPI_Session_finalize(&shandle);
	MPI_Finalize();
	return 0; 
}
