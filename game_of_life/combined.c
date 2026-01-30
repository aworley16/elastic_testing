
#include <pmix.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>

int main(){
	//initialize PMix
	
	pmix_proc_t myproc;
	pmix_status_t ps;
	PMIx_Init(&myproc, NULL, 0);
	if(ps != PMIX_SUCCESS){
		fprintf(stderr, "PMIx_Init failed: %s\n", PMIx_Error_string(ps));return 1;
	}

	//start initial session to get ranks
	MPI_Session world_sess;
	int si = MPI_Session_init(MPI_INFO_NULL, MPI_ERRORS_RETURN, &world_sess);
	if (si != MPI_SUCCESS) {
		fprintf(stderr, "MPI_Session_init failed\n");return 1;
    }
	MPI_Group world_group;
	MPI_Group_from_session_pset(world_sess, "test", &world_group);
	MPI_Comm world_comm;
	MPI_Comm_create_from_group(world_group, "phase_comm", MPI_INFO_NULL, MPI_ERRORS_RETURN, &world_comm);
	
	int size;
	int rank;
	MPI_Comm_size(world_comm, &size);
	MPI_Comm_rank(world_comm, &rank);
	
	printf("HELLO from rank %d of %d\n", rank, size);
	fflush(stdout);
	
	//Check if current session is big enough
	//if so skip to computation
	//else kill session, spawn processes and recreate session
	MPI_Session_finalize(&world_sess);
	
/* 	
	MPI_Session_init(info, MPI_ERRORS_RETURN, &session);
	char pset_name[] = "user://phase_group"
	MPI_Group phase_group;
	
	MPI_Session_get_group(session, pset_name, &phase_group);
	MPI_Comm_create_from_group(phase_group, "phase_comm", MPI_INFO, &comm); */
	PMIx_Finalize(NULL, 0);
	return 0; 
}


/* int spawn(){
	
}

int pset_link(){	
	
	pmix_query_t query;
	PMIX_QUERY_CONSTRUCT(&query);
	query.keys = (char **)malloc(2 * sizeof(char *));
	query.keys[0] = strdup(PMIX_JOB_SIZE);
	query.keys[1] = NULL;
	
	PMIX_INFO_CREATE(query.qualifiers, 1);
    PMIX_INFO_LOAD(&query.qualifiers[0],PMIX_NSPACE, target_nspace, PMIX_STRING);
	

}

std::system("date /t"); */