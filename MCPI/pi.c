#include <mpi.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <stdio.h>

int make_sample(double* hit, double* samples){
	double x = (double)rand() / (RAND_MAX + 1.0);
	double y = (double)rand() / (RAND_MAX + 1.0);
	double dist = sqrt(pow(x,2) + pow(y,2)); 
	if(dist <= 1){(*hit)++;}
	(*samples)++;
	return 0; 
}

MPI_Comm spawn(char* program, int expand_num){
	MPI_Comm children; 
	MPI_Comm_spawn("./pi.out", MPI_ARGV_NULL, expand_num, MPI_INFO_NULL, 0, MPI_COMM_SELF, &children, MPI_ERRCODES_IGNORE);
	return children;
}

int main(int argc, char** argv){
	
	double hits = 0;  
	double samples = 0;
	double total_hits = 0; 
	double total_samples = 0;
	double start_time = MPI_Wtime();
	double time_elasped = MPI_Wtime();
	int expand = 0; 
	MPI_Init(&argc, &argv);
	MPI_Comm universe;
	MPI_Comm parent;
	MPI_Comm_get_parent(&parent);	
	
	if(parent == MPI_COMM_NULL){
		MPI_Comm_dup(MPI_COMM_WORLD, &universe);
		expand = atoi(argv[1]);
	}else{
		MPI_Comm_dup(parent, &universe);
	}
	
	int flag = 0;
	int rank;
	int size;
	MPI_Comm_size(universe, &size);     
	MPI_Comm_rank(universe, &rank);
	MPI_Request flag_request;
	MPI_Status status;
	MPI_Comm children = MPI_COMM_SELF;
	srand(time(NULL)+rank);

	//if not root, continue sampling until flag tripped. 
	if(rank != 0){
		MPI_Ibcast(&flag, 1, MPI_INT, 0, universe, &flag_request);
		while(!flag){
			make_sample(&hits, &samples);
			MPI_Test(&flag_request, &flag, &status);	
		}
	}else{
		//sample for 30s
		while(time_elasped < 10){
			make_sample(&hits, &samples);
			time_elasped = MPI_Wtime();
		}
		//after 60s spawn additional processes. 
		printf("rank 0 spawning kid\n");
		//MPI_Comm_spawn("./pi.out", MPI_ARGV_NULL, expand, MPI_INFO_NULL, 0, MPI_COMM_SELF, &children, MPI_ERRCODES_IGNORE);
		printf("kid spawned\n");
		//children = spawn(argv[0], expand);
		//sample for additional 30s
		while(time_elasped < 30){
			make_sample(&hits, &samples);
			time_elasped = MPI_Wtime();
		}
		printf("sending flag\n");
		//after time has elapsed signal other processes to stop
		flag = 1;
		MPI_Ibcast(&flag, 1, MPI_INT, 0, universe, &flag_request);
		MPI_Ibcast(&flag, 1, MPI_INT, 0, children, &flag_request);
	}

	printf("rank %d after loop\n", rank);
	//initial spawn results -- children also call this. 	
	MPI_Reduce(&hits, &total_hits, 1, MPI_DOUBLE, MPI_SUM, 0, universe);
	MPI_Reduce(&samples, &total_samples, 1, MPI_DOUBLE, MPI_SUM, 0, universe);
	
	printf("rank %d after initial reduce\n", rank);
	double base = total_hits/(total_samples-total_hits);
	double base_hits = total_hits; 
	double base_samples = total_samples;
	double expanded = 0; 
	double expanded_hits = 0; 
	double expanded_samples = 0;
	MPI_Barrier(children);
	//add in results from the children
	if(rank == 0){
/* 		double child_hits = 0; 
		double child_samples = 0; 
		MPI_Reduce(&expanded_hits, &child_hits, 1, MPI_DOUBLE, MPI_SUM, 0, children);
		MPI_Reduce(&expanded_hits, &child_samples, 1, MPI_DOUBLE, MPI_SUM, 0, children);
		expanded_hits = base_hits + child_hits;
		expanded_samples = base_samples + child_samples;
		expanded = expanded_hits/expanded_samples; */
		printf("original-hits samples pi:, %f, %f, %f, \n", base_hits, base_samples, base);
		printf("expanded-hits samples pi:, %f, %f, %f, \n", expanded_hits, expanded_samples, expanded);
	}
	
	MPI_Finalize();
}

