#define _DEFAULT_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <mpi.h>
#include <time.h>
#include <string.h>


int*Step(int**local, int** local_new, int N, int rows);
void square(int N);
void floater(int N);
void gosper_glider_gun(int N);
void rand_grid(int N, int fillRatio);
void print_matrix(int* A, int n, int rows);

void Halo(int* local, int N, int rows, int rank, MPI_Comm comm);

int *boardState;
int *sendcounts, *disp;

void cleanup(){
	if(boardState!=NULL){free(boardState);}
	if(sendcounts!=NULL){free(sendcounts);}
	if(disp!=NULL){free(disp);}
	return;
}

int* initalize_root_board(int N, char seed){
	switch (seed)
    {	
		
    case 's': square(N);
        break;

    case 'f': floater(N);
        break;

    case 'g':
        if(N<38)
        {
            printf("Too small grid, running square instead\n");
            square(N);
        }
        else
        {
            gosper_glider_gun(N);
        }
        break;
    
    case 'r': rand_grid(N, 0);
        break;

    default:
        printf("Unrecognised command, defaulting to square\n");
        square(N);
    }
	return boardState;
}

int read_file(char* filename, int x, int y){
	FILE* ptr = fopen(filename, "r");
	if(ptr==NULL){return -1;}
	char line[(x+2)*2];
	char* token;
	int offset = 0; 
	int k=0;
	while(fgets(line, (x+2)*2, ptr)!=NULL){		
		if(line[0] == '\n'){continue;} // c is stupid when it comes to \n in files. 
		token = strtok(line, ",");
		while(token!=NULL){
			boardState[offset] = atoi(token);
			offset++;
			token = strtok(NULL, ",");
		}
		k++;
    }
	fclose(ptr);
	return 0; 
}

int write_file(char* filename, int x, int y){
	FILE* ptr = fopen(filename, "w");
	if(ptr==NULL){exit(1);}

	for(int i=0; i<x+2; i++){
		for(int j=0; j<y+1; j++){
			fprintf(ptr, "%d,", boardState[i*(x+2)+j]);
		}
		fprintf(ptr, "%d\n", boardState[i*(x+2)+(y+1)]); //no comma at end of line
	}
	fclose(ptr);
	return 0;
}

//update universe, phase_comm, color, and makes sure head proc is known everywhere if it changed
int setup_comms(int N, int phase, int phase_size, MPI_Comm* universe, MPI_Comm* phase_comm, int* color, char** argv){
	
	int uni_size = -1;
	int uni_rank = -1;

	//If phase_comm exists and is same size as previous then just use previous setup
	if(*phase_comm != MPI_COMM_NULL)
	{
		int curr_size = -1;
		MPI_Comm_size(*phase_comm, &curr_size);
		if(phase_size == curr_size){
			return 1;	
		}
		else{
			MPI_Comm_free(phase_comm);
		}
	}
	
	//if not figure out what needs to happen & clear room for new phase_comm
	MPI_Comm_size(*universe, &uni_size);
	MPI_Comm_rank(*universe, &uni_rank);
	*color = 0; 
	
	//if additional processes needed, expand universe	
	if(phase_size > uni_size){
		//calculate and spawn processes as needed.
		int expand_num = phase_size - uni_size;
		MPI_Comm bridge;
		MPI_Comm new_uni;

		printf("Spawning %d processes \n", expand_num);
		MPI_Comm_spawn("./gol.exe", &argv[1], expand_num, MPI_INFO_NULL, 0, *universe, &bridge, MPI_ERRCODES_IGNORE);
		MPI_Bcast(&phase, 1, MPI_INT, 0, bridge);
		MPI_Intercomm_merge(bridge, 0, &new_uni); //create new universe
		MPI_Comm_free(universe);
		MPI_Comm_dup(new_uni, universe);
		MPI_Comm_size(*universe, &uni_size);
		MPI_Comm_rank(*universe, &uni_rank);
		*color=1;
	}
	//if all processes will be used in phase, dupe universe and mark all;
	if(phase_size == uni_size){
		MPI_Comm_dup(*universe, phase_comm);
		*color = 1;
	}
	
	//if only some processes will be used split phase_comm from universe
	if(phase_size < uni_size){
		if(uni_rank < phase_size){*color = 1;}
		MPI_Comm_split(*universe, *color, uni_rank, phase_comm);
		int test_size = 0;
		MPI_Comm_size(*phase_comm, &test_size);
	}
	printf("rank %d -- color %d -- %d\n", uni_rank, *color, uni_size);
	return 0;
}

//given a phase_comm, allocate space for local rows
int setup_grids(int** local, int** local_new, int N, MPI_Comm* phase_comm, int change){

	//if grids are setup and no change in phase, use existing grids
	if(*local != NULL && *local_new != NULL && change!=0){return 0;}

	int size;
	int rank;
	MPI_Comm_size(*phase_comm, &size);
	MPI_Comm_rank(*phase_comm, &rank);
	
	//allocate space for displacements and counts;	
	if(sendcounts!=NULL){free(sendcounts);}
	if(disp!=NULL){free(disp);}
	if(*local!=NULL){free(*local);}
	if(*local_new!=NULL){free(*local_new);}
	
	sendcounts = malloc(sizeof(int)*size);
	disp = malloc(sizeof(int)*size);

	//determine the number of rows required.
	int local_rows = N/size;
	int remand     = N%size;
	
	//the first x ranks get an additional row for balance. 
	disp[0] = 0;
	sendcounts[0] = local_rows*(N+2);
	if(remand>0){sendcounts[0] += (N+2);}
	for(int i=1; i<size; i++){
		disp[i] = disp[i-1]+sendcounts[i-1];
		sendcounts[i] = local_rows*(N+2);
		if(i < remand){sendcounts[i] += (N+2);}
	}
	
	//allocate enough space to hold expected data and ghost cells. 
	//int local_size = (sendcounts[rank]+2*(N+2)* sizeof(int));
	
	//int* temp_ptr = (int*)realloc(*local, local_size);
	int* temp_ptr = NULL;
	int* temp_ptr2 = NULL;
	temp_ptr = calloc(sendcounts[rank]+2*(N+2),sizeof(int));
	if(temp_ptr==NULL){printf("ERROR WITH REALLOC\n"); free(*local); exit(EXIT_FAILURE); }
	*local = temp_ptr;
	
	//temp_ptr = (int*)realloc(*local_new,local_size);
	temp_ptr2 = calloc(sendcounts[rank]+2*(N+2),sizeof(int));
	if(temp_ptr2==NULL){printf("ERROR WITH REALLOC\n"); free(*local_new); exit(EXIT_FAILURE); }
	*local_new= temp_ptr2;

	MPI_Scatterv(boardState+(N+2), sendcounts, disp, MPI_INT, *local+(N+2), sendcounts[rank], MPI_INT, 0, *phase_comm);
	
	return 0; 
}

int main(int argc, char *argv[])
{
    
	MPI_Init(&argc, &argv);	
	int start_time = MPI_Wtime();
	char type_of_matrix = 's';  // inital state

	int phase = 0;
	int num_phases = 1; 
	if(argc >= 5){num_phases = atoi(argv[4]);}
	int* phase_sizes = malloc(sizeof(int)*num_phases);
	for(int i=0; i<num_phases; i++){phase_sizes[i]=atoi(argv[5+i]);} 
	
	int nsteps = atoi(argv[2]);
	int N = atoi(argv[1]);         
	int phase_size;
	int color = 1; //assume everyone partipates in first stage, then trim.                       

    int local_rows;
	int global_rank;
    int uni_size;
	
	int*local = NULL; 
	int*local_new = NULL;
	
	//local timing variables
	double phase_start = 0;
	double phase_end   = 0; 
	
	//double setup_start = 0;
	double setup_time  = 0;
	double gather_time = 0; 
	
	 

	MPI_Comm universe = MPI_COMM_NULL;
	MPI_Comm phase_comm = MPI_COMM_NULL;
	MPI_Comm parent;
    MPI_Comm_get_parent(&parent); //check if this is a spawned child processes
	
	int original = uni_size;
	int previous = uni_size;
	 
	char* filename = NULL;

	filename = argv[3];
	double read_start=0, read_end=0;
	//if child process go ahead a merge into universe 
	if(parent != MPI_COMM_NULL) 
	{
		MPI_Bcast(&phase, 1, MPI_INT, MPI_PROC_NULL, parent); //Null due to being outside world. 
		MPI_Intercomm_merge(parent, 0, &universe); //merge with parent comm (current universe)	

		MPI_Comm_size(universe, &uni_size);     
		MPI_Comm_rank(universe, &global_rank);   
	}
    else //if original dup MPI_COMM_WORLD so we have a handle that we can manipulate 
	{
		MPI_Comm_dup(MPI_COMM_WORLD, &universe);
		MPI_Comm_size(universe, &uni_size);     
		MPI_Comm_rank(universe, &global_rank);    
		original = uni_size; 
		
		if(global_rank==0){
			boardState = (int*) calloc((N+2) * (N+2), sizeof(int));
			int read = -1;
			if(argc >= 3 && filename != NULL){
				read_start = MPI_Wtime();
				read = read_file(filename, N, N);
				read_end = MPI_Wtime();
			}
			if(read < 0){
				boardState = initalize_root_board(N, type_of_matrix);	
			}
		}
	}	
	
	printf("rank %d at phase loop\n", global_rank);
	for(; phase < num_phases; phase++)
	{
		phase_start = MPI_Wtime();
		color = 0;
		phase_size = phase_sizes[phase]; 
		int change = setup_comms(N, phase, phase_size, &universe, &phase_comm, &color, argv);
		printf("rank %d after comm_set -- %d\n", global_rank, color);
		if(color == 1){
			setup_grids(&local, &local_new, N, &phase_comm, change);
			printf("rank %d after grid_set\n", global_rank);
			local_rows = sendcounts[global_rank]/(N+2);
			setup_time = MPI_Wtime();
			
			//Do iterations for this phase
			printf("rank %d at iter loop\n", global_rank);
			for (int i = 0; i < nsteps; i++)
			{
				Halo(local, N, local_rows, global_rank, phase_comm);
				Step(&local, &local_new, N, local_rows);
			}
			printf("rank %d after iter loop\n", global_rank);
			phase_end = MPI_Wtime();
			//Gather data back to main board	
			MPI_Gatherv(local+(N+2), sendcounts[global_rank], MPI_INT, boardState+(N+2), sendcounts, disp, MPI_INT, 0, phase_comm);		
			gather_time = MPI_Wtime();  
			printf("rank %d after gather\n", global_rank);
			if(global_rank == 0){
				printf("%d, %d, %d, %d,",N, nsteps, previous, phase_size);
				printf("%f, %f, %f, %f, %f, ",
				        read_end-read_start,
						setup_time-phase_start,
						setup_time-start_time,
						phase_end-setup_time,
						gather_time-phase_end
			    );
			
				previous = phase_size;
				fflush(stdout);
			}	
		}
	}	
	if(global_rank ==0){
		double write_start = MPI_Wtime();
		if(filename==NULL){filename="out--2.csv";}
		write_file(filename, N, N);
		double write_end = MPI_Wtime();
		printf("%f, %f\n",write_end-write_start,write_end-start_time);
	}
	cleanup();
	if(local_new!=NULL){free(local_new);local_new=NULL;}
	if(local!=NULL){free(local);local=NULL;}
    MPI_Finalize();

    return 0;
}

// The stepping function
int*Step(int** local, int** local_new, int cols, int rows)
{
    int i, j; 
    int neighbours = 0;
	int TL, TC, TR;
	int L, R;
	int BL, BC, BR;

    // i and j are used to cycle through the pixels of a matrix, while k and l are used to look at each pixel's neighbours
	//print_matrix(*local,cols,rows);
    for (i = 1; i < rows+1; i++)
    {
        for (j = 1; j < cols+1; j++)
        {
			//get neighbor values
			TL = (*local)[(i-1)*(cols+2)+(j-1)];
			TC = (*local)[(i-1)*(cols+2)+(j+0)];
			TR = (*local)[(i-1)*(cols+2)+(j+1)];
			
			L = (*local)[(i)*(cols+2)+(j-1)];
			//C = (*local)[(i)*(cols+2)+(j+0)];
			R = (*local)[(i)*(cols+2)+(j+1)];
            
			BL = (*local)[(i+1)*(cols+2)+(j-1)];
			BC = (*local)[(i+1)*(cols+2)+(j+0)];
			BR = (*local)[(i+1)*(cols+2)+(j+1)];
			
			neighbours = TL + TC + TR + L + R + BL + BC + BR;
			
			// Stating the rules of Game of Life
 			if (neighbours == 2)
			{
				(*local_new)[i*(cols+2)+j] = (*local)[i*(cols+2)+j];
			}
			else if (neighbours == 3)
			{
				(*local_new)[i*(cols+2)+j] = 1;
			}
			else
			{
				(*local_new)[i*(cols+2)+j] = 0;
			}
			neighbours = 0; 
		}	
	}

	int* temp_grid = *local;
	*local = *local_new;
	*local_new = temp_grid;
	
    return MPI_SUCCESS;
}

/* The matrix has the size of N*N given as input, with 1's along the boundary and 0's elsewhere.
   The matrix is also zero-padded, so neighbourchecks in the step function don't go out of bounds.*/

void square(int N)
{
    for (int i = 1; i < N+1; i++)
    {
        for(int j = 1; j < N+1; j++)
        {
            if (i > 1 && i < N+1 && (j == 1 || j == N))
            {
                boardState[i*(N+2)+j] = 1;
            }
            else if (i==1 || i == N)
            {
                boardState[i*(N+2)+j] = 1;
            }
			else{boardState[i*(N+2)+j] = 0;}
        }
    }
}

// A single floater (or glider) starting in the lower left corner
void floater(int N)
{
    boardState[N+2+3] = 1;

    boardState[2*(N+2)+1] = 1;
    boardState[2*(N+2)+3] = 1;
    
    boardState[3*(N+2)+3] = 1;
    boardState[3*(N+2)+2] = 1;
}

// A glider gun in the lower left corner, shooting gliders diagonally up/right
void gosper_glider_gun(int N)
{
    boardState[N+2+25] = 1;

    boardState[2*(N+2)+23] = 1;
    boardState[2*(N+2)+25] = 1;
    
    boardState[3*(N+2)+13] = 1;
    boardState[3*(N+2)+14] = 1;
    boardState[3*(N+2)+21] = 1;
    boardState[3*(N+2)+22] = 1;
    boardState[3*(N+2)+35] = 1;
    boardState[3*(N+2)+36] = 1;

    boardState[4*(N+2)+12] = 1;
    boardState[4*(N+2)+16] = 1;
    boardState[4*(N+2)+21] = 1;
    boardState[4*(N+2)+22] = 1;
    boardState[4*(N+2)+35] = 1;
    boardState[4*(N+2)+36] = 1;

    boardState[5*(N+2)+1] = 1;
    boardState[5*(N+2)+2] = 1;
    boardState[5*(N+2)+11] = 1;
    boardState[5*(N+2)+17] = 1;
    boardState[5*(N+2)+21] = 1;
    boardState[5*(N+2)+22] = 1;

    boardState[6*(N+2)+1] = 1;
    boardState[6*(N+2)+2] = 1;
    boardState[6*(N+2)+11] = 1;
    boardState[6*(N+2)+15] = 1;
    boardState[6*(N+2)+17] = 1;
    boardState[6*(N+2)+18] = 1;
    boardState[6*(N+2)+23] = 1;
    boardState[6*(N+2)+25] = 1;

    boardState[7*(N+2)+11] = 1;
    boardState[7*(N+2)+17] = 1;
    boardState[7*(N+2)+25] = 1;
    
    boardState[8*(N+2)+12] = 1;
    boardState[8*(N+2)+16] = 1;

    boardState[9*(N+2)+13] = 1;
    boardState[9*(N+2)+14] = 1;    
}

// Grid randomly initialised
void rand_grid(int N, int fillRatio)
{
    srand(time(0));

    if(fillRatio == 0) fillRatio = rand() % 100;

    for (int i = 1; i < N+1; i++)
    {
        for(int j = 1; j < N+1; j++)
        {
            if (rand() % 100 < fillRatio % 100) 
            {
                boardState[i*(N+2)+j] = 1;
            }
        }
    }
}

// For printing the matrices in the terminal
void print_matrix(int* A, int n, int rows){
    int i;
    for (i=0; i<(rows+2)*(n+2); i++) {
        printf("%d ", A[i]);
        if ((i+1)%(n+2) == 0) printf("\n");
    }
    printf("\n");
}

void Halo(int* local, int N, int rows, int rank, MPI_Comm comm)
{
	int size;
	MPI_Comm_size(comm, &size);
	MPI_Request send_next;
	MPI_Request send_prev;
	MPI_Status status_next;
	MPI_Status status_prev;
	int next = rank+1;
	int prev = rank-1;
	if(next >= size){next = MPI_PROC_NULL;}
	if(prev < 0)   {prev = MPI_PROC_NULL;}
	
	MPI_Isend(local+(N+2)*(rows), (N+2), MPI_INT, next, 0, comm, &send_next);
	
	MPI_Recv(local, (N+2), MPI_INT, prev, 0, comm, &status_prev);

	MPI_Isend(local+(N+2), (N+2), MPI_INT, prev, 0, comm, &send_prev);
	
	MPI_Recv(local+(N+2)*(rows+1), (N+2), MPI_INT, next, 0, comm, &status_next);
}

