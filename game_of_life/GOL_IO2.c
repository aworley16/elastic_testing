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
int *local, *local_new; 
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
			if(offset == (x+2)*(x+2)){fclose(ptr); return 0;}
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

//expands or shrinks universe/phase_comm to match requested phase
//returns time that universe changed. 
double setup_comms(int N, int phase, int phase_size, MPI_Comm* universe, int* color, char** argv, int* change){
	//check if universe is big enough for phase. 
	int uni_size = -1;
	int uni_rank = -1;
	MPI_Comm_size(*universe, &uni_size);
	MPI_Comm_rank(*universe, &uni_rank);
	
	//if exactly enough processes exist for phase, just use current universe. 
	if(uni_size == phase_size){
		//printf("%d SKIP \n", uni_rank);
		if(local != NULL){*change = 0; *color=1;}else{*change=1;} 
		return 0;
	}
	// if too many processes exist, split and send kill flag (color == 0) to unneeded processes. 
	else if(uni_size > phase_size){
		if(uni_rank >= phase_size){*color = 0;}
		MPI_Comm phase_comm;
		MPI_Comm_split(*universe, *color, uni_rank, &phase_comm);
		MPI_Comm_free(universe);  //free old universe
		*universe = phase_comm;
		*change = 1;
		return 0; 
	}
	//if not enough processes exist, spawn up new processes and merge into universe
	else{
		MPI_Comm old_uni = *universe;
		MPI_Comm new_uni;
		MPI_Comm bridge; 
		int expand_num = phase_size - uni_size;
		double spawn_time = MPI_Wtime();
		MPI_Comm_spawn(argv[0], &argv[1], expand_num, MPI_INFO_NULL, 0, *universe, &bridge, MPI_ERRCODES_IGNORE); 
		if(uni_rank==0){MPI_Bcast(&phase, 1, MPI_INT, MPI_ROOT, bridge);} //broadcast current phase number.
		MPI_Intercomm_merge(bridge, 0, &new_uni); //create new universe
		MPI_Comm_free(&old_uni);                  //free old handle
		*universe = new_uni;              //assign new uni to uni handle -- check if scoping issue occurs. 
/* 		printf("AT NEW BCAST\n"); fflush(stdout);
		MPI_Bcast(&phase, 1, MPI_INT, 0, *universe);
		printf("AFTER NEW BCAST\n"); fflush(stdout); */
		spawn_time = MPI_Wtime()-spawn_time;
		*change = 2;
		return spawn_time; 
	}
}

//given a phase_comm, allocate space for local rows
int setup_grids(int** local, int** local_new, int N, MPI_Comm* phase_comm, int change){

	//if grids are setup and no change in phase, use existing grids
	if(*local != NULL && *local_new != NULL && change!=0){
		printf("NULL GRID change\n");fflush(stdout);
		return 0;
	}

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
	
	return 0; 
}

int main(int argc, char *argv[])
{
	MPI_Init(&argc, &argv);	
	double start_time = MPI_Wtime();
	
	//variables for selecting process ranks for phases
	int color = 1;
	int uni_size;
	int uni_rank;
	MPI_Comm_size(MPI_COMM_WORLD, &uni_size);     
	MPI_Comm_rank(MPI_COMM_WORLD, &uni_rank);
	int change = 0; //flag to see if local grid needs to be updated.
	
	//parse arguments
	int N = atoi(argv[1]);     //size of grid N x N 
	int nsteps = atoi(argv[2]);//iterations per phase
	
	//IO file names
	char* filename = NULL;
	char type_of_matrix = 's';  //initial board state no file provided. 	
	int file = 0;
	if(argc >= 4){filename = argv[3]; file=1;} //get I/O file if provided.
	
	//parse the phase diagram if phases given or go default if not provided.
	int phase = 0;
	int num_phases = 1; 
	int* phase_sizes;
	int phase_size;
	
	if(argc >= 5){
		num_phases = atoi(argv[4]);
		phase_sizes = (int*)malloc(sizeof(int)*num_phases);
		for(int i=0; i<num_phases; i++){
			phase_sizes[i]=atoi(argv[5+i]);
		} 
	}else{ //If no phases provided use entire world as single phase. 
		MPI_Comm_size(MPI_COMM_WORLD, &uni_size); 
		phase_sizes = (int*)malloc(sizeof(int));
		*phase_sizes = uni_size;
	}
	
	//local size variables
	int local_rows; 	

	//generate placeholders for various comms
	MPI_Comm universe   = MPI_COMM_NULL;
	MPI_Comm parent     = MPI_COMM_NULL;
    MPI_Comm_get_parent(&parent); //check if this is a spawned child processes

	//local timing variables
	double read_time = 0;
	double comm_time;
	double scatter_time;
	double spawn_time;
	double grid_time;
	double work_time;
	double gather_time;
	//double write_time;
	//double total_time;

	//if child process go ahead a merge into existing universe 
	if(parent != MPI_COMM_NULL) 
	{
		//printf("CHILD PROCESS\n"); fflush(stdout);
		MPI_Comm new_uni;
		MPI_Bcast(&phase, 1, MPI_INT, 0, parent); //receive current phase via bcast
		MPI_Intercomm_merge(parent, 0, &new_uni); //merge with parent comm (current universe)	;
		universe = new_uni;                       //apply universe handle to the new universe
		
		MPI_Comm_size(universe, &uni_size);      
		MPI_Comm_rank(universe, &uni_rank);
	}
    else //if original dup MPI_COMM_WORLD so we have a handle that we can manipulate 
	{
		//printf("MAIN PROCESS\n"); fflush(stdout);
		MPI_Comm_dup(MPI_COMM_WORLD, &universe); //universe starts as dupe of WORLD
		//if root setup initial board (either reading or by generating). 
		if(uni_rank==0){
			boardState = (int*) calloc((N+2) * (N+2), sizeof(int));
			int read = -1;
			if(file){
				read_time = MPI_Wtime();
				read = read_file(filename, N, N);
				read_time = MPI_Wtime()-read_time;
			}
			if(read < 0){
				boardState = initalize_root_board(N, type_of_matrix);	
			}
			printf("N, nsteps, phase_size, uni_size,read_time,spawn_time,comm_time,grid_time, scatter_time,work_time, gather_time \n");
		}
	}
	//printf("AT PHASE\n"); fflush(stdout);
	//start at phase described in bcast or 0 if initial 
	for(; phase < num_phases; phase++)
	{
		//modify universe to fit phase
		comm_time = MPI_Wtime();
		color = 1;
		change = 0; 
		phase_size = phase_sizes[phase]; 
		//printf("%d AT comm\n", uni_rank); fflush(stdout);
		spawn_time = setup_comms(N, phase, phase_size, &universe, &color, argv, &change);

		//printf("%d AFTER comm\n", uni_rank); fflush(stdout);
		if(color == 0){  //if process kill signal -- kill process
		//      		 universe should have been updated in setup_comms
			MPI_Finalize();
			return 0; 
		}
		comm_time = MPI_Wtime()-comm_time;
		
		//if participating in phase
		
		if(color == 1){
			
			//setup local grids if universe changed. 
			grid_time=MPI_Wtime();
			//printf("%d AT GRIDS = 1\n", uni_rank); fflush(stdout);
			setup_grids(&local, &local_new, N, &universe, change);
			//printf("%d AFTER GRIDS = 1\n", uni_rank); fflush(stdout);
			local_rows = sendcounts[uni_rank]/(N+2);
			grid_time=MPI_Wtime()-grid_time;
			
			scatter_time=MPI_Wtime();
			MPI_Comm_size(universe, &uni_size);
			printf("%d AT Scatter size %d -- %d, %d\n", uni_rank, uni_size, sendcounts[0], sendcounts[1]); fflush(stdout);
			MPI_Scatterv(boardState+(N+2), sendcounts, disp, MPI_INT, local+(N+2), sendcounts[uni_rank], MPI_INT, 0, universe);
			printf("%d AFTER Scatter = 1\n", uni_rank); fflush(stdout);
			scatter_time=MPI_Wtime()-scatter_time;
			
			//Do iterations for this phase
			work_time = MPI_Wtime();
			printf("%d AT Steps = 1\n", uni_rank); fflush(stdout);
			for (int i = 0; i < nsteps; i++)
			{
				Halo(local, N, local_rows, uni_rank, universe);
				Step(&local, &local_new, N, local_rows);
			}
			printf("%d AFTER Steps = 1\n", uni_rank); fflush(stdout);
			work_time = MPI_Wtime()-work_time;
			
			//Gather data back to main board at end of phase	
			gather_time = MPI_Wtime();  
			MPI_Gatherv(local+(N+2), sendcounts[uni_rank], MPI_INT, boardState+(N+2), sendcounts, disp, MPI_INT, 0, universe);		
			gather_time = MPI_Wtime()-gather_time;  
			
			//print phase statistics
			if(uni_rank == 0){
				MPI_Comm_size(universe, &uni_size); 
				printf("%d, %d, %d, %d,", N, nsteps, phase_size, uni_size);
				printf("%f, %f, %f, %f, %f, %f, %f, \n",
				        read_time,
						spawn_time,
						comm_time,
						grid_time,
						scatter_time,
						work_time,
						gather_time
			    );
				read_time=0; //reset read time for next phase. 
			}	
		}fflush(stdout);
	}
	//print write time and overall time
	if(uni_rank ==0){
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

