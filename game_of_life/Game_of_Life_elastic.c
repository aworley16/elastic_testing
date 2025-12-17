#define _DEFAULT_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include "graphics.h"
#include <math.h>
#include <sys/time.h>
#include <mpi.h>
#include<time.h>

/*
In Game Of Life, every cell on the grid is either alive or dead, the rules are:
1. A new cell is born if it has exactly three neighbours
2. A cell dies if it has more than three neighbours
3. A cell dies if it has less than two neighbours
4. For all others cases, the cells remain unchanged

Game Of Life is deterministic as all future states depend on the initial state.
*/

/*
When testing out parameters to run with graphics, mind that smaller grids will be much quicker to compute
and it is therefore advisable to run with a small waitTime,  e.g. 0.5. 
*/

char *Step(char **local, char** local_new, int N, int rows);
void square(int N);
void floater(int N);
void gosper_glider_gun(int N);
void rand_grid(int N, int fillRatio);
void drawGraphics(int N, char * boardState, float waitTime);
void print_matrix(char* A, int n, int rows);

void Halo(char* local, int N, int rows, int rank, MPI_Comm comm);

void print_grid(char* local, int N, int rows){
	for(int i=0; i<rows; i++)
	{
		for(int j=0; j<N; j++)
			printf("%d ", local[i*(N)+j]);
		printf("\n");
	}
	printf("\n");
}

const float Color = 0;
const int windowWidth = 600, windowHeight = 600, W = 1, H = 1;
char *boardState, *newBoardState;

char* initalize_root_board(int N, char seed)
{
	boardState = (char*) calloc((N+2) * (N+2), sizeof(char));
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

int MPI_Comm_spawn(const char* command,
                   char** arguments,
                   int max_process_number,
                   MPI_Info info,
                   int root,
                   MPI_Comm intracommunicator,
                   MPI_Comm* intercommunicator,
                   int* error_codes);

MPI_Comm expand(int expand_num, MPI_Comm everyone){
	   
	   MPI_Comm universe;
	   MPI_Comm_spawn(gol.exe, MPI_ARGV_NULL, expand_num, MPI_INFO_NULL, 0, MPI_COMM_WORLD, &everyone, MPI_ERRCODES_IGNORE);
	   MPI_Intercomm_merge(everyone, 0, &universe);
	   return universe;
}


int setup_phase(char** localw, char** local_neww, int N, int old_size, int new_size,  int global_rank, MPI_Comm* phase_comm, MPI_Comm* universe){
	char* local = *localw;
	char* local_new = *local_neww;
	int color = 0;
	int rows = 0;
	
	//spawn additional processes to match requirement
	if(
	{
		universe = MPI_Comm expand(int expand_num, universe);
	}
	
	//label ranks to be used this phase, and make phase_comm;
	if(global_rank < size){color = 1;}
	MPI_Comm_split(universe, color, global_rank, phase_comm);
	 
	//calculate number of local  if participating in phase
	if(color == 1){
		rows = N/size;
		//realloc board and new_board
		char* temp_ptr = (char *)realloc(local, (N+2) * (rows+2) * sizeof(char));
	    if(temp_ptr==NULL){printf("ERROR WITH REALLOC\n"); free(local); exit(EXIT_FAILURE); }
		local = temp_ptr;
		*localw = local;
		
		temp_ptr = (char *)realloc(local_new, (N+2) * (rows+2) * sizeof(char));
	    if(temp_ptr==NULL){printf("ERROR WITH REALLOC\n"); free(local_new); exit(EXIT_FAILURE); }
		local_new = temp_ptr;
		*local_neww= local_new;
		
		//Scatter the data from rank 0;
		MPI_Scatter(boardState, (N+2)*(rows+2), MPI_CHAR, local, (N+2)*(rows+2), MPI_CHAR, 0, *phase_comm);
	}

	//printf("rank %d local size %d at %p \n", global_rank, (N+2) * (rows+2), local);
	return color;
}

int main(int argc, char *argv[])
{
    if (argc < 2)
    {
        printf("Wrong number of input arguments\n");
        printf("To run, enter 'mpirun -np [number of processes] ./gol [steps in each phase]'\n");
        return -1;
    }
    
	int N = 27720;                     // Evenly Divisible by 1-16, 32, 64, & 128
    //char type_of_matrix = argv[2][0]; // The type of initial states.
	char type_of_matrix = 's';          //
    int nsteps = atoi(argv[1]);         // The number of iterations per phase
	
	int num_phases = 3;
	int phases[] = {2,3,4,5,6,7,8};
		
	double waitTime = 0;
    int  rows;
    int size, rank;
	int phase_size;
	int global_rank;
    int i, j;
    char *local = NULL;
	char *local_new = NULL;
    double starttime, time;
	//timing variables
	double setup_start = 0;
	double setup_end = 0;
	double setup_time = 0;
	double phase_end = 0;
	
	double local_halo_start= 0;
	double local_halo_time = 0;
	double total_halo_time = 0;
   
	double local_calc_start = 0;
	double local_calc_time  = 0;
	double total_calc_time  = 0;
	double min_calc_time;
    int color = 0; 
	
    /* Initialize MPI */
    MPI_Init(&argc, &argv);       
	MPI_Request request;
    MPI_Status  status;
    
	/* Global IDs */
    MPI_Comm_size(MPI_COMM_WORLD, &size);     
    MPI_Comm_rank(MPI_COMM_WORLD, &global_rank);    
	MPI_Comm phase_comm;
	MPI_Comm universe;
	
	//check if this is a spawned child processes
    MPI_Comm_get_parent(&parent);
	//if it is the original world then have root setup initial grid
    if(parent == MPI_COMM_NULL && global_rank==0){
		boardState = initalize_root_board(N, type_of_matrix);
	}
	else //else is spawned sync with parent and help recreate universal comm
	{
		
	}	
 
    starttime = MPI_Wtime();	

	for(int phase = 0; phase < num_phases; phase++)
	{
		setup_start = MPI_Wtime();
		color = setup_phase(&local, &local_new, N, phase, phases[phase], global_rank, &phase_comm);
		setup_end = MPI_Wtime();
		
		MPI_Comm_size(&phase_size, phase_comm)
		rows = N/phases_size;
		if(color == 1){
			//execute phase
			for (i = 0; i < nsteps; i++)
			{
				local_calc_start = MPI_Wtime();
				Step(&local, &local_new, N, rows);
				local_calc_time += MPI_Wtime() - local_calc_start;
				
				local_halo_start = MPI_Wtime();
				Halo(local, N, rows, global_rank, phase_comm);
				local_halo_time += MPI_Wtime()- local_halo_start;				
			}
			
			phase_end=MPI_Wtime();
			//Gather data back to main board	
			MPI_Gather(local+(N+2), (N+2)*(rows), MPI_CHAR, boardState+(N+2), (N+2)*(rows), MPI_CHAR, 0, phase_comm);
			MPI_Reduce(&local_calc_time, &total_calc_time, 1, MPI_DOUBLE, MPI_SUM, 0, phase_comm);
            MPI_Reduce(&local_calc_time, &min_calc_time, 1, MPI_DOUBLE, MPI_MIN, 0, phase_comm); 
			
			//log times 
			if(global_rank==0){
				//     id size  as  ac  mc halo, total
				printf("%d, %d, %f, %f, %f, %f, %f \n", phase, phases[phase] ,setup_end-setup_start, total_calc_time/phases[phase], min_calc_time, local_halo_time/phases[phase], phase_end-setup_start);
			}
		}
		//reset local variables --- TODO remove debug barrier
		local_calc_time =0;
		local_halo_time =0;
		MPI_Barrier(MPI_COMM_WORLD);
	}
	time = MPI_Wtime() - starttime;
	//MPI_Reduce(&local_calc_time, &total_calc_time, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	//MPI_Reduce(&local_halo_time, &total_halo_time, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	
	//Have root print the results
	/* 	if(global_rank == 0 ){
		printf("Exeution time: %f\n", time);
    } */
    MPI_Finalize();

    return 0;
}

// The stepping function
char *Step(char**localw, char** local_neww, int N, int rows)
{
    char* local = *localw;
	char* local_new = *local_neww;
    int i, j, k, l; 
    int neighbours = 0;
 
    // i and j are used to cycle through the pixels of a matrix, while k and l are used to look at each pixel's neighbours
    for (i = 1; i < rows+1; i++)
    {
        for (j = 1; j < N+1; j++)
        {
            for (k = i - 1; k < i + 2; k++)
            {
                for (l = j - 1; l < j + 2; l++)
                {
					
                        if (!(i == k && j == l))
                        {
							//printf("offset %d \n", k*(N+2)+l);
                            neighbours += local[k*(N+2)+l];
                        }
                }
            }
			// Stating the rules of Game of Life
			if (neighbours == 2)
			{
				local_new[i*(N+2)+j] = local[i*(N+2)+j];
			}
			else if (neighbours == 3)
			{
				local_new[i*(N+2)+j] = 1;
			}
			else
			{
				local_new[i*(N+2)+j] = 0;
			}
			neighbours = 0;
		}	
		//printf("ROW %d complete \n", i);
	}
	//printf("AT SWAP\n");
	//after calculating new grid swap
	//printf("before swap local:%p local_new: %p\n", local, local_new);
	char* temp_grid = local;
	local = local_new;
	//printf("mid swap local:%p local_new: %p\n", local, local_new);
	local_new = temp_grid;
	//printf("after swap local:%p local_new: %p\n", local, local_new);
	*localw = local;
	*local_neww = local_new;
	//printf("after map local:%p local_new: %p\n", *localw, *local_neww);
	
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

// The function to draw the graphics 
void drawGraphics(int N,char * boardState, float waitTime)
{
    float oneDivN = 1 / (float)(N);
    ClearScreen();
    for (int i = 1; i < N + 2; i++)
    {
        for (int j = 1; j < N + 2; j++)
        {
            if (boardState[(N+2)*j+i] == 1) DrawRectangle((float)(i - 1) * oneDivN, (float)(j - 1) * oneDivN, W, H, oneDivN, oneDivN, Color);
        }
    }
    Refresh();
    usleep(waitTime);
}

// For printing the matrices in the terminal
void print_matrix(char* A, int n, int rows){
    int i,j;
    for (i=0; i<(rows+2)*(n+2); i++) {
        printf("%c ", A[i]);
        if ((i+1)%(n+2) == 0) printf("\n");
    }
    printf("\n");
}

void Halo(char* local, int N, int rows, int rank, MPI_Comm comm)
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
	//printf("%d -- comm_size %d \n", rank, size);
    //printf("%d -- sending %d chars to %d tag %d\n", rank, N+2, next, 0);
	//printf("%d -- sending %d chars to %d  tag %d\n", rank, N+2, prev, 0); 
	MPI_Isend(local+(N+2)*(rows), (N+2), MPI_CHAR, next, 0, comm, &send_next);
	MPI_Isend(local+(N+2), (N+2), MPI_CHAR, prev, 0, comm, &send_prev);
		
	MPI_Barrier(comm);
	//printf("%d -- waiting on %d chars from %d tag %d\n", rank, N+2, next,0);
	//printf("%d -- waiting on %d chars from %d  tag %d\n", rank, N+2, prev,0); 		
	MPI_Recv(local, (N+2), MPI_CHAR, prev, 0, comm, &status_prev);
	//printf("%d -- after prev wait\n", rank);
	
	MPI_Recv(local+(N+2)*(rows+1), (N+2), MPI_CHAR, next, 0, comm, &status_next);
	//printf("%d -- after next wait\n", rank);
	
	//MPI_Barrier(comm);
}

