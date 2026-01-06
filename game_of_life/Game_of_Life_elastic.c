#define _DEFAULT_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <mpi.h>
#include<time.h>

char *Step(char **local, char** local_new, int N, int rows);
void square(int N);
void floater(int N);
void gosper_glider_gun(int N);
void rand_grid(int N, int fillRatio);
void print_matrix(char* A, int n, int rows);

void Halo(char* local, int N, int rows, int rank, MPI_Comm comm);

char *boardState, *newBoardState;

void print_grid(char* local, int N, int rows){
	for(int i=0; i<rows; i++)
	{
		for(int j=0; j<N; j++)
			printf("%d ", local[i*(N)+j]);
		printf("\n");
	}
	printf("\n");
}

char* initalize_root_board(int N, char seed){
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

//update universe, phase_comm, color, and makes sure head proc is known everywhere if it changed
int setup_comms(int N, int phase, int* phase_sizes, MPI_Comm* universe, MPI_Comm* phase_comm, int* color){
	
	int uni_size = -1;
	int uni_rank = -1;
	int phase_size = phase_sizes[phase];
	
	//If phase_comm exists and is same size as previous then just use previous setup
	if(*phase_comm != MPI_COMM_NULL && phase>0 && phase_size== phase_sizes[phase-1]){return 1;}
	
	//if not figure out what needs to happen & clear room for new phase_comm
	MPI_Comm_size(*universe, &uni_size);
	MPI_Comm_rank(*universe, &uni_rank);
	if(*phase_comm != MPI_COMM_NULL){MPI_Comm_free(phase_comm);} 
	*color = 0; 
	
	//if additional processes needed, expand universe	
	if(phase_size > uni_size){
		//calculate and spawn processes as needed.
		int expand_num = phase_size - uni_size;
		MPI_Comm bridge;
		MPI_Comm new_uni;
		char gsize[10];
		char p[10];
		sprintf(gsize, "%d", N);
		sprintf(p, "%d", phase);
	    char *args[] = {gsize, p, NULL};
		MPI_Comm_spawn("./gol.exe", args, expand_num, MPI_INFO_NULL, 0, *universe, &bridge, MPI_ERRCODES_IGNORE);
		MPI_Intercomm_merge(bridge, 0, &new_uni); //create new universe
		MPI_Comm_free(universe);
		MPI_Comm_dup(new_uni, universe);
		
		MPI_Comm_size(*universe, &uni_size);
		MPI_Comm_rank(*universe, &uni_rank);
		
	}
	//printf("CHECK --- rank %d  size %d phase_size %d \n", uni_rank, uni_size, phase_size);
	//if all processes will be used in phase, dupe universe and mark all;
	if(phase_size == uni_size){
		MPI_Comm_dup(*universe, phase_comm);
		*color = 1;
	}
	
	//if only some processes will be used split phase_comm from universe
	if(phase_size < uni_size){
		if(uni_rank < phase_size){*color = 1;}
		MPI_Comm_split(*universe, *color, uni_rank, phase_comm);
	}
	
	return 0;
}

//given a phase_comm, allocate space for local rows
int setup_grids(char** localw, char** local_neww, int N, MPI_Comm phase_comm, int change){
	//assume comm_size = phase_size;
	int size;
	int rank;
	MPI_Comm_size(phase_comm, &size);
	MPI_Comm_rank(phase_comm, &rank);
	
	//if grids are setup and no change in phase, use existing grids
	if(*localw != NULL && *local_neww != NULL && !change){return 0;}
	
	char* local = *localw;
	char* local_new = *local_neww;
	
	
	
	//determine the number of rows required. 
	int rows = 0;
	int remander = N%size;
	rows = N/size;
	if(rank < remander){rows++;}
	int local_size = (N+2) * (rows+2) * sizeof(char);
	
	char* temp_ptr = (char *)realloc(local, local_size);
	if(temp_ptr==NULL){printf("ERROR WITH REALLOC\n"); free(local); exit(EXIT_FAILURE); }
	*localw = temp_ptr;
		
	temp_ptr = (char *)realloc(local_new,local_size);
	if(temp_ptr==NULL){printf("ERROR WITH REALLOC\n"); free(local_new); exit(EXIT_FAILURE); }
	*local_neww= temp_ptr;
	
	return 0; 
}

int main(int argc, char *argv[])
{
	//printf("NEW PROC!!! \n");
    
	int N = 2772;               // Evenly Divisible by 1-16, 32, 64, & 128
	char type_of_matrix = 's';  // inital state
    //int nsteps = atoi(argv[1]);          // The number of iterations per phase
	int nsteps = 500;
	int num_phases = 3;
	int phases[] = {2,3,4,5,6,7,8};
	int phase = 0; 
	int phase_size;
	int color = 0;             //color == 1 if process is participating in phase
	
    int  rows;
	int global_rank;
    int uni_size;

	char *local = NULL; 
	char *local_new = NULL;
	
	//local timing variables
	double phase_start = 0;
	double phase_time = 0; 
	
	double setup_start = 0;
	double setup_time = 0;

	double local_halo_start= 0;
	double local_halo_time = 0;

   
	double local_calc_start = 0;
	double local_calc_time  = 0;
	
	//combined timing variables
	double starttime; // endtime;
	double total_calc_time  = 0;	
	//double total_halo_time = 0;
	double min_calc_time;
    
	
    /* Initialize MPI */
    MPI_Init(&argc, &argv);       
	//MPI_Request request;
    //MPI_Status  status;
    
	MPI_Comm universe = MPI_COMM_NULL;
	MPI_Comm phase_comm = MPI_COMM_NULL;
	MPI_Comm parent;
	
	//check if this is a spawned child processes
    MPI_Comm_get_parent(&parent);
	
	//if child process go ahead a merge into universe 
	if(parent != MPI_COMM_NULL) 
	{
		MPI_Comm new_uni;
		phase = atoi(argv[2]); //use second argument to skip to the current phase in loop. 
		MPI_Intercomm_merge(parent, 0, &new_uni); //merge with parent comm (current universe)
		MPI_Comm_dup(new_uni, &universe);  //dupe new_uni with parent to maintain context with handle.
		MPI_Comm_free(&new_uni);             //cleanup unnecessary comm??
		
		//Update ID
		MPI_Comm_size(universe, &uni_size);     
		MPI_Comm_rank(universe, &global_rank);   
		//printf("SPAWNED process %d of %d\n", global_rank, uni_size); 
	}
    else //if original dup MPI_COMM_WORLD so we have a handle that we can manipulate 
	{
		MPI_Comm_dup(MPI_COMM_WORLD, &universe);
		MPI_Comm_size(universe, &uni_size);     
		MPI_Comm_rank(universe, &global_rank);    
		
		//if original root setup initial board state 
		if(global_rank==0){boardState = initalize_root_board(N, type_of_matrix);}
	}
 
    starttime = MPI_Wtime();	

	for(; phase < num_phases; phase++)
	{
		color = 0; 
		phase_size = phases[phase];
		phase_start = MPI_Wtime();
		
		//spawn processes and create phase_comm for phase;
		int change = setup_comms(N, phase, phases, &universe, &phase_comm, &color);
		
		//reallocate size for local grids
		setup_grids(&local, &local_new, N, phase_comm, change);
		
		setup_time = MPI_Wtime()-phase_start;
		
		//sanity check
		int check = -1;
		MPI_Comm_size(phase_comm,&check);
		
		if(phase_size != check && color == 0){printf("ERROR PHASE COMM MISMATCH!!!!  %d -- %d -- %d \n", phase_size, check, color);}
		
		rows = N/phase_size;
		
		if(color == 1){
			//execute phase
			for (int i = 0; i < nsteps; i++)
			{
				//printf("rank %d at calc  \n", global_rank);
				local_calc_start = MPI_Wtime();
				Step(&local, &local_new, N, rows);
				local_calc_time += MPI_Wtime() - local_calc_start;
			
				
				local_halo_start = MPI_Wtime();
				Halo(local, N, rows, global_rank, phase_comm);
				local_halo_time += MPI_Wtime()- local_halo_start;	

			}
			
			phase_time=MPI_Wtime()-phase_start;
			
			//Gather data back to main board	
			MPI_Gather(local+(N+2), (N+2)*(rows), MPI_CHAR, boardState+(N+2), (N+2)*(rows), MPI_CHAR, 0, phase_comm);
			MPI_Reduce(&local_calc_time, &total_calc_time, 1, MPI_DOUBLE, MPI_SUM, 0, phase_comm);
            MPI_Reduce(&local_calc_time, &min_calc_time, 1, MPI_DOUBLE, MPI_MIN, 0, phase_comm); 
			
			//log times 
			if(global_rank==0){
				//     id size  as  ac  mc halo, total
				printf("%d, %d, %f, %f, %f, %f, %f\n", phase, phases[phase] ,setup_time, total_calc_time/phases[phase], min_calc_time, local_halo_time/phases[phase], phase_time); fflush(stdout);
			}
		}
		//reset local variables --- TODO remove debug barrier
		local_calc_time =0;
		local_halo_time =0;
/* 		for(int j=0; j<phase_size; j++)
		{
			if(global_rank == j){
				printf("%d done with phase %d \n", global_rank, phase);
				
			}
			MPI_Barrier(phase_comm);	
		}
		//printf("-- \n");	
		MPI_Barrier(phase_comm);	 */	
	}
	MPI_Barrier(universe);	
	//endtime = MPI_Wtime() - starttime;
	//MPI_Reduce(&local_calc_time, &total_calc_time, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	//MPI_Reduce(&local_halo_time, &total_halo_time, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	
	//Have root print the results
		//if(global_rank == 0 ){
		//printf("Exeution time: %f\n", endtime);
    } 
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

// For printing the matrices in the terminal
void print_matrix(char* A, int n, int rows){
    int i;
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
    //printf("%d -- sending to %d tag %d\n", rank, next, 0);
	//printf("%d -- sending to %d  tag %d\n", rank, prev, 0); 
	
	//if(next > -1){printf("%d --- sending message to %d \n", rank, next);}
	MPI_Isend(local+(N+2)*(rows), (N+2), MPI_CHAR, next, 0, comm, &send_next);
	
	//if(prev > -1){printf("%d --- waiting on message from %d \n", rank, prev);}
	MPI_Recv(local, (N+2), MPI_CHAR, prev, 0, comm, &status_prev);
	
	//printf("%d --- message to next(%d) complete\n", rank, next);
	
	//if(prev > -1){printf("%d --- sending message to %d \n", rank, prev);}
	MPI_Isend(local+(N+2), (N+2), MPI_CHAR, prev, 0, comm, &send_prev);
	
	//if(prev > -1){printf("%d --- waiting on message from %d \n", rank, next);}
	MPI_Recv(local+(N+2)*(rows+1), (N+2), MPI_CHAR, next, 0, comm, &status_next);
}

