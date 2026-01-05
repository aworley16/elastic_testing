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
int setup_comms(int* head_proc, int phase_size, int* phase, MPI_Comm universe, MPI_Comm phase_comm, int* color){
	int old_phase_size = -1;
	int old_uni_rank = -1;
	int uni_size = -1;
	int uni_rank = -1;

	MPI_Comm_size(phase_comm, &old_phase_size);
	MPI_Comm_rank(universe, &old_uni_rank);
	MPI_Comm_size(universe, &uni_size);
	MPI_Comm new_uni;
	
	
	//if phase is same size then just use pre-existing comms
	if(phase_size == old_phase_size){
		printf("NO CHANGE --%d sees %d\n", old_uni_rank, phase_size); 
		*color= 1; 
		MPI_Comm_split(universe, *color, uni_rank, &phase_comm);
		return 0;
	}
	
	//if needed spawn additional processes expand universe	
	else if(phase_size > uni_size){
		int expand_num = phase_size - uni_size;
		MPI_Comm bridge;
		
		printf("%d at spawning -- new to spawn %d \n", old_uni_rank, expand_num);
		
		MPI_Comm_spawn("./gol.exe", MPI_ARGV_NULL, expand_num, MPI_INFO_NULL, 0, MPI_COMM_WORLD, &bridge, MPI_ERRCODES_IGNORE);
		if(old_uni_rank == 0){
			printf("sending %d to spawn\n", *phase);
			MPI_Send(phase, 1, MPI_INT, 0, 1, bridge);
		}
	
		MPI_Intercomm_merge(bridge, 0, &new_uni);
		
		//printf(" %d after merge \n", old_uni_rank);
		//fflush(stdout);
		//make sure everyone knows who root is. should be 0, but being paranoid here.
		MPI_Comm_rank(new_uni, &uni_rank);
		
		//printf("old rank %d is now rank %d \n", old_uni_rank, uni_rank);
		//fflush(stdout);
		//if(*head_proc == old_uni_rank){*head_proc = uni_rank;}
		//MPI_Bcast(head_proc, 1, MPI_INT, *head_proc, new_uni); //broadcast so everyone knows who root is. 
		//MPI_Bcast(phase, 1, MPI_INT, *head_proc, new_uni);     //broadcast so that the newbies can skip ahead to the correct phase;
		universe = new_uni;
		MPI_Comm_split(universe, *color, uni_rank, &phase_comm);
	}
 	
	//determine what processes will be active this phase.	
	MPI_Comm_size(universe, &uni_size);
	MPI_Comm_rank(universe, &uni_rank);
	if(uni_rank > phase_size){*color = 1;}

	//delete old phase_comm and create new phase_comm
	//MPI_Comm_free(&phase_comm)
	
	
	MPI_Comm_size(phase_comm, &old_phase_size);
	printf("%d sees phase_comm of size %d \n", old_uni_rank, old_phase_size);
	fflush(stdout);
	return 0;
}

//given a phase_comm, allocate space for local rows
int setup_grids(char** localw, char** local_neww, int N, MPI_Comm phase_comm){
	char* local = *localw;
	char* local_new = *local_neww;
	
	//assume comm_size = phase_size;
	int size;
	int rank;
	MPI_Comm_size(phase_comm, &size);
	MPI_Comm_rank(phase_comm, &rank);
	
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
	printf("NEW PROC!!! \n");
    
	int N = 60;//27720;               // Evenly Divisible by 1-16, 32, 64, & 128
	char type_of_matrix = 's';  // inital state
    //int nsteps = atoi(argv[1]);          // The number of iterations per phase
	int nsteps = 3;
	int num_phases = 2;
	int phases[] = {2,3,4,5,6,7,8};
	int phase = 0; 
	int phase_size;
	int color = 0;             //color == 1 if process is participating in phase
	
    int  rows;
	int global_rank;
    int size, rank;
	int head_proc = 0;   //check to see if root changes during comm merges 
    
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
	double starttime, endtime;
	double total_calc_time  = 0;	
	double total_halo_time = 0;
	double min_calc_time;
    
	
    /* Initialize MPI */
    MPI_Init(&argc, &argv);       
	MPI_Request request;
    MPI_Status  status;
    
	MPI_Comm universe;
	MPI_Comm phase_comm;
	MPI_Comm parent;
	MPI_Comm_dup(MPI_COMM_WORLD, &universe);
	
	
	/* Global IDs */
    MPI_Comm_size(universe, &size);     
    MPI_Comm_rank(universe, &global_rank);    
	MPI_Comm_dup(universe, &phase_comm);
	//check if this is a spawned child processes
    MPI_Comm_get_parent(&parent);
	
	//if it is the original world then have root setup initial grid
    if(parent == MPI_COMM_NULL && global_rank==0){
		boardState = initalize_root_board(N, type_of_matrix);
	}
	else if(parent != MPI_COMM_NULL) //else if spawned sync with parent and help recreate universal comm
	{
		printf("Spawned process at recv  -- world size %d \n", size);
		MPI_Recv(&phase, 1, MPI_INT, 0, 1, parent, &status); 
		//printf("spawned starting at phase %d \n", phase);
		MPI_Intercomm_merge(parent, 0, &universe);
		
		MPI_Comm_size(universe, &size);  
		MPI_Comm_rank(universe, &global_rank);  
		
		printf("Spawned new process rank -- %d \n", global_rank);
		phase_comm = universe;
	}
	
	
    starttime = MPI_Wtime();	

	for(; phase < num_phases; phase++)
	{
		printf("rank %d starting phase %d \n", global_rank, phase);
		fflush(stdout);
		phase_size = phases[phase];
		phase_start=MPI_Wtime();
		setup_start = MPI_Wtime();
		setup_comms(&head_proc, phase_size, &phase, universe, phase_comm, &color);
		printf("rank %d after comm setup \n", global_rank);
		fflush(stdout);
		
		int phase_check;
		MPI_Comm_size(phase_comm, &phase_check);
		
		printf("comm_split %d of %d \n", global_rank, phase_check);
		fflush(stdout);
		
		setup_grids(&local, &local_new, N, phase_comm);
		
		setup_time = MPI_Wtime()-setup_start;
		
		MPI_Comm_size(phase_comm,&phase_size);
		rows = N/phase_size;
		
		printf("rank %d at phase branch -- color %d \n", global_rank, color);
		fflush(stdout);
		if(color == 1){
			//execute phase
			for (int i = 0; i < nsteps; i++)
			{
				printf("rank %d at calc  \n", global_rank);
				fflush(stdout);
				local_calc_start = MPI_Wtime();
				Step(&local, &local_new, N, rows);
				local_calc_time += MPI_Wtime() - local_calc_start;
				
				printf("rank %d at Halo \n", global_rank);
				fflush(stdout);
				
				local_halo_start = MPI_Wtime();
				Halo(local, N, rows, global_rank, phase_comm);
				local_halo_time += MPI_Wtime()- local_halo_start;				
			}
			
			phase_time=MPI_Wtime()-phase_start;
			
			printf("%d AFTER HALO \n", global_rank);
			//Gather data back to main board	
			MPI_Gather(local+(N+2), (N+2)*(rows), MPI_CHAR, boardState+(N+2), (N+2)*(rows), MPI_CHAR, 0, phase_comm);
			MPI_Reduce(&local_calc_time, &total_calc_time, 1, MPI_DOUBLE, MPI_SUM, 0, phase_comm);
            MPI_Reduce(&local_calc_time, &min_calc_time, 1, MPI_DOUBLE, MPI_MIN, 0, phase_comm); 
			
			//log times 
			if(global_rank==0){
				//     id size  as  ac  mc halo, total
				printf("%d, %d, %f, %f, %f, %f, %f \n\n\n", phase, phases[phase] ,setup_time, total_calc_time/phases[phase], min_calc_time, local_halo_time/phases[phase], phase_time);
			}
		}
		//reset local variables --- TODO remove debug barrier
		local_calc_time =0;
		local_halo_time =0;
		for(int j=0; j<phase_size; j++)
		{
			if(global_rank == j){
				printf("%d done with phase %d \n", global_rank, phase);
				fflush(stdout);
			}
			MPI_Barrier(universe);	
		}
		printf("\n\n\n");		
	}
	endtime = MPI_Wtime() - starttime;
	//MPI_Reduce(&local_calc_time, &total_calc_time, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	//MPI_Reduce(&local_halo_time, &total_halo_time, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	
	//Have root print the results
		if(global_rank == 0 ){
		printf("Exeution time: %f\n", endtime);
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

