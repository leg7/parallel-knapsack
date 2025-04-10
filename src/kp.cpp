#include <cstring>
#include <iostream>
#include <fstream>
#include <vector>
#include <cstdlib>
#include <cassert>
#include <string>

#include <mpi.h>

using namespace std;

// --- Forward declarations

struct Instance;

using fn = void (*)(Instance& instance, int world_size, int world_rank, bool distribute_backtrack);

void solver1(Instance& instance, int world_size, int world_rank, bool distribute_backtrack); 	// Schema 1
void solver2(Instance& instance, int world_size, int world_rank, bool distribute_backtrack);	// Schema 2
void solver3(Instance& instance, int world_size, int world_rank, bool distribute_backtrack);	// Schema 3 with Bsend and Reicv
void solver4(Instance& instance, int world_size, int world_rank, bool distribute_backtrack);	// Schema 3 with Isend and Ireicv

void copy_last_k_elements(unsigned int* source, unsigned int* destination, int total_size, int k);
void sequential_backtrack(Instance const& instance, unsigned const** matrixDP, double& timer);
void gather_results_and_backtrack(Instance const& instance, unsigned const** matrixDP, int world_size, int local_weight, bool we_are_the_master, double& timer, bool distribute_backtrack);

// Declarations

enum Error_code {
	ERROR_CODE_TOO_MANY_RANKS,
	ERROR_CODE_BAD_ARGUMENTS,
	ERROR_CODE_INVALID_SOLVER
};

struct Instance {
	struct {
		int* weights;
		int* values;
		int count = 0;
	} items;
	int max_weight = 0;
};

enum Solver_key {
	SOLVER_1,
	SOLVER_2,
	SOLVER_3,
	SOLVER_4,

	NUMBER_OF_SOLVERS,
};

fn SOLVER_LUT[] = {
	solver1,
	solver2,
	solver3,
	solver4
};

bool VERBOSE_MODE = false;

void abort_bad_argument()
{
	cerr << "Usage: knapsack --input-file [path] --schema [1-4] --verbose --distribute-backtrack" << endl;
	MPI_Abort(MPI_COMM_WORLD, ERROR_CODE_BAD_ARGUMENTS);
}

int main(int argc, char** argv) {
	MPI_Init(&argc, &argv);

	int world_size; MPI_Comm_size(MPI_COMM_WORLD, &world_size);
	int world_rank; MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

	Solver_key chosen_solver = SOLVER_1;

	const char* instance_file = nullptr;
	bool distribute_backtrack = false;

	for (int i = 1; i < argc; ++i) {
		char const* option = argv[i];

		if (strcmp(option, "--input-file") == 0) {
			++i;

			if (i >= argc) {
				cerr << "Please provide a file to the --input-file option" << endl;
				abort_bad_argument();
			}

			char const* argument = argv[i];

			instance_file = argument;
		} else if (strcmp(option, "--schema") == 0) {
			++i;

			if (i >= argc) {
				cerr << "Please provide a schema number between 1 and 4 to the option --schema\n";
				abort_bad_argument();
			}

			char const* argument = argv[i];
			chosen_solver = (Solver_key)(std::atoi(argument) - 1);

			if (chosen_solver < SOLVER_1 || chosen_solver >= NUMBER_OF_SOLVERS) {
				cerr << "The schema number for the option --schema must be between " << SOLVER_1 + 1 << " and " << NUMBER_OF_SOLVERS << endl;
				cerr << "You provided schema number " << argument << endl;
				MPI_Abort(MPI_COMM_WORLD, ERROR_CODE_INVALID_SOLVER);
			}

		} else if (strcmp(option, "--verbose") == 0) {
			VERBOSE_MODE = true;
		} else if (strcmp(option, "--distributed-backtrack") == 0) {
			distribute_backtrack = true;
		} else {
			cerr << "Argument " << option << " not supported\n";
			abort_bad_argument();
		}
	}

	if (chosen_solver <= SOLVER_2 && distribute_backtrack) {
		cerr << "Schemas 1 and 2 don't support distributed backtracking\n";
		abort_bad_argument();
	}

	/* ----- Initialisation variables et conteneurs ----- */

	Instance instance;

	// Read and parse the instance file
	ifstream infile;
	infile.exceptions(ifstream::failbit | ifstream::badbit);
	infile.open(instance_file);

	infile >> instance.items.count;

	instance.items.weights = (int*)calloc(instance.items.count, sizeof(instance.items.weights));
	for (int i = 0; i < instance.items.count; i++) infile >> instance.items.weights[i];

	instance.items.values = (int*)calloc(instance.items.count, sizeof(instance.items.values));
	for (int i = 0; i < instance.items.count; i++) infile >> instance.items.values[i];

	infile >> instance.max_weight;
	infile.close();

	if (world_size > instance.max_weight) {
		MPI_Abort(MPI_COMM_WORLD, ERROR_CODE_TOO_MANY_RANKS);
	}

	SOLVER_LUT[chosen_solver](instance, world_size, world_rank, distribute_backtrack);

	MPI_Finalize();

	free(instance.items.values);
	free(instance.items.weights);

	return EXIT_SUCCESS;
}

// DP matrix owned by each MPI rank
// Communication: Use MPI_Allreduce to combine each partial result on each step
// Ram use : HIGH O(world_size * sizeof(matrixDP))
void solver1(Instance& instance, int world_size, int world_rank, bool distribute_backtrack)
{
	// Initialise instance data
	bool const we_are_the_master = world_rank == 0;
	int subdomain_size = instance.max_weight / world_size;
	int subdomain_start = world_rank * subdomain_size + world_rank;
	bool const we_are_the_last_rank = world_size - 1 == world_rank;
	int subdomain_end = we_are_the_last_rank ? instance.max_weight : subdomain_start + subdomain_size;


	unsigned int** matrixDP = new unsigned int* [instance.items.count];
	for (int i = 0; i < instance.items.count; i++){
		matrixDP[i] = new unsigned int [instance.max_weight+1];
		for (int j = 0; j <= instance.max_weight; j++) {
			matrixDP[i][j] = 0;
		}
	}

	/* ----- Solver ----- */
	double timer = MPI_Wtime();

	// --- Phase propagation formule de recurrence pour contruire la matrice de programmation dynamique
	for (int m = 0; m <= instance.max_weight; m++) {
		bool const objectDoesNotFit = instance.items.weights[0] > m;
		if (objectDoesNotFit) {
			matrixDP[0][m] = 0;
		} else {
			matrixDP[0][m] = instance.items.values[0];
		}
	}

	unsigned* copy = (unsigned*)calloc(instance.max_weight + 1, sizeof(unsigned));
	for (int i = 1; i < instance.items.count; i++) {
		for (int m = subdomain_start; m <= subdomain_end; m++) {
			bool const objectFits = instance.items.weights[i] <= m;
			if (objectFits) {
				matrixDP[i][m] = max(instance.items.values[i] + matrixDP[i-1][m - instance.items.weights[i]],  matrixDP[i-1][m]);
			} else {
				matrixDP[i][m] = matrixDP[i-1][m];
			}
		}

		memcpy(copy, matrixDP[i], (instance.max_weight + 1) * sizeof(unsigned));
		MPI_Allreduce(copy, matrixDP[i], instance.max_weight + 1, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);
	}
	free(copy);


	if (we_are_the_master) {
		// On connait alors le cout optimal
		sequential_backtrack(instance, (unsigned const**)matrixDP, timer);
	}
	/* ----- Cleanup ----- */

	for (int i = 0; i < instance.items.count; i++) {
		delete [] matrixDP[i];
	}
	delete [] matrixDP;
}


// Each rank will have its own part of the DP matrix
// process 0 will have part 0 to max_weight / world_size
// process 1 will have part max_weight / world_size + 1 to 2 * max_weight / world_size
// etc.. last process will have little small part, but it's nothing for large instance.
// At last gather all to root to build the final matrix and backtrack.
// Communication: Use MPI_Allgather to know each partial matrix size on each process
// 					then use MPI_Allgatherv to share result on a global tab (last matrix line made).
// Ram use : Lower O(world_size * (sizeof(matrixDP) / P) + (world_size * sizeof(row) + sizeof(matrixDP)).
void solver2(Instance& instance, int world_size, int world_rank, bool distribute_backtrack)
{

	// Initialise instance data
	bool const we_are_the_master = world_rank == 0;
	int subdomain_size = instance.max_weight / world_size;
	int subdomain_start = world_rank * subdomain_size + world_rank;
	bool const we_are_the_last_rank = world_size - 1 == world_rank;
	int subdomain_end = we_are_the_last_rank ? instance.max_weight : subdomain_start + subdomain_size;
	int local_weight = subdomain_end - subdomain_start;


	unsigned int** matrixDP = new unsigned int* [instance.items.count];
	for (int i = 0; i < instance.items.count; i++){
		matrixDP[i] = new unsigned int [local_weight + 1];
		for (int j = 0; j <= local_weight; j++) {
			matrixDP[i][j] = 0;
		}
	}

	/* ----- Solver ----- */
	double timer = MPI_Wtime();

	// --- Phase propagation formule de recurrence pour contruire la matrice partielle de programmation dynamique dans chaque processus
	for (int i = 0; i <= local_weight; i++) {
		int m = subdomain_start + i;
		bool const objectDoesNotFit = instance.items.weights[0] > m;
		if (objectDoesNotFit) {
			matrixDP[0][i] = 0;
		} else {
			matrixDP[0][i] = instance.items.values[0];
		}
	}

 	// Gather all partial last line local matrixDP down to all the processes
	int* recv_counts = (int *)malloc(world_size * sizeof(int));
	int* displs = (int *)malloc(world_size * sizeof(int));
  	assert(recv_counts != NULL);
  	assert(displs != NULL);

	int send_count = local_weight + 1;
	// All process need to know how much elements on each process
	MPI_Allgather(&send_count, 1, MPI_INT, recv_counts, 1, MPI_INT, MPI_COMM_WORLD);

	// Indicate which place each process data start on reception buffer.
	displs[0] = 0;
	for (int i = 1; i < world_size; i++) {
		displs[i] = displs[i - 1] + recv_counts[i - 1];
	}

	unsigned* tab = (unsigned*)malloc((instance.max_weight + 1) * sizeof(unsigned));
  	assert(tab != NULL);
	for (int i = 1; i < instance.items.count; i++) {
  		MPI_Allgatherv(matrixDP[i - 1], send_count, MPI_UNSIGNED,
				   		tab, recv_counts, displs, MPI_UNSIGNED, MPI_COMM_WORLD);

		for (int j = 0; j <= local_weight; j++) {
			int m = subdomain_start + j;
			bool const objectFits = instance.items.weights[i] <= m;
			if (objectFits) {
				matrixDP[i][j] = max(instance.items.values[i] + tab[m - instance.items.weights[i]],  tab[m]);
			} else {
				matrixDP[i][j] = tab[m];
			}
		}

	}
	free(tab);

	unsigned int** final_matrixDP = nullptr;
	if (we_are_the_master) {
		final_matrixDP = new unsigned int * [instance.items.count];
		for (int i = 0; i < instance.items.count; i++) {
			final_matrixDP[i] = new unsigned int [instance.max_weight + 1];
		}
	}

	for (int i = 0; i < instance.items.count; i++) {
		if (we_are_the_master) {
			MPI_Gatherv(matrixDP[i], send_count, MPI_UNSIGNED,
			   final_matrixDP[i], recv_counts, displs, MPI_UNSIGNED,
			   0, MPI_COMM_WORLD);
		} else {
			MPI_Gatherv(matrixDP[i], send_count, MPI_UNSIGNED,
			   NULL, NULL, NULL, MPI_UNSIGNED,
			   0, MPI_COMM_WORLD);
		}
	}

	if (we_are_the_master) {
		sequential_backtrack(instance, (unsigned const**)final_matrixDP, timer);

		for (int i = 0; i < instance.items.count; i++) {
			delete [] final_matrixDP[i];
		}
		delete [] final_matrixDP;
	}
	/* ----- Cleanup ----- */

	for (int i = 0; i < instance.items.count; i++) {
		delete [] matrixDP[i];
	}
	delete [] matrixDP;

	free(recv_counts);
	free(displs);
}


// Each rank will have its own part of the DP matrix
// process 0 will have part 0 to max_weight / world_size
// process 1 will have part max_weight / world_size + 1 to 2 * max_weight / world_size
// etc.. last process will have little small part, but it's nothing for large instance.
// At least gather all to root to build final matrix.
// Communication: In this schema we will send only neccessary date on a buffer by MPI_Bsend and reveive by MPI_Recv,
// 					the master will only send data, the last rand will only receive data,
// 					on the first elements (row) we will only send data,
// 					on the last element (row) we will only receive data.
// Ram use : Lower O(world_size * (sizeof(matrixDP) / P) + sizeof(neccessary(row)) + sizeof(matrixDP)).
void solver3(Instance& instance, int world_size, int world_rank, bool distribute_backtrack)
{

	// Initialise instance data
	bool const we_are_the_master = world_rank == 0;
	int subdomain_size = instance.max_weight / world_size;
	int subdomain_start = world_rank * subdomain_size + world_rank;
	bool const we_are_the_last_rank = world_size - 1 == world_rank;
	int subdomain_end = we_are_the_last_rank ? instance.max_weight : subdomain_start + subdomain_size;
	int local_weight = subdomain_end - subdomain_start;


	unsigned int** matrixDP = new unsigned int* [instance.items.count];
	for (int i = 0; i < instance.items.count; i++){
		matrixDP[i] = new unsigned int [local_weight + 1];
		for (int j = 0; j <= local_weight; j++) {
			matrixDP[i][j] = 0;
		}
	}

	/* ----- Solver ----- */
	double timer = MPI_Wtime();


	// --- Phase propagation formule de recurrence pour contruire la matrice partielle de programmation dynamique dans chaque processus
	for (int i = local_weight; i >=0 ; i--) {
		int m = subdomain_start + i;
		bool const objectDoesNotFit = instance.items.weights[0] > m;
		if (objectDoesNotFit) {
			matrixDP[0][i] = 0;
		} else {
			matrixDP[0][i] = instance.items.values[0];
		}
	}

 	// send and receive with buffer needed result down to next processes (masse of next element)
	unsigned* buffer = NULL;

	// row 0, process 0 to n - 1 send data
	// row 1 to n - 1, process 0 to n - 1 send data, process 1 to n receive data
	// row n, process 1 to n receive data
	for (int i = 1; i < instance.items.count; i++) {

		if (! we_are_the_last_rank) {
			int next_item_weight = instance.items.weights[i];
			int buffer_size = (next_item_weight) * sizeof(unsigned) + MPI_BSEND_OVERHEAD;
			buffer = (unsigned*)malloc(buffer_size);
			MPI_Buffer_attach(buffer, buffer_size);

			unsigned* send_data = (unsigned*)malloc(next_item_weight * sizeof(unsigned));
			copy_last_k_elements(matrixDP[i - 1], send_data, local_weight, next_item_weight);

			MPI_Bsend(send_data, next_item_weight, MPI_UNSIGNED, world_rank + 1, i, MPI_COMM_WORLD);
			free(send_data);

			MPI_Buffer_detach(&buffer, &buffer_size);
			free(buffer);
		}


		unsigned* recv_data = NULL;
		if (! we_are_the_master) {
			int item_weight = instance.items.weights[i];
			recv_data = (unsigned*)malloc(item_weight * sizeof(unsigned));

			MPI_Recv(recv_data, item_weight, MPI_UNSIGNED, world_rank - 1, i, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

		}

		for (int j = local_weight; j >= 0; j--) {
			int m = subdomain_start + j;
			bool const objectFits = instance.items.weights[i] <= m;
			if (objectFits) {
				int diff_mass = m - instance.items.weights[i];
				bool isOnCurrentProcess = diff_mass >= subdomain_start;
				if (isOnCurrentProcess) {
					// On suppose que la masse m d'un objet est toujours <= à la taille d'un subdomain
					matrixDP[i][j] = max(instance.items.values[i] + matrixDP[i-1][j - instance.items.weights[i]],  matrixDP[i-1][j]);
				} else {
					matrixDP[i][j] = max(instance.items.values[i] + recv_data[m - subdomain_start],  matrixDP[i-1][j]);
				}
			} else {
				matrixDP[i][j] = matrixDP[i-1][j];
			}
		}

		if (recv_data != NULL) {
			free(recv_data);
		}
	}

	gather_results_and_backtrack(instance, (unsigned const**)matrixDP, world_size, local_weight, we_are_the_master, timer, distribute_backtrack);

	/* ----- Cleanup ----- */

	for (int i = 0; i < instance.items.count; i++) {
		delete [] matrixDP[i];
	}
	delete [] matrixDP;
}


// Like solver3, but using Isend and Irecv
void solver4(Instance& instance, int world_size, int world_rank, bool distribute_backtrack)
{

	// Initialise instance data
	bool const we_are_the_master = world_rank == 0;
	int subdomain_size = instance.max_weight / world_size;
	int subdomain_start = world_rank * subdomain_size + world_rank;
	bool const we_are_the_last_rank = world_size - 1 == world_rank;
	int subdomain_end = we_are_the_last_rank ? instance.max_weight : subdomain_start + subdomain_size;
	int local_weight = subdomain_end - subdomain_start;

	unsigned int** matrixDP = new unsigned int* [instance.items.count];
	for (int i = 0; i < instance.items.count; i++){
		matrixDP[i] = new unsigned int [local_weight + 1];
		for (int j = 0; j <= local_weight; j++) {
			matrixDP[i][j] = 0;
		}
	}

	/* ----- Solver ----- */
	double timer = MPI_Wtime();


	// --- Phase propagation formule de recurrence pour contruire la matrice partielle de programmation dynamique dans chaque processus
	for (int i = local_weight; i >=0 ; i--) {
		int m = subdomain_start + i;
		bool const objectDoesNotFit = instance.items.weights[0] > m;
		if (objectDoesNotFit) {
			matrixDP[0][i] = 0;
		} else {
			matrixDP[0][i] = instance.items.values[0];
		}
	}


	MPI_Request send_request = MPI_REQUEST_NULL;
	MPI_Request recv_request = MPI_REQUEST_NULL;
	unsigned* send_data = NULL;
	unsigned* recv_data = NULL;
	unsigned* next_recv_data = NULL;

	// i = 0 -> i = 1
	if (!we_are_the_last_rank) {
		int next_item_weight = instance.items.weights[1];
		send_data = (unsigned*)malloc(next_item_weight * sizeof(unsigned));
		copy_last_k_elements(matrixDP[0], send_data, local_weight, next_item_weight);

		MPI_Isend(send_data, next_item_weight, MPI_UNSIGNED, world_rank + 1, 1, MPI_COMM_WORLD, &send_request);
	}

   if (! we_are_the_master) {
        int item_weight = instance.items.weights[1];
        next_recv_data = (unsigned*)malloc(item_weight * sizeof(unsigned));
        MPI_Irecv(next_recv_data, item_weight, MPI_UNSIGNED, world_rank - 1, 1, MPI_COMM_WORLD, &recv_request);
    }

	for (int i = 1; i < instance.items.count; i++) {

		// wait receive data from previous process
		if (!we_are_the_master) {
			MPI_Wait(&recv_request, MPI_STATUS_IGNORE);
			recv_data = next_recv_data;
			next_recv_data = NULL;
		}

		// prepare receive data from previous process
		if (! we_are_the_master && i < instance.items.count - 1) {
			int item_weight = instance.items.weights[i+1];
			next_recv_data = (unsigned*)malloc(item_weight * sizeof(unsigned));

			MPI_Irecv(next_recv_data, item_weight, MPI_UNSIGNED, world_rank - 1, i + 1, MPI_COMM_WORLD, &recv_request);
		}

		// calculate DP matrix
		for (int j = local_weight; j >= 0; j--) {
			int m = subdomain_start + j;
			bool const objectFits = instance.items.weights[i] <= m;
			if (objectFits) {
				int diff_mass = m - instance.items.weights[i];
				bool isOnCurrentProcess = diff_mass >= subdomain_start;
				if (isOnCurrentProcess) {
					// On suppose que la masse m d'un objet est toujours <= à la taille d'un subdomain
					matrixDP[i][j] = max(instance.items.values[i] + matrixDP[i-1][j - instance.items.weights[i]],  matrixDP[i-1][j]);
				} else {
					matrixDP[i][j] = max(instance.items.values[i] + recv_data[m - subdomain_start],  matrixDP[i-1][j]);
				}
			} else {
				matrixDP[i][j] = matrixDP[i-1][j];
			}
		}

		if (recv_data != NULL) {
			free(recv_data);
			recv_data = NULL;
		}

		if (send_data != NULL) {
            MPI_Wait(&send_request, MPI_STATUS_IGNORE);
            free(send_data);
			send_data = NULL;
        }

		// prepare data to send to next process
		if (! we_are_the_last_rank && i < instance.items.count - 1) {
			int next_item_weight = instance.items.weights[i+1];
			send_data = (unsigned*)malloc(next_item_weight * sizeof(unsigned));
			copy_last_k_elements(matrixDP[i], send_data, local_weight, next_item_weight);

			MPI_Isend(send_data, next_item_weight, MPI_UNSIGNED, world_rank + 1, i + 1, MPI_COMM_WORLD, &send_request);
		}

	}

	if (send_data != NULL) {
        MPI_Wait(&send_request, MPI_STATUS_IGNORE);
        free(send_data);
    }

	gather_results_and_backtrack(instance, (unsigned const**)matrixDP, world_size, local_weight, we_are_the_master, timer, distribute_backtrack);
	/* ----- Cleanup ----- */

	for (int i = 0; i < instance.items.count; i++) {
		delete [] matrixDP[i];
	}
	delete [] matrixDP;
}

void copy_last_k_elements(unsigned int* source, unsigned int* destination, int total_size, int k) {
    if (k > total_size) {
        k = total_size;
    }
    //
    // memcpy(destination, source + (total_size - k), k * sizeof(unsigned int));
	std::copy(source + (total_size - k), source + total_size, destination);
}

void sequential_backtrack(Instance const& instance, unsigned const** matrixDP, double& timer)
{
	int solution_cost = matrixDP[instance.items.count-1][instance.max_weight];

	vector<bool> solution;
	solution.clear();
	solution.resize(instance.items.count);

	int m = instance.max_weight;
	for (int i = instance.items.count-1; i >= 1 ; i--){
		if (m < instance.items.weights[i] || matrixDP[i][m] == matrixDP[i-1][m]) {
			solution[i] = false;
		} else {
			solution[i] = true;
			m -= instance.items.weights[i];
		}
	}

	solution[0] = m >= instance.items.weights[0];

	timer += MPI_Wtime();

	if (VERBOSE_MODE) { // Print knapsack composition
		cout << "solution cost by DP: "  << solution_cost << endl;
		cout << "print DP matrix : " << endl;
		for (int i = 0; i < instance.items.count; i++){
			for (int j = 0; j <= instance.max_weight; j++) cout <<  matrixDP[i][j] << " "  ;
			cout << endl;
		}

		cout << "knapsack composition  : ";
		for (std::vector<bool>::iterator it = solution.begin() ; it != solution.end(); ++it)
			std::cout << ' ' << *it;

		cout  << endl;
	}

	cout << "solution optimale trouvee de cout " << solution_cost << " en temps: " << timer << "s" << endl<< endl;
}

void gather_results_and_backtrack(Instance const& instance, unsigned const** matrixDP, int world_size, int local_weight, bool we_are_the_master, double& timer, bool distribute_backtrack)
{
	// Merge each matrixDP on final_matrixDP
	int* recv_counts = (int *)malloc(world_size * sizeof(int));
	int* displs = (int *)malloc(world_size * sizeof(int));
  	assert(recv_counts != NULL);
  	assert(displs != NULL);

	int send_count = local_weight + 1;
	// All process need to know how much elements on each process
	MPI_Allgather(&send_count, 1, MPI_INT, recv_counts, 1, MPI_INT, MPI_COMM_WORLD);

	// Indicate which place each process data start on reception buffer.
	displs[0] = 0;
	for (int i = 1; i < world_size; i++) {
		displs[i] = displs[i - 1] + recv_counts[i - 1];
	}

	unsigned int ** final_matrixDP = nullptr;
	if (we_are_the_master) {
		final_matrixDP = new unsigned int * [instance.items.count];
		for (int i = 0; i < instance.items.count; i++) {
			final_matrixDP[i] = new unsigned int [instance.max_weight + 1];
		}
	}

	for (int i = 0; i < instance.items.count; i++) {
		if (we_are_the_master) {
			MPI_Gatherv(matrixDP[i], send_count, MPI_UNSIGNED,
			   final_matrixDP[i], recv_counts, displs, MPI_UNSIGNED,
			   0, MPI_COMM_WORLD);
		} else {
			MPI_Gatherv(matrixDP[i], send_count, MPI_UNSIGNED,
			   NULL, NULL, NULL, MPI_UNSIGNED,
			   0, MPI_COMM_WORLD);

		}
	}

	if (we_are_the_master) {
		if (distribute_backtrack) {
			// TODO: Distribute backtrack
			sequential_backtrack(instance, (unsigned const**)final_matrixDP, timer);
		} else {
			sequential_backtrack(instance, (unsigned const**)final_matrixDP, timer);
		}

		for (int i = 0; i < instance.items.count; i++) {
			delete [] final_matrixDP[i];
		}
		delete [] final_matrixDP;
	}

	free(recv_counts);
	free(displs);
}

