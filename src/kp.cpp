#include <cstring>
#include <iostream>
#include <fstream>
#include <chrono>
#include <vector>
#include <cstdlib>

#include <mpi.h>
#include <stdio.h>
#include <assert.h>

using namespace std;

enum ERROR_CODE {
	ERROR_CODE_TOO_MANY_RANKS,
	ERROR_CODE_MISSING_INSTANCE_ARGUMENT,
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

using fn = void (*)(Instance& instance, int world_size, int world_rank, bool verbose_mode);

void solver1(Instance& instance, int world_size, int world_rank, bool verbose_mode); 	// Schema 1
void solver2(Instance& instance, int world_size, int world_rank, bool verbose_mode);		// Schema 2

fn solvers[] = {
	nullptr,
	solver1,
	solver2,
};



int main(int argc, char** argv) {
	MPI_Init(&argc, &argv);

	int world_size; MPI_Comm_size(MPI_COMM_WORLD, &world_size);
	int world_rank; MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

	int method = 1;
	bool verbose_mode = false;
	Instance instance;


	if (argc < 3) {
		cerr << "Usage: knapsack solver[1|2|3|4|5] inputFile  [verbose] " << endl;
		cerr << "A thrid optional allows to disable the verbose mode for debugging" << endl;
		MPI_Abort(MPI_COMM_WORLD, ERROR_CODE_MISSING_INSTANCE_ARGUMENT);
	}


	if (argc > 1) {
		method = std::atoi(argv[1]);
	}

	if (method < 1 or method > 5) {
		cerr << "Le numéro du solver doit etre entre 1 et 5" << endl;
		MPI_Abort(MPI_COMM_WORLD, ERROR_CODE_INVALID_SOLVER);
	}

	const char* instance_file = argv[2];

	if (argc == 4) {
		verbose_mode = true;
	}

	/* ----- Initialisation variables et conteneurs ----- */

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

	solvers[method](instance, world_size, world_rank, verbose_mode);

	MPI_Finalize();

	free(instance.items.values);
	free(instance.items.weights);

	return EXIT_SUCCESS;
}

// DP matrix owned by each MPI rank
// Communication: Use MPI_Allreduce to combine each partial result on each step
// Ram use : HIGH O(world_size * sizeof(matrixDP)) for each process.
void solver1(Instance& instance, int world_size, int world_rank, bool verbose_mode) {
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
		int solution_cost = matrixDP[instance.items.count-1][instance.max_weight];

		vector<bool> solution;
		solution.clear();
		solution.resize(instance.items.count);

		int m = instance.max_weight;
		for (int i = instance.items.count-1; i >=1 ; i--){
			if (m < instance.items.weights[i] || matrixDP[i][m] == matrixDP[i-1][m]) {
				solution[i] = false;
			} else {
				solution[i] = true;
				m -= instance.items.weights[i];
			}
		}

		solution[0] = m >= instance.items.weights[0];

		timer += MPI_Wtime();

		if (verbose_mode) { // Print knapsack composition
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
// At least gather all to root to build final matrix.
// Communication: Use MPI_Allgather to know each partial matrix size on each process
// 					then use MPI_Allgatherv to share result on a global tab (last matrix line made).
// Ram use : Lower O(world_size * (sizeof(matrixDP) / P)) for each process.
void solver2(Instance& instance, int world_size, int world_rank, bool verbose_mode) {

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
 	// Gather all partial averages down to all the processes
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
		int solution_cost = final_matrixDP[instance.items.count-1][instance.max_weight];

		vector<bool> solution;
		solution.clear();
		solution.resize(instance.items.count);

		int m = instance.max_weight;
		for (int i = instance.items.count-1; i >=1 ; i--){
			if (m < instance.items.weights[i] || final_matrixDP[i][m] == final_matrixDP[i-1][m]) {
				solution[i] = false;
			} else {
				solution[i] = true;
				m -= instance.items.weights[i];
			}
		}

		solution[0] = m >= instance.items.weights[0];

		timer += MPI_Wtime();

		if (verbose_mode) { // Print knapsack composition
			cout << "solution cost by DP: "  << solution_cost << endl;
			cout << "print DP matrix : " << endl;
			for (int i = 0; i < instance.items.count; i++){
				for (int j = 0; j <= instance.max_weight; j++) cout <<  final_matrixDP[i][j] << " "  ;
				cout << endl;
			}

			cout << "knapsack composition  : ";
			for (std::vector<bool>::iterator it = solution.begin() ; it != solution.end(); ++it)
				std::cout << ' ' << *it;
			cout  << endl;
		}

		cout << "solution optimale trouvee de cout " << solution_cost << " en temps: " << timer << "s" << endl<< endl;

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
