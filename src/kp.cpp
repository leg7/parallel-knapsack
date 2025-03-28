#include <cstring>
#include <iostream>
#include <fstream>
#include <chrono>
#include <vector>
#include <cstdlib>

#include <mpi.h>
#include <stdio.h>

using namespace std;

enum ERROR_CODE {
	ERROR_CODE_TOO_MANY_RANKS,
	ERROR_CODE_MISSING_INSTANCE_ARGUMENT,
};

struct Instance {
	struct {
		int* weights;
		int* values;
		int count = 0;
	} items;
	int max_weight = 0;
};

int main(int argc, char** argv) {
	MPI_Init(&argc, &argv);

	int world_size; MPI_Comm_size(MPI_COMM_WORLD, &world_size);
	int world_rank; MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

	bool verbose_mode = false;

	Instance instance;

	// Initialise instance data
	bool const we_are_the_master = world_rank == 0;

	if (argc == 3) {
		verbose_mode = true;
	} else if (argc < 2) {
		cerr << "Usage: knapsack inputFile  [verbose] " << endl;
		cerr << "A second optional allows to disable the verbose mode for debugging" << endl;
		MPI_Abort(MPI_COMM_WORLD, ERROR_CODE_MISSING_INSTANCE_ARGUMENT);
	}

	const char* instance_file = argv[1];

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

	// Annoying that we have to do so much work before knowing if there are too many workers
	// Maybe don't crash and just print a warning and don't use all the workers?
	if (world_size > instance.max_weight) {
		MPI_Abort(MPI_COMM_WORLD, ERROR_CODE_TOO_MANY_RANKS);
	}

	int subdomain_size = instance.max_weight / world_size;
	int subdomain_start = world_rank * subdomain_size + world_rank;
	bool const we_are_the_last_rank = world_size - 1 == world_rank;
	int subdomain_end = we_are_the_last_rank ? instance.max_weight : subdomain_start + subdomain_size;

	// DP matrix owned by each MPI rank
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

	MPI_Finalize();

	free(instance.items.values);
	free(instance.items.weights);

	return EXIT_SUCCESS;
}
