#include <iostream>
#include <fstream>
#include <chrono>
#include <vector>

#include <mpi.h>
#include <stdio.h>

using namespace std;

int main(int argc, char** argv) {
	bool verboseMode = false;
	if (argc == 3) {
		verboseMode = true;
	} else if (argc < 2) {
		cerr << "Usage: knapsack inputFile  [verbose] " << endl;
		cerr << "A second optional allows to disable the verbose mode for debugging" << endl;
		return 1;
	}

	const char* instanceFile = argv[1];

	/* ----- Initialisation variables et conteneurs ----- */
	vector<int> weights;
	vector<int> values;
	int knapsackBound=0;
	int nbItems;
	int costSolution=0;
	vector<bool> solution;

	// Read and parse the instance file
	ifstream infile;
	infile.exceptions(ifstream::failbit | ifstream::badbit);
	infile.open(instanceFile);

	infile >> nbItems;

	weights.resize(nbItems);
	for (int i = 0; i < nbItems; i++) infile >> weights[i];

	values.resize(nbItems);
	for (int i = 0; i < nbItems; i++) infile >> values[i];

	infile >> knapsackBound;
	infile.close();

	/* ----- Solver ----- */
	auto start = std::chrono::steady_clock::now();

	unsigned int** matrixDP = new unsigned int* [nbItems];
	for (int i = 0; i < nbItems; i++){
		matrixDP[i] = new unsigned int [knapsackBound+1];
		for (int j = 0; j <= knapsackBound; j++) {
			matrixDP[i][j] = 0;
		}
	}

	// Phase propagation formule de recurrence pour contruire la matrice de programmation dynamique
	for (int m = 0; m <= knapsackBound; m++) {
		if (m <weights[0]) {
			matrixDP[0][m] = 0;
		} else {
			matrixDP[0][m] = values[0];
		}
	}

	for (int i = 1; i < nbItems; i++) {
		for (int m = 1; m <= knapsackBound; m++) {
			if (weights[i] <= m) {
				matrixDP[i][m] = max(values[i] + matrixDP[i-1][m - weights[i]],  matrixDP[i-1][m]);
			} else {
				matrixDP[i][m] =  matrixDP[i-1][m];
			}
		}
	}

	// On connait alors le cout optimal
	costSolution = matrixDP[nbItems-1][knapsackBound];

	if (verboseMode) {
		cout << "solution cost by DP: "  << costSolution << endl;
		cout << "print DP matrix : " << endl;
		for (int i = 0; i < nbItems; i++){
			for (int j = 0; j <= knapsackBound; j++) cout <<  matrixDP[i][j] << " "  ;
			cout << endl;
		}
	}

	solution.clear();
	solution.resize(nbItems);

	int m = knapsackBound;
	for (int i = nbItems-1; i >=1 ; i--){
		if (m < weights[i] || matrixDP[i][m] == matrixDP[i-1][m]) {
			solution[i] = false;
		} else {
			solution[i] = true;
			m -= weights[i];
		}
	}

	solution[0] = m >= weights[0];

	auto end = std::chrono::steady_clock::now();
	std::chrono::duration<double> elapsed_seconds = end-start;

	cout << "solution optimale trouvee de cout " << costSolution << " en temps: " << elapsed_seconds.count() << "s" << endl<< endl;
	if (verboseMode) { // Print knapsack composition
		cout << "knapsack composition  : ";
		for (std::vector<bool>::iterator it = solution.begin() ; it != solution.end(); ++it)
			std::cout << ' ' << *it;
		cout  << endl;
	}

	// MPI_Init(NULL, NULL);
	//
	// int world_size;
	// MPI_Comm_size(MPI_COMM_WORLD, &world_size);
	// int world_rank;
	// MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
	// char processor_name[MPI_MAX_PROCESSOR_NAME];
	// int name_len;
	// MPI_Get_processor_name(processor_name, &name_len);
	//
	// printf("Hello world from processor %s, rank %d out of %d processors\n",
	//        processor_name, world_rank, world_size);
	//
	// MPI_Finalize();


	// Deallocations

	for (int i = 0; i < nbItems; i++) {
		delete [] matrixDP[i];
	}
	delete [] matrixDP;

	return 0;
}
