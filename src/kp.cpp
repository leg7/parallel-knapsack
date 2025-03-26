#include <iostream>
#include <fstream>
#include <chrono>
#include <vector>

#include <mpi.h>
#include <stdio.h>

using namespace std;

bool verboseMode = true;

//lit une instance donnée par le nom du fichier  fileName, et crée le sac à dos correspondant dans les tableaux donnés par pointeur
void readInstance(const string& fileName, vector<int> & weights, vector<int> &values,  int & knapsackBound,  int & nbItems) {

	ifstream infile;
	infile.exceptions(ifstream::failbit | ifstream::badbit);
	infile.open(fileName.c_str());

	infile >> nbItems;

    weights.resize(nbItems);
    for (int i = 0; i < nbItems; i++) infile >> weights[i];

    values.resize(nbItems);
    for (int i = 0; i < nbItems; i++) infile >> values[i];

	infile >> knapsackBound;
	infile.close();
/*
	cout << "print knapsack: " << endl;
	cout << "maximal weight: " << knapsackBound << endl;
	for (int i = 0; i < nbItems; i++) cout << "item: " << i << " weight: " << weights[i] << " value: " << values[i] << endl;
*/

}


// resout par programmation dynamique une instance donnée par les caractéristiques du sac à dos (masses, valeurs, masse maximale)
// remplit la solution  dans le vector solution, et le cout costSolution
void solveDP(vector<int> & weights, vector<int> & values, int knapsackBound,  int  nbItems, int & costSolution, vector<bool> & solution) {



	//creation matrice programmation dynamique, initialisation avec des 0
	unsigned int** matrixDP;

	matrixDP = new unsigned int* [nbItems];
	for(int i = 0; i < nbItems; i++){
		matrixDP[i] = new unsigned int [knapsackBound+1];
		for(int j = 0; j <= knapsackBound; j++) matrixDP[i][j] =0;
	}


	//phase propagation formule de recurrence pour contruire la matrice de programmation dynamique

	for(int m = 0; m <= knapsackBound; m++) if (m <weights[0]) matrixDP[0][m] =0; else matrixDP[0][m] =values[0];

	for (int i = 1; i < nbItems; i++) {
	    for (int m = 1; m <= knapsackBound; m++) {
	             if (weights[i] <= m)
	            	 matrixDP[i][m] = max(values[i] + matrixDP[i-1][m - weights[i]],  matrixDP[i-1][m]);
	            else
	            	matrixDP[i][m] =  matrixDP[i-1][m];
	    }
	}

	//on connait alors le cout optimal
	costSolution = matrixDP[nbItems-1][knapsackBound];

	if (verboseMode) cout << "solution cost by DP: "  << costSolution << endl;
	if (verboseMode) {
		cout << "print DP matrix : " << endl;
		for(int i = 0; i < nbItems; i++){
			for(int j = 0; j <= knapsackBound; j++) cout <<  matrixDP[i][j] << " "  ;
			cout << endl;
		}
	}
//	if (verboseMode) cout << "backtrack operations:" << endl;
	solution.clear();
	solution.resize(nbItems);
	//for (int i = 0; i < nbItems; i++)  solution[i]=false; // commente car inutile, resize affecte des false par defaut

	int m = knapsackBound;
	for (int i = nbItems-1; i >=1 ; i--){
    		if (m < weights[i] || matrixDP[i][m] == matrixDP[i-1][m]) solution[i]=false;
    			else {
    				solution[i]=true;
    				m -= weights[i];
    			}
    	}


	if (m < weights[0]) solution[0]=false;
	else solution[0]=true;
//	if (verboseMode) cout << "backtrack operations achieved:" << endl;


	//destruction de la matrice de programmation dynamique
	for(int i = 0; i < nbItems; i++) delete [] matrixDP[i];
	delete [] matrixDP;

}

void printKnapsackSolution(vector<bool> & solution) {
	cout << "knapsack composition  : ";
	for (std::vector<bool>::iterator it = solution.begin() ; it != solution.end(); ++it)
		    std::cout << ' ' << *it;
	cout  << endl;
}



int main(int argc, char** argv) {
    // if (argc < 2) {
    //     cerr << "Usage: knapsack inputFile  [verbose] " << endl;
    //     cerr << "A second optional allows to disable the verbose mode for debugging" << endl;
    //     return 1;
    // }
    //
    //
    // // initialisation variables et conteneurs
    // vector<int> weights;
    // vector<int> values;
    // int knapsackBound=0;
    // int nbItems;
    // int costSolution=0;
    // vector<bool> solution;
    //
    //
    // //if (argc < 3) nbMaxItems = atoi(argv[3]);
    // if (argc < 3) verboseMode = false;
    // const char* instanceFile = argv[1];
    //
    // readInstance(instanceFile, weights, values, knapsackBound, nbItems);
    //
    // auto start = std::chrono::steady_clock::now();
    // solveDP(weights, values, knapsackBound, nbItems, costSolution, solution);
    // auto end = std::chrono::steady_clock::now();
    // std::chrono::duration<double> elapsed_seconds = end-start;
    // cout << "solution optimale trouvee de cout " << costSolution << " en temps: " << elapsed_seconds.count() << "s" << endl<< endl;
    // if (verboseMode) printKnapsackSolution(solution);
    //
	    // Initialize the MPI environment
    MPI_Init(NULL, NULL);

    // Get the number of processes
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    // Get the name of the processor
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    int name_len;
    MPI_Get_processor_name(processor_name, &name_len);

    // Print off a hello world message
    printf("Hello world from processor %s, rank %d out of %d processors\n",
           processor_name, world_rank, world_size);

    // Finalize the MPI environment.
    MPI_Finalize();

    return 0;

}
