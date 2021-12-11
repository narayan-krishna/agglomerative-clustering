#include <mpi.h>
#include <iostream>
#include <string>
#include <vector>
#include <unistd.h>
#include <sstream>
#include <fstream>
#include <cmath>

using namespace std;

//read string from file into a vector -> translate chars to ints
void read_csv_coords(vector<float> &x, vector<float> &y, const string &path){
    ifstream input_stream (path);

    if (!input_stream.is_open()) {
      cerr << "coudn't find/open file..." << endl;
      exit(1);
    }

    bool alternate = 0;
    for(string line; getline(input_stream, line);) {
      stringstream ss(line);

      string float_string;
      while(getline(ss, float_string, ',')) {
        if (alternate == 0) {
          x.push_back( (float)atof(float_string.c_str()) );
        } else {
          y.push_back( (float)atof(float_string.c_str()) );
        }
        alternate = !alternate; 
      }
    }
}

//print a character vector
template <class T>
void print_vector(const vector<T> &vec) {
  for(auto i : vec) {
    cout << i;
  }
}

//run mpi while checking errors, take an error message
void check_error(int status, const string message="MPI error") {
  if ( status != 0 ) {    
    cerr << "Error: " << message << endl;
    exit(1);
  }
}

void get_2d_coords(const int &index, const int &rows, 
                   int &x_coord, int &y_coord) {
  x_coord = (index % rows);
  y_coord = (index / rows); 
}

template <class T>
vector<T> flatten_2d(const vector<vector<T>> &vec_2d) {
  int rows = vec_2d[0].size();
  int columns = vec_2d.size();

  vector<T> flattened(rows*columns);

  int k = 0;
  for(int i = 0; i < rows; i++) {
    for(int j = 0; j < rows; j++) {
      flattened[k] = vec_2d[i][j];
      k++;
    }
  }

  return flattened;
}

void compute_distance_matrix (const int &particle_count, const vector<float> &x, 
                              const vector<float> &y, vector<float> &distance_matrix) {
  int distance_matrix_loc = 0;
  for(int q = 0; q < particle_count; q ++) {
    for(int k = 0; k < particle_count; k ++) {
      // if (k != q) {
      float x_diff = x[q] - x[k];
      float y_diff = y[q] - y[k];
      distance_matrix[distance_matrix_loc] = sqrt((x_diff * x_diff) + 
                                                  (y_diff * y_diff));
      distance_matrix_loc ++;
      // }
    }
  } 
}

int main (int argc, char *argv[]) {
  //initialize ranks and amount of processes 
  int rank;
  int p;

  // Initialized MPI
  check_error(MPI_Init(&argc, &argv), "unable to initialize MPI");
  check_error(MPI_Comm_size(MPI_COMM_WORLD, &p), "unable to obtain p");
  check_error(MPI_Comm_rank(MPI_COMM_WORLD, &rank), "unable to obtain rank");
  cout << "Starting process " << rank << "/" << "p\n";

  //info necessary to perform task on separate processes
  vector<float> x;
  vector<float> y;
  int points;

  if(rank == 0) {
    //read csv
    read_csv_coords(x, y, "input.csv");
    int points = x.size();
    vector<int> distances(pow(points, 2));
    // compute_distance_matrix(
  }

  //compute distance matrix

  //broadcast cut_size so that processes can resize to hold enough data
  //check_error(MPI_Bcast(&cut_size, 1, MPI_INT, 0, MPI_COMM_WORLD));  
  //cut.resize(cut_size);

  //scatter input string
  //check_error(MPI_Scatter(&sequence[0], cut_size, MPI_CHAR, &cut[0], cut_size, 
                           //MPI_CHAR, 0, MPI_COMM_WORLD));  

  //count cut sequences and add to each result array
  //count_sequence(results, cut);

  //sum results using mpi reduce to an array final results
  //check_error(MPI_Reduce(&results[0], &final_results[0], 4, MPI_INT, MPI_SUM, 
              //0, MPI_COMM_WORLD));

  //print results from rank 0
  if (rank == 0) {
    print_vector(x); cout << endl;
    print_vector(y); cout << endl;
  }

  //finalize and quit mpi, ending processes
  check_error(MPI_Finalize());
  cout << "Ending process " << rank << "/" << "p\n";

  return 0;
}  
