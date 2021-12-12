#include <mpi.h>
#include <iostream>
#include <string>
#include <vector>
#include <unistd.h>
#include <sstream>
#include <fstream>
#include <cmath>

using namespace std;

//run mpi while checking errors, take an error message
void check_error(int status, const string message="MPI error") {
  if ( status != 0 ) {    
    cerr << "Error: " << message << endl;
    exit(1);
  }
}

void check_input_count(const int &argc) {
  if ( argc != 2 ) {    
    cerr << "Error: requires desired cluster count -> "; 
    cerr << "./a.out <CLUSTER COUNT>" << endl;
    exit(1);
  }
}


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
    cout << i << " ";
  }
}

inline void get_2d_coords(const int &index, const int &rows, 
                   int &x_coord, int &y_coord) {
  x_coord = (index % rows);
  y_coord = (index / rows); 
}

inline void average_points(vector<float> &x, vector<float> &y, int a, int b) {
  x[a] = (x[a] + x[b]) / 2.0;
  y[a] = (y[a] + y[b]) / 2.0;
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

void compute_distance_matrix (const int &points, const vector<float> &x, 
                              const vector<float> &y, 
                              vector<float> &distance_matrix) {
  int distance_matrix_loc = 0;
  for(int q = 0; q < points; q ++) {
    for(int k = 0; k < points; k ++) {
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

void visualize_distance_matrix (const vector<float> &distance_matrix, 
                                const int &points) {
  for(int i = 0; i < points; i ++) {
    for(int j = 0; j < points; j ++) {
      cout << distance_matrix[(i*(points))+j] << " ";
    }
    cout << endl;
  }
}

void visualize_clusters (const vector<vector<int>> &clusters) {
  for(auto i : clusters) {
    cout << "< ";
    for(auto j : i) {
      cout << j + 1 << " ";
    }
    cout << "> ";
  }
  cout << endl;
}

void compute_min_distance_between_clusters(int &cluster1, int &cluster2,
                                           const vector<float> &distance_matrix, 
                                           const int &points) {
  float min_distance = -1;
  for(int i = 0; i < pow(points, 2); i ++) {
    float curr_distance = distance_matrix[i];
    // cout << curr_distance << endl;
    if (curr_distance != 0) {
      if (min_distance == -1 || 
         (min_distance != -1 && curr_distance < min_distance)) {
        get_2d_coords(i, points, cluster1, cluster2);
        min_distance = curr_distance;
      }
    }
  }
  cout << cluster1 << ", " << cluster2 << " (" << min_distance << ")" << endl;
}

inline void update_clusters(const int &cluster1, const int &cluster2, 
                     vector<vector<int>> &clusters) {
  //what does update clusters do? add cluster 2 to cluster 1, get rid of cluster 2
  for(int i : clusters[cluster2]) {
    cout << i << endl;
    clusters[cluster1].push_back(i);
  }
  clusters.erase(clusters.begin() + cluster2);
}

void update_dist_matrix () {}

int main (int argc, char *argv[]) {
  //initialize ranks and amount of processes 
  int rank;
  int p;

  // Initialized MPI
  check_error(MPI_Init(&argc, &argv), "unable to initialize MPI");
  check_error(MPI_Comm_size(MPI_COMM_WORLD, &p), "unable to obtain p");
  check_error(MPI_Comm_rank(MPI_COMM_WORLD, &rank), "unable to obtain rank");
  cout << "Starting process " << rank << "/" << "p\n";

  check_input_count(argc);

  //info necessary to perform task on separate processes
  int points;
  vector<float> x;
  vector<float> y;
  vector<vector<int>> clusters;
  // int expected_cluster_count = atoi(argv[2]);
  int expected_cluster_count = 2;

  //have main keep track of clusters?

  if(rank == 0) {
    //read csv
    read_csv_coords(x, y, "input.csv");
    //i have a bunch of points ...
      //vector of x coords
      //vector of y coords
      //vector of vectors of point nums as clusters
    int points = x.size();
    int starting_points = points;
    clusters.resize(points);

    for(int i = 0; i < points; i ++) {
      clusters[i].push_back(i);
    }
    vector<float> distance_matrix(pow(points, 2));

    //i compute the distance matrix...
    for (int i = 0; i < starting_points - expected_cluster_count; i ++) {

      compute_distance_matrix(points, x, y, distance_matrix);

      cout << " ----distance---- " << endl;
      visualize_distance_matrix(distance_matrix, points);
      cout << " ----points 2---- " << endl;

      int cluster1, cluster2;
      //i find the minimum distance... 
      compute_min_distance_between_clusters(cluster1, cluster2, distance_matrix, points);
      //i add whatever is in the index of the second cluster to whatever is in the first cluster
      //i remove the 4 coordinates from the x,y vector and add back the average of the points

      print_vector(x); cout << endl;
      print_vector(y); cout << endl;
      cout << " ----points 2---- " << endl;

      average_points(x, y, cluster1, cluster2);
      x.erase(x.begin() + cluster2);
      y.erase(y.begin() + cluster2);

      print_vector(x); cout << endl;
      print_vector(y); cout << endl;

      cout << " ----clusters---- " << endl;

      visualize_clusters(clusters);
      update_clusters(cluster1, cluster2, clusters);
      cout << "adding from " << cluster2 << " to " << cluster1 << endl;
      visualize_clusters(clusters);

      cout << " ----updates---- " << endl;
      
      points = x.size();
      cout << "new point count: " << points << endl;
    }
  }

  if (rank == 0) {
    // print_vector(x); cout << endl;
    // print_vector(y); cout << endl;
  }

  //finalize and quit mpi, ending processes
  check_error(MPI_Finalize());
  cout << "Ending process " << rank << "/" << "p\n";

  return 0;
}  
