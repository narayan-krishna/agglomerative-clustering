/**
 * Hierarchical Agglomerative Clustering using MPI
 * Final Project for CPSC445 (High Performance Computing)
 * 
 * Krishna Narayan 
 * December 13th, 2021
 */

#include <mpi.h>
#include <iostream>
#include <string>
#include <vector>
#include <unistd.h>
#include <sstream>
#include <fstream>
#include <cmath>

using namespace std;


/*run mpi while checking errors, take an error message*/
void check_error(int status, const string message="MPI error") {
  if ( status != 0 ) {    
    cerr << "Error: " << message << endl;
    exit(1);
  }
}

/*assigns rows to process based on rank, and whether there are more processes
then points*/
void acquire_partition_rows(vector<int> &partition_rows,
                                   const int points, const int processes, 
                                   const int rank) {
  //if more processes than points
  if(processes >= points) {
    // cout << "processes greater or equal to points" << endl;
    if (rank >= points || rank == 0) {
      partition_rows.push_back(-1);    
    } else {
      partition_rows.push_back(rank);
    }

    return;
  }

  int rows_per_process = (points/processes) + ((rank < points%processes)?1:0);

  size_t index_start = rank*(points/processes) + min(rank, points%processes);
  size_t index_end = index_start + rows_per_process;

  for (int i = index_start; i < index_end; i++) {
    partition_rows.push_back(i);
  }
}

/*read string from file into a vector -> translate chars to ints*/
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

/*print a vector  of type t, with an endline if signaled*/
template <class T>
inline void print_vector(const vector<T> &vec, const int signal) {
  for(auto i : vec) {
    cout << i << " ";
  }
  if (signal == 1) cout << endl;
}

/*print a vector of type t (no signal)*/
template <class T>
inline void print_vector(const vector<T> &vec) {
  for(auto i : vec) {
    cout << i << " ";
  }
}

/*for testing -- see what value a process holds in a given variable. takes a 
  tag for readability*/
template <class T>
inline void process_has(int process, T data, string info) {
  cout << "process " << process << " has: " << data;
  cout << " " << info << endl;
}

/*for utility, get the 2d coordinates given an 1d index*/
inline void get_2d_coords(const int index, const int rows, 
                   int &x_coord, int &y_coord) {
  x_coord = (index % rows);
  y_coord = (index / rows); 
}

/*for utility, get the 1d index given 2d coordinates*/
inline int get_1d_coord(const int &x_coord, const int &y_coord, const int columns) {
  return (x_coord*columns)+y_coord;
}

/*with a vector of x coords and y coords and two pointers, get the avg of the 
  points being pointed at. place avg coords into first pointer position*/
inline void average_points(vector<float> &x, vector<float> &y, const int a, const int b) {
  x[a] = (x[a] + x[b]) / 2.0;
  y[a] = (y[a] + y[b]) / 2.0;
}

/*for utility, flatten a 2d vector into a 1d vector*/
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

/*visualize 1d distance matrix as 2d (for testing) -- visualize only bottom
  corner (distance matrix is mirrored)*/
void visualize_distance_matrix (const vector<float> &distance_matrix, 
                                const int points) {
  for(int i = 0; i < points; i ++) {
    for(int j = 0; j < i; j ++) {
      cout << distance_matrix[(i*points)+j] << " ";
    }
    cout << endl;
  }
}

/*visualize 1d distance matrix as 2d (for testing) -- signal to visualize 
  entire grid*/
void visualize_distance_matrix (const vector<float> &distance_matrix, 
                                const int points, const int signal) {
  for(int i = 0; i < points; i ++) {
    for(int j = 0; j < points; j ++) {
      cout << distance_matrix[(i*points)+j] << " ";
    }
    cout << endl;
  }
}

/*visualize point clusters*/
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

/*given a vector of rows, calculate distances at points in rows specified by
  vector. processes as many points as row count (distance marix mirror)*/
void compute_distance_matrix (const int points, 
                              const vector<int> &partition_rows, 
                              const vector<float> &x, const vector<float> &y, 
                              vector<float> &distance_matrix) {
  //iterate through the top row of matrix
  for (int i : partition_rows) {
    for (int q = 0; q < i; q++) {
        float x_diff = x[q] - x[i];
        float y_diff = y[q] - y[i];
        distance_matrix[(i*(points))+q] = sqrt((x_diff * x_diff) + 
                                              (y_diff * y_diff));
    }
  }
}

/*find the minimum distance in distance matrix, within specified partition rows*/
void compute_min_distance_between_clusters(vector<float> &min_cluster,
                                           const vector<int> &partition_rows,
                                           const vector<float> &distance_matrix, 
                                           const int points, const int rank) {
  float min_distance = -1; 
  int cluster1, cluster2;
  bool marker = false;

  for (int i : partition_rows) {
    for (int j = 0; j < i; j ++) {

      if (i > 0 && i < points) {
        float curr_distance = distance_matrix[(i*points) + j];
        if (curr_distance != 0) {
          marker = true;
          if (min_distance == -1 || 
            (min_distance != -1 && curr_distance < min_distance)) {
            cluster1 = i; 
            cluster2 = j;
            min_distance = curr_distance;
          }
        }
      }
    }
  }

  if(marker) {
    min_cluster[0] = min_distance;
    min_cluster[1] = cluster1;
    min_cluster[2] = cluster2;
  }
}

/*all processes return a minimum. the overall minimum (champion) is then found
  by rank 0 in this function*/
void extract_champion_minimum(vector<float> &cluster_candidates, int points) {
  int size = cluster_candidates.size();
  int current_min_index = -1;
  for(int i = 0; i < size; i += 3) {
    if (cluster_candidates[i] > 0) {
      if (current_min_index == -1) {
        current_min_index = i;
      } else {
        if (cluster_candidates[current_min_index] > cluster_candidates[i]) {
          current_min_index = i;
        }
      }
    }
  }

  cluster_candidates[0] = cluster_candidates[current_min_index + 1];
  cluster_candidates[1] = cluster_candidates[current_min_index + 2];
  cluster_candidates.resize(2);
}

/*update clusters based on points. (cluster 1 and cluster 2 are merged into
  cluster1, and cluster2 is erased)*/
inline void update_clusters(const int &cluster1, const int &cluster2, 
                            vector<vector<int>> &clusters) {
  for(int i : clusters[cluster2]) {
    // cout << i << endl;
    clusters[cluster1].push_back(i);
  }
  clusters.erase(clusters.begin() + cluster2);
}

int main (int argc, char *argv[]) {
  //initialize ranks and amount of processes 
  int rank;
  int p;

  // Initialized MPI
  check_error(MPI_Init(&argc, &argv), "unable to initialize MPI");
  check_error(MPI_Comm_size(MPI_COMM_WORLD, &p), "unable to obtain p");
  check_error(MPI_Comm_rank(MPI_COMM_WORLD, &rank), "unable to obtain rank");
  cout << "Starting process " << rank << "/" << p << "\n";

  //info necessary to perform task on separate processes
  int starting_points, points;
  vector<int>partition_rows;
  vector<float> x;
  vector<float> y;
  vector<float> distance_matrix;
  // vector<float> partition;
  vector<vector<int>> clusters;
  vector<float> cluster_candidates;
  vector<float> min_cluster;
  int expected_cluster_count = 4;


  /*rank 0 has things to do intitally, and then in the clustering loop*/
  if(rank == 0) {
    //read csv
    read_csv_coords(x, y, "input.csv");
    starting_points = x.size();
    clusters.resize(starting_points);

    for(int i = 0; i < starting_points; i ++) {
      clusters[i].push_back(i);
    }
    visualize_clusters(clusters);
  }

  /*start clustering loop*/
  for (int i = 0; i < 5; i ++) {

    if(rank == 0) {
      points = x.size(); //acquire points
      cluster_candidates.resize(3*p); //resize cluster candidate vector
      if (clusters.size() == 1) {
        cerr << "only one cluster left..." << endl;
        exit(1);
      }
    }

    /*broadcast point count*/
    check_error(MPI_Bcast(&points, 1, MPI_INT, 0, MPI_COMM_WORLD));  

    if(rank != 0) {
      x.resize(points);
      y.resize(points);
    }

    /*broadcast point information*/
    check_error(MPI_Bcast(&x[0], points, MPI_INT, 0, MPI_COMM_WORLD));  
    check_error(MPI_Bcast(&y[0], points, MPI_INT, 0, MPI_COMM_WORLD));  

    distance_matrix.resize(pow(points, 2));
    min_cluster.resize(3, 0);

    acquire_partition_rows(partition_rows, points, p, rank);

    /*compute distance matrix for partition (if assigned)*/
    if(partition_rows[0] != -1) {
      compute_distance_matrix(points, partition_rows, x, y, distance_matrix);
    }

    /*compute distance min distance for partition (if assigned)*/
    if(partition_rows[0] != -1) {
      compute_min_distance_between_clusters(min_cluster, partition_rows, 
                                            distance_matrix, points, rank);
    } 

    /*ensure all processes have finished calculations*/
    check_error(MPI_Barrier(MPI_COMM_WORLD));

    /*gather candidates together*/
    check_error(MPI_Gather(&min_cluster[0], 3, MPI_FLOAT, &cluster_candidates[0],
                3, MPI_FLOAT, 0, MPI_COMM_WORLD));

    if (rank == 0) {
      /*find the minimum among candidate mins*/
      extract_champion_minimum(cluster_candidates, points); 

      /*rearrange clusters, coordinates*/
      average_points(x, y, cluster_candidates[0], cluster_candidates[1]);
      x.erase(x.begin() + cluster_candidates[1]);
      y.erase(y.begin() + cluster_candidates[1]);

      update_clusters(cluster_candidates[0], cluster_candidates[1], clusters);

      cout << "cluster iteration " << i + 1 << ": " << endl; 
      visualize_clusters(clusters);
    }

    if(rank != 0) {
      x.clear();
      y.clear();
    }

    /*reset for next iteration*/
    distance_matrix.clear();
    cluster_candidates.clear();
    min_cluster.clear();
  }

  /*finalize and quit mpi, ending processes*/
  check_error(MPI_Finalize());
  cout << "Ending process " << rank << "/" << "p\n";

  return 0;
}  
