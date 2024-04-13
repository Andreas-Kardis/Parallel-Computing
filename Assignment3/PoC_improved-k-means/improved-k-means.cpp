#include <math.h>
#include "Point.h"
#include "Cluster.h"
#include <stdio.h>
#include <omp.h>
#include <vector>

int num_point = 231735;
int num_cluster = 20;
int max_iterations = 20;
int num_threads = 8;

std::vector<Point> init_point();

void add_point_cluster(std::vector<Point> &points, std::vector<Cluster> &clusters);
void set_centroids(std::vector<Point> &points, std::vector<Cluster> &clusters);

double euclidean_dist(Point point, Cluster cluster);
double manhattan_dist(Point point, Cluster cluster);

void compute_distance_euclidean(std::vector<Point> &points, std::vector<Cluster> &clusters);
void compute_distance_manhattan(std::vector<Point> &points, std::vector<Cluster> &clusters);

void mean_sd(std::vector<Point> &points, std::vector<Cluster> &clusters, double *point_mean, double *point_sd);
double compute_q(std::vector<Point> &points, double *point_mean, double *point_sd);

bool update_clusters(std::vector<Cluster> &clusters);

int main() {
  printf("Number of points %d\n", num_point);
  printf("Number of clusters %d\n", num_cluster);
 
  double time_point1 = omp_get_wtime();
  printf("Initialising %d points...\n", num_point);
  std::vector<Point> points;
  points = init_point(); // initialise the data points to be added to a cluster
 
  printf("Initialising %d clusters...\n", num_cluster);
  std::vector<Cluster> clusters(num_cluster);
 
  double time_point2 = omp_get_wtime();
  double duration = time_point2 - time_point1;
  printf("\nPoints and clusters generated in: %f seconds\n", duration);

  add_point_cluster(points, clusters); // add points to clusters

  printf("Calculating initial centroids...\n");
  set_centroids(points, clusters);

  bool centroids_changed = true;
  int iterations = 0;

  double mean, sd;
  mean_sd(points, clusters, &mean, &sd);

  double q = compute_q(points, &mean, &sd);

  if (q <= 8.4595)
    {
        printf("Computing with Euclidean Distance\n");
    } else {
        printf("Computing with Manhattan Distance\n");
    }

  printf("Starting main k-means loop...\n");
  while (centroids_changed && iterations < max_iterations) {
    iterations++;

    if (q <= 8.4595) {
      compute_distance_euclidean(points, clusters);
      set_centroids(points, clusters);
    } else {
      compute_distance_manhattan(points, clusters);
      set_centroids(points, clusters);
    }

    centroids_changed = update_clusters(clusters);
  }
  
  double time_point3 = omp_get_wtime();
  duration = time_point3 - time_point2;
  double total_time = time_point3 - time_point1;

  printf("Number of iterations: %d, total time: %f seconds, time per iteration: %f seconds\n", iterations, total_time, duration/iterations);
}

// Initialisation of all data points to vector<Point>
std::vector<Point> init_point() {
  std::vector<Point> points(num_point);

  FILE *file = fopen("Employee_Payroll.txt", "r");

  if (file != NULL) {
    char *line = NULL;
    size_t len = 0;
    ssize_t read;
    int i = 0;
    double read_x_coord;

    while ((read = getline(&line, &len, file)) != -1 && num_point > i) {
      sscanf(line, "%lf", &read_x_coord);
      points[i].set_x_coord(read_x_coord);
      i++;
    }

    fclose(file);
  } else {
    printf("Error opening the file\n");
  }

  return points;
}

void add_point_cluster(std::vector<Point> &points, std::vector<Cluster> &clusters) {
  int init_cluster_size = floor(num_point/num_cluster);
  int loop_cluster = init_cluster_size; // variable to loop through points in a specific cluster
  int j = 0;

  for (int i = 0; i < num_cluster; i++) {
    Cluster &cluster = clusters[i];

    for (; j < loop_cluster; j++) {
      Point &point = points[j];
      points[j].set_cluster_id(i);
      clusters[i].add_point(points[j]);
    }

    j = loop_cluster;
    loop_cluster += init_cluster_size; // move to next cluster

    if (i == num_cluster-2) { // approaching alst iteration
      loop_cluster = num_point;
    }

  }
}

// set initial centroids by finding the sum of the points' x coords and dividing it by
// the amount of points considered. Giving the mean of the points.
void set_centroids(std::vector<Point> &points, std::vector<Cluster> &clusters) {
  for (int i = 0; i < num_cluster; i++) {
    Cluster &cluster = clusters[i];
    double point_sum = 0;
    int point_counter = 0;
    double cluster_mean = 0;

    for (int j = 0; j < num_point; j++) {
      Point &point = points[j];

      if (point.get_cluster_id() == i) {
        point_sum += point.get_x_coord();
        point_counter++;
      }
    }

    cluster_mean = point_sum/point_counter;
    clusters[i].set_x_coord(cluster_mean);
  }
}

double euclidean_dist(Point point, Cluster cluster) {
  double distance = sqrt(pow(point.get_x_coord() - cluster.get_x_coord(), 2));
  return distance;
}

double manhattan_dist(Point point, Cluster cluster) {
  double distance = abs(point.get_x_coord() - cluster.get_x_coord());
  return distance;
}

void compute_distance_euclidean(std::vector<Point> &points, std::vector<Cluster> &clusters) {
  double min_distance;
  int min_index;

  #pragma omp parallel for schedule(static) private(min_distance, min_index) num_threads(num_threads)
  for (int i = 0; i < num_point; i++) {
    Point &point = points[i];
    min_distance = euclidean_dist(point, clusters[0]);
    min_index = 0;

    for (int j = 1; j < num_cluster; j++) {
      Cluster &cluster = clusters[j];

      double distance = euclidean_dist(point, cluster);

      if (distance < min_distance) {
        min_distance = distance;
        min_index = j;
      }
    }

    point.set_cluster_id(min_index);
    clusters[min_index].add_point(point);
  }
}


void compute_distance_manhattan(std::vector<Point> &points, std::vector<Cluster> &clusters) {
  double min_distance;
  int min_index;

  #pragma omp parallel for schedule(static) private(min_distance, min_index) num_threads(num_threads)
  for (int i = 0; i < num_point; i++) {
    Point &point = points[i];
    min_distance = manhattan_dist(point, clusters[0]);
    min_index = 0;

    for (int j = 1; j < num_cluster; j++) {
      Cluster &cluster = clusters[j];

      double distance = manhattan_dist(point, cluster);

      if (distance < min_distance) {
        min_distance = distance;
        min_index = j;
      }
    }

    point.set_cluster_id(min_index);
    clusters[min_index].add_point(point);
  }

}

void mean_sd(std::vector<Point> &points, std::vector<Cluster> &clusters, double *point_mean, double *point_sd) {
  double point_sum = 0;
  int point_counter = 0;
  *point_mean = 0;

  for (int i = 0; i < num_cluster; i++) {
    point_sum += clusters[i].get_sum_cluster();
    point_counter += clusters[i].get_size();
  }

  *point_mean = point_sum/point_counter;

  for (int i = 0; i < point_counter; i++) {
    Point &point = points[i];
    *point_sd += pow(point.get_x_coord() - *point_mean, 2);
  }

  *point_sd = sqrt(*point_sd/point_counter);
}

double compute_q(std::vector<Point> &points, double *point_mean, double *point_sd) {
  double q = 0;
  double numerator = 0;
  double denominator = 0;

  for (int i = 0; i < num_point; i++) {
    Point &point = points[i];

    numerator += pow(point.get_x_coord() - *point_mean, 4);
  }

  numerator = numerator/num_point;
  denominator = pow(*point_sd, 2);

  q = numerator/denominator;

  return q;
}

bool update_clusters(std::vector<Cluster> &clusters) {
  bool cluster_centroid_moved = false;

  for (int i = 0; i < num_cluster; i++) {
    cluster_centroid_moved = clusters[i].update_coords();
    clusters[i].free_points();
  }

  return cluster_centroid_moved;
}
