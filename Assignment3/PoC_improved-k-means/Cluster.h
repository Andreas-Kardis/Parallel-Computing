#ifndef CLUSTER_H
#define CLUSTER_H

#include "Point.h"
#include <omp.h>

class Cluster {
  private:
    double sum_cluster;
    double x_coord;
    int size;

  public:
    Cluster() {
      sum_cluster = 0;
      size = 0;
      this->x_coord = 0;
    }

    void add_point(Point point) {
      #pragma omp atomic
      sum_cluster += point.get_x_coord();

      #pragma omp atomic
      size++;
    }

    void set_x_coord(double x_coord) {
      this->x_coord = x_coord;
    }

    double get_x_coord() {
      return this->x_coord;
    }

    double get_size() {
      return this->size;
    }

    double get_sum_cluster() {
      return this->sum_cluster;
    }

    void free_points() {
      this->size = 0;
      this->sum_cluster = 0;
    }

    bool update_coords() {
      if (this->x_coord == sum_cluster/this->size) {
        return false;
      } else {
        this->x_coord = sum_cluster/this->size;
        return true;
      }
    }
};

#endif

