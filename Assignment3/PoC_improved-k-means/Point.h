
#ifndef POINT_H
#define POINT_H

class Point {
  private:
    double x_coord;
    int cluster_id;

  public:
    Point() {
      x_coord = 0;
      cluster_id = 0;
    }

    void set_x_coord(double x_coord) {
      this->x_coord = x_coord;
    }

    void set_cluster_id(int cluster_id) {
      this->cluster_id = cluster_id;
    }

    double get_x_coord() {
      return this->x_coord;
    }

    int get_cluster_id() {
      return this->cluster_id;
    }
};

#endif
