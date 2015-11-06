/**
    Purpose: some experiments on lifting line calculations
    @author L Lechner
*/


#ifndef lifting_line_H
#define lifting_line_H


#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include <Eigen/Core>

#include "element_influence.h"
#include "vector.h"


using namespace std;

class LineSegment{
private:
public:
  LineSegment(Vector3& v1, Vector3& v2, Vector3& v_inf);
  double gamma = 1;
  double best_gamma = 1;
  double lift_factor = 1;
  double ca();
  double b();
  Vector3* v1;
  Vector3* v2;
  Vector3 py_v1(); //for the python module
  Vector3 py_v2();
  Vector3 mid;
  Vector3 n;
  Vector3 t = Vector3(1,0,0);
  Vector3 ind_influence(Vector3& point, double gamma=1);
  double ind_influence(LineSegment& ls, double gamma=1);
  double lift(Vector3& v_inf, double gamma);
};

class LiftingLine{
public:
  vector<LineSegment*> segments;
  vector<Vector3> lifting_line;
  Vector3 v_inf = Vector3(1, 0, 0);
  void append_point(Vector3 v);
  void initialize(Vector3 v_inf);
  void best_gamma(double cL);			// least square methode for minimizing cl/cdi for a given cL
};


#endif