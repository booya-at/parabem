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
public:
    LineSegment(Vector3& v1, Vector3& v2, Vector3& v_inf);
    ~LineSegment(){};
    double gamma = 1;
    double best_gamma = 1;
    double lift_factor = 1;
    double ca();
    double b();
    Vector3* v1;
    Vector3* v2;
    Vector3 mid;
    Vector3 n;
    Vector3 t = Vector3(1,0,0);
    Vector3 ind_influence(Vector3& point, double gamma=1);
    double ind_influence(LineSegment& ls, double gamma=1);
    double lift(Vector3& v_inf, double gamma);
};

class LiftingLine{
private:
    vector<LineSegment> segments;
    vector<Vector3*> points;
public:
    LiftingLine(vector<Vector3*>);
    vector<LineSegment> get_segments();
    vector<Vector3*> get_points();
    Vector3 v_inf = Vector3(1, 0, 0);
    void solve_for_best_gamma(double cL);	 // least square method for minimizing cl/cdi for a given cL
};


#endif