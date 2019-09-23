/**
    Purpose: some 2D Panel-methods
    @author L Lechner
*/


#ifndef case2_H
#define case2_H

#include <math.h>
#include <iostream>
#include <vector>
#include <Eigen/Core>
#include <Eigen/Dense>

#include "element_influence.h"
#include "vector.h"
#include "panel2.h"

using namespace std;


/**
 * Case2 is the base class of all Panelmethods in 2d space. it provides a mask 
 * for the different approaches. 
 */

class Case2
{
public:
    Case2(vector<Panel2*>);
    Eigen::MatrixXd matrix;
    Eigen::VectorXd rhs;
    Eigen::VectorXd result;
    vector<Panel2*> panels;
    vector<PanelWakeVertex> wake_vertices;
    Vector2 v_inf = Vector2(1, 0);

    vector<Panel2*> get_all_panels();
    vector<PanelVector2*> get_all_points();
    
    int mat_size;

    double farfield = 5;
    double A_ref = 0.;
    Vector2 mom_ref_point = Vector2(0.25, 0.);
    void solve_system();
    void sum_forces();
    
    Vector2 force;
    Vector2 center_of_pressure;
    double ca;
    double moment;
    
    //######################### To implement by subclasses ##################
    virtual void run();							
    virtual double off_body_potential(Vector2&);			
    virtual Vector2 off_body_velocity(Vector2&);				
    virtual vector<Vector2> flow_path(Vector2 start, double inc, int count);
    //#######################################################################
};


/**
 * DirichletDoublet0Case2 is a PanelelMethode for 2d objects which uses constant doublet panels.
 * The boundarycondition is applied in the boundaryintegral equation in form of an zero inner 
 * potential. -> Low-Speed Aerodynamics 11.3.2 (p. 294)
 */

class DirichletDoublet0Case2: public Case2{
protected:
    double freeflow_influence(Vector2&);
    double freeflow_influence(Panel2&);
    double doublet_influence(Vector2&, Panel2&);
    double doublet_influence(Panel2&, Panel2&);
    double wake_influence(Panel2&, Vector2&);
public:
    DirichletDoublet0Case2(vector< Panel2* >);
    virtual void run();
    virtual double off_body_potential(Vector2&);
    virtual Vector2 off_body_velocity(Vector2&);
};


/**
 * DirichletDoublet0Source0Case2 is a PanelelMethode for 2d objects which uses constant doublet and source panels.
 * The boundarycondition is applied in the boundaryintegral equationin. The inner 
 * potential is chosen to be equal to the onset flow. -> Low-Speed Aerodynamics 11.3.2 (p. 294)
 */

class DirichletDoublet0Source0Case2: public DirichletDoublet0Case2{
    double src_influence(Panel2&, Panel2&);
    double src_influence(Vector2&, Panel2&);
public:
    DirichletDoublet0Source0Case2(vector<Panel2*>);
    void run();
    double off_body_potential(Vector2&);
    Vector2 off_body_velocity(Vector2&);
};


/**
 * NeumannDoublet0Case2 is a PanelelMethode for 2d objects which uses constant doublet panels and
 * a Neumann BC. -> Low-Speed Aerodynamics 11.2.2 (p. 280)
 */

class NeumannDoublet0Case2: public Case2{
private:
    double surface_influence(Panel2&, Panel2&);
    double wake_influence(Panel2&, Vector2&);
public:
    NeumannDoublet0Case2(vector<Panel2*>);
    void run();
    double off_body_potential(Vector2&);
    Vector2 off_body_velocity(Vector2&);
};

/**
 * DirichletDoublet0Case2 is a PanelelMethode for 2d objects which uses linear doublet panels.
 * The boundarycondition is directly applied in the boundaryintegral equation and the inner 
 * potential is chosen to zero. -> Low-Speed Aerodynamics 11.3.2 (p. 294)
 * The unknowns live on the edges of the panels. The boundary integral equation is applied to gauss points 
 * This results in a m >n matrix and is solved withleast square method.
 */
//Not yet ready

class DirichletDoublet1Case2: public Case2{
private:
    double freeflow_influence(Vector2& target);
    double surface_influence(Vector2& target, PanelVector2* source);
    Vector2 surface_velocity_influence(Vector2& target, PanelVector2* source);
    double wake_influence(Vector2& target, Vector2& source);
    vector<PanelVector2> collocation_points();                   //all points where the boundary integral equation should hold
public:
    vector<PointWakeVertex> linear_wake_vertices;
    DirichletDoublet1Case2(std::vector< Panel2* >);
    void run();
    double off_body_potential(Vector2&);
    Vector2 off_body_velocity(Vector2&);
};


/**
 * NeumannSource0Case2 is a PanelelMethode for 2d objects which uses constant source panels and
 * a Neumann BC. -> Low-Speed Aerodynamics 11.2.2 (p. 279) (no lift!!!)
 */

class NeumannSource0Case2: public Case2{
private:
    double surface_influence(Panel2& target, Panel2& source);
    double wake_influence(Panel2& target, Vector2& source);
    double freeflow_influence(Panel2& target);
public:
    NeumannSource0Case2(vector<Panel2*>);
    void run();
    double off_body_potential(Vector2&);
    Vector2 off_body_velocity(Vector2&);
};
#endif
