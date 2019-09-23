/**
    Purpose: some 3D Panel-methods
    @author L Lechner
*/


#ifndef case3_H
#define case3_H

#include <math.h>
#include <iostream>
#include <vector>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <omp.h>

#include "element_influence.h"
#include "vector.h"
#include "panel3.h"

using namespace std;

class CoordSys
{
public:
    CoordSys(Vector3 v_inf, Vector3 lift_ref);
    CoordSys();
    Vector3 lift_direction;
    Vector3 drag_direction;
    Vector3 side_lift_direction;
};

class Case3;

class AeroCoef3
{
public:
    AeroCoef3(Case3& c, Vector3 v_inf);
    Vector3 v_inf;
    Vector3 force;
    Vector3 cop;
    double cL;
    double cD;
    double cS;
    double cR;
    double cP;
    double cY;
    double alpha();
    vector<double> as_vector();
    vector<string> labels();
};

class Polar3
{
public:
    vector<AeroCoef3> values;
    vector<vector<double>> as_matrix();
    vector<string> labels();
};

/**
 * Case3 is the base class of all Panelmethods in 3d space. it provides a mask 
 * for the different approaches. 
 */

class Case3{
private:
    vector<PanelVector3*> cut_line;

public:
    Case3(vector<Panel3*> panels, vector<PanelVector3*> trailing_edge);
    Case3(vector<Panel3*> panels);
    ~Case3();
    
    Eigen::MatrixXf matrix;
    Eigen::MatrixXf rhs;
    Eigen::MatrixXf result;
    
    vector<Panel3*> panels;                                // containing all panels make private + getter
    vector<PanelVector3*> trailing_edge;                   // make private + getter
    vector<PanelVector3*> vertices;                        // containing only the non symmetric vertices make private + getter 
    vector<SymmetricPanel3*> sym_panels;                   // containing all symmetric panels
    vector<WakePanel3*> wake_panels;                       // make private + getter 
    vector<SymmetricWakePanel3*> sym_wake_panels;          // make private + getter 
    vector<WakePanel3*> first_wake_row;                    // delete + getter  from wake_panels
    vector<vector<PanelVector3*>> wake_streams;            // delete + getter from wake_panels
    
    vector<Panel3*> get_all_panels();
    vector<WakePanel3*> get_all_wake_panels();
    vector<PanelVector3*> get_all_points();
    
    Vector3 v_inf = Vector3(1,0,0);
    Vector3 mom_ref_point = Vector3(0,0,0);
    Vector3 lift_ref = Vector3(0,0,1);
    Vector3 trefftz_cut_pos = Vector3(10, 0, 0);
    double farfield = 5.;
    double A_ref = 1;
    string drag_calc = "trefftz";                          // "body_sum", trefftz, trailing_edge

    // PARAMETERS FOR SYMMETRIC CASE, n is normalized!!!
    Vector3 symmetric_plane_n = Vector3(0, 1, 0);
    Vector3 symmetric_plane_p = Vector3(0, 0, 0);
    Vector3 mirror(Vector3& vec);
    
    void calc_geo();

    int mat_size;
    int vert_size = 0;
    double get_volume();
    double get_surface_area();
    double get_projected_area();
    
    void append_panel(Panel3*);
    void append_panel(WakePanel3*);
    void append_wakepoint(PanelVector3*);
    void write_to_verts();
    
    void create_wake(double length, int count){this->create_wake(length, count, this->v_inf);}
    void create_wake(double length, int count, Vector3& direction);
    void relax_wake(int iterations=1, double smoothening=1);
    
    
    void sum_forces(Vector3 v_inf_);
    vector<Edge> trefftz_cut();
    vector<Edge> get_trailing_edge();
    
    Vector3 force;
    Vector3 center_of_pressure;
    double cL;			//resulting force projected to z
    double cD;			//resulting projected to v_inf
    double cS;			//resulting force projected to the side_lift_direction
    Vector3 cM;		        //resulting Moment Roll, Pitch, Yaw
    
    Vector3 off_body_velocity_vec(Vector3& point);
    Vector3 off_body_wake_velocity_vec(Vector3& point);
    void run();
    
    
    //######################### OVERWRITE ################################
    virtual Polar3 polars(vector<Vector3*> v_inf_range);
    virtual void off_body_potential(PanelVector3&);
    virtual void off_body_velocity(PanelVector3& point);
    virtual void off_body_wake_velocity(PanelVector3& point);
    virtual vector<Vector3> flow_path(Vector3& start, double inc, int count);
    virtual vector<Vector3> body_flow_path(Panel3* start, int num_traverse);
    //####################################################################
};


/**
 * DirichletDoublet0Case3 is a PanelelMethode for 2d objects which uses constant doublet panels.
 * The boundarycondition is applied in the boundaryintegral equation in form of an zero inner 
 * potential. -> Low-Speed Aerodynamics 11.3.2 (p. 294)
 */

class DirichletDoublet0Case3: public Case3{
protected:
    void surface_influence(Panel3* pan_i, Panel3* pan_j, vector<Vector3*> v_inf_range);
    void wake_influence(Panel3* pan_i, WakePanel3* pan_j, vector<Vector3*> v_inf_range);
    double freeflow_influence(Panel3* pan_i, Vector3 v_inf_);
    double freeflow_influence(Panel3* pan_i);
    double freeflow_influence(Vector3 vec);

public:
    DirichletDoublet0Case3(vector<Panel3*> panels, vector<PanelVector3*> trailing_edge);
    DirichletDoublet0Case3(vector<Panel3*> panels);
    virtual Polar3 polars(vector<Vector3*> v_inf_range);
    virtual void off_body_potential(PanelVector3&);
    virtual void off_body_velocity(PanelVector3&);
    virtual void off_body_wake_velocity(PanelVector3& point);
};

/**
 * DirichletDoublet0Source0Case3 is a PanelelMethode for 3d objects which uses constant doublet and source panels.
 * The boundarycondition is applied in the boundaryintegral equationin. The inner 
 * potential is chosen to be equal to the onset flow. -> Low-Speed Aerodynamics 11.3.2 (p. 294)
 */

class DirichletDoublet0Source0Case3: public DirichletDoublet0Case3{
protected:
    void surface_influence(Panel3* pan_i, Panel3* pan_j, vector<Vector3*> v_inf_range);
public:
    DirichletDoublet0Source0Case3(vector<Panel3*> panels, vector<PanelVector3*> trailing_edge);
    DirichletDoublet0Source0Case3(vector<Panel3*> panels);
    virtual Polar3 polars(vector<Vector3*> v_inf_range);
    virtual void off_body_potential(PanelVector3&);
    virtual void off_body_velocity(PanelVector3&);
};

#endif //case_H
