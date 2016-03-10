/**
    Purpose: Geometric properties of 3d-Panels
    @author L Lechner
*/


#ifndef panel3_H
#define panel3_H
#include <iostream>
#include <vector>
#include <sstream>
#include <iostream>
#include <Eigen/Dense>
#include <math.h>

#include "vector.h"
using namespace std;

// forward decleration
class SymmetricPanel3;
class Edge;


class Panel3{
private:
    void ls_gradient();
    PanelVector3 project_to_panel(Panel3&);
    int nr;
    bool _is_symmetric = false;
    double sigma = 0.;
    Vector3 velocity = Vector3(0,0,0);
    double cp = 0.;
    double potential = 0.;
    double mue = 0.;

protected:
    bool neighbours_set = false;
    
public:
    
    Panel3();
    ~Panel3(){};
    Panel3(vector<PanelVector3*>);
    virtual string type() {return "Panel3";}
    vector<PanelVector3*> points;                       //Maybe protected
    virtual vector<PanelVector3*> get_points();
    
    vector<Edge> get_edges();
    Panel3* get_edge_neighbour(Edge e);
    
    virtual bool is_symmetric();
    void set_symmetric();

    void append_point(PanelVector3*);
    virtual void calc_geo();
    void set_nr(int nr){this->nr = nr;}
    virtual int get_nr(){return this->nr;}
    
    vector<Panel3*> neighbours;
    virtual void set_neighbours();
    double area = 0.;
    double side_sum = 0.;

    Vector3 center, n, m, l;
    Vector3 loc_coord(Vector3);
    
    void set_sigma(double sigma){this->sigma = sigma;}
    void set_mue(double mue){this->mue = mue;}
    virtual double get_sigma(){return this->sigma;}
    virtual double get_mue(){return this->mue;}
    
    //  FLUID PROPERTIES GETTER & SETTER
    virtual Vector3 get_force();
    virtual Vector3 get_velocity(){return this->velocity;};
    virtual double get_cp(){return this->cp;}
    virtual double get_potential(){return this->potential;}
    void add_velocity(Vector3 vel){this->velocity += vel;}
    void reset_properties(){this->velocity.setZero(); this->cp=0; this->mue=0; this->potential=0;this->sigma=0;}
    void set_potential(double pot){this->potential = pot;}
    
    void compute_gradient(Vector3& v_inf);
    void compute_pressure(Vector3& v_inf);

    string tostring();
    int aligned(PanelVector3*, PanelVector3*);
    bool is_inside(Vector3&);
};


class WakePanel3: public virtual Panel3{
private:
    Panel3* upper_operating_panel;
    Panel3* lower_operating_panel;
public:
    virtual bool is_symmetric();
    virtual string type() {return "WakePanel3";}
    virtual Panel3* get_lower_operating_panel();
    virtual Panel3* get_upper_operating_panel();
    void set_lower_operating_panel(Panel3*);
    void set_upper_operating_panel(Panel3*);
};


class SymmetricPanel3: public virtual Panel3{
private:
    Vector3 plane_n, plane_p;
    PanelVector3* mirror_point(PanelVector3* point);
    Panel3* parent;
public:
    SymmetricPanel3(Panel3* parent, Vector3 plane_n, Vector3 plane_p);
    virtual string type() {return "SymmetricPanel3";}
    virtual double get_sigma(){return this->parent->get_sigma();}
    virtual double get_cp(){return this->parent->get_cp();}
    virtual Vector3 get_velocity();
    virtual double get_mue(){return this->parent->get_mue();}
    virtual double get_potential(){return this->parent->get_potential();}
    virtual vector< PanelVector3* > get_points();
    int get_nr();
    void calc_geo();
    virtual void set_neighbours();
};


class SymmetricWakePanel3: public SymmetricPanel3, public WakePanel3{
private:
    WakePanel3* parent;
public:
    string type() {return "SymmetricWakePanel3";}
    SymmetricWakePanel3(WakePanel3* parent, Vector3 plane_n, Vector3 plane_p);;
    void calc_geo();
    int get_nr();
    virtual vector<PanelVector3*> get_points();
    Panel3* get_lower_operating_panel();
    Panel3* get_upper_operating_panel();
};


class Edge{
private:
    void setup_Edge();
public:
    Edge(){};
    Edge(PanelVector3*, PanelVector3*);
    PanelVector3* v1;
    PanelVector3* v2;
    Panel3* p1;                                         // make private + getter 
    Panel3* p2;                                         // make private + getter
    double vorticity;
    double lambda_cut(Vector3 &p, Vector3 &n);		//v1 + lambda_cut * (v2 - v1) = cut_pos
    Vector3 lambda_pos(double lambda);

    Vector3 center();
    Vector3 diff_vector();
    Vector3 t();
    double l();
    bool is_symmetric = false;
    bool is_wake_edge();
};

ostream &operator<<(ostream &output, Panel3&);

#endif 