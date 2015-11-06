/**
    Purpose: Geometric properties of 2d-Panels
    @author L Lechner
*/


#ifndef panel2_H
#define panel2_H
#include <iostream>
#include <sstream>
#include <vector>

#include "element_influence.h"
#include "vector.h"
using namespace std;


/**
 * The Panel2 is a line segment created from two points. It has some geometric properties which are set in
 * the calc_geo methode. 
 */
class Panel2{
private:
    Panel2* get_neigbour(int, Panel2*);
public:
    Panel2(){};
    ~Panel2(){};
    Panel2(PanelVector2*, PanelVector2*);
    vector<PanelVector2*> points;
    virtual vector<PanelVector2*> get_points();
    void append_point(PanelVector2*);
    double l = 0.;				//panel length
    int nr = 0;
    Vector2 center;
    Vector2 n;
    Vector2 t;
    double potential;
    Vector2 velocity;
    double cp = 0;
    double mue = 0;
    double sigma = 0;
    void calc_geo();
    void compute_gradient();
    void compute_pressure(Vector2& v_inf);
    Vector2 force();
};

/**
 * A WakeVertex acting on a upper and lower PanelVector.
 */
class PanelWakeVertex 
{
public:
    PanelWakeVertex(PanelVector2* point);
    ~PanelWakeVertex(){};
    PanelVector2* panel_vector;
    Panel2* upper_operating_panel;
    Panel2* lower_operating_panel;
};

/**
 * A WakeVertex acting on a upper and lower PanelVector. Used for methodes where the unknowns
 * are located at the corners of the panel eg higher order methodes.
 */
class PointWakeVertex
{
public:
    PointWakeVertex(PanelVector2* point1, PanelVector2* point2);
    ~PointWakeVertex(){};
    PanelVector2* upper_point;
    PanelVector2* lower_point;
};


#endif 