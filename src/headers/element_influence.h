/**
    Purpose: functions to calculate the influence of a singularity-panel
    @author L Lechner
*/


#ifndef element_influence_H
#define element_influence_H

#include "vector.h"
#include "panel3.h"
#include "panel2.h"

/************************-NAMING CONVENTION-************************/

/*
the influence function
    (1) element (doublet, monopol, vortex), 
    (2) dimension(2,3),
    (3) order (0,1,2), if no value -> point element
    (4) integration_sheme,
    [5] v... velocity

the arguments
    first argument is always the target object (point in the space)
    second argument is the perturbation object (panel, line, point)

variable names
    vector in space                    point        or  pnt
    a Panel                            panel        or  pan
    the influence object               source       or  src
    the target object                  target       or  trg
    a monopol or monopol panel         monopol      or  mop
    a doublet or doublet panel         doublet      or  dip
    a velocity vector                  velocity     or  vel
    the potential                      potential    or  pot
    the jump in potential                               mue
    the jump of the pot derivative     sigma        or  sig
    a influence strength               .._influence or  .._infl

*/


/************************-3D ELEMENT INFLUENCE-*********************/
void doublet_3_0_vsaero(Vector3& target, Panel3* source, double& dip_infl);
double doublet_3_0_sphere(Vector3& target, Panel3* source, double& dip_infl);
void doublet_3_0_n0(Vector3& target, Panel3* source, double& dip_infl);
void doublet_3_0_vsaero_v(Vector3& target, Panel3* source, Vector3& vel_infl);
void doublet_src_3_0_vsaero(Vector3 & target, Panel3* source, double& dip_infl, double& src_infl);
void doublet_src_3_0_n0(Vector3 & target, Panel3* source, double& dip_infl, double& src_infl);
double doublet_3(Vector3& target, Vector3& source, Vector3& normal);

Vector3 src_3_0_vsaero_v(Vector3 & target, Panel3* source, double sigma=1);
Vector3 doublet_3_0_vsaero_v(Vector3& target, Panel3* source);
Vector3 doublet_src_3_0_vsaero_v(Vector3 & target, Panel3* source, double mue=1, double sigma=1);


Vector3 vortex_3_0_v(Vector3& target, Vector3& source_point_0, Vector3& source_point_1);
Vector3 vortex_3_0_v(Vector3& target, Edge& e);
Vector3 vortex_3_0_half_infinity_v(Vector3& target, Vector3& source_direction, Vector3& source_point);	
Vector3 doublet_3_v(Vector3& target, Vector3& source, Vector3& normal);


/************************-2D ELEMENT INFLUENCE-*********************/

double source_2(Vector2& target, Vector2& source);
double vortex_2(Vector2& target, Vector2& source, Vector2& direction);
double doublet_2(Vector2& target, Vector2& source, Vector2& direction);
double source_2_0(Vector2& target, Panel2& source);
double doublet_2_0(Vector2& target, Panel2& source);
double doublet_2_1(Vector2& target, Panel2& source, bool left);
double source_2_1(Vector2& target, Panel2& source, bool left);

Vector2 source_2_v(Vector2& target, Vector2& source);
Vector2 doublet_2_v(Vector2& target, Vector2& source, Vector2& direction);
Vector2 vortex_2_v(Vector2& target, Vector2& source);
Vector2 source_2_0_v(Vector2& target, Panel2& source);
Vector2 doublet_2_0_v(Vector2& target, Panel2& source);
Vector2 doublet_2_1_v(Vector2& target, Panel2& source, bool left);
Vector2 source_2_1_v(Vector2& target, Panel2& source, bool left);


#endif //element_influence_H