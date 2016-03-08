#include "element_influence.h"
#include "panel2.h"
#include <stdlib.h> 
#include <iostream>
#include <math.h>

using namespace std;
double coresize = 0.0000000000000000000001;


void doublet_3_0_n0(Vector3& target, Panel3* source, double& dip_infl)
{
  Vector3 p_v = target - source->center;
  double pn = p_v.dot(source->n);
  dip_infl = -pn * source->area / pow(p_v.norm(), 3) / 4 / M_PI;
}

void doublet_3_0_vsaero_v(Vector3& target, Panel3* source, Vector3& vel_infl)
{
  Vector3 p_v = target - source->center;
  int l = source->points.size();

  for (int i = 0; i < l; i++){
    int j = (i == l-1? 0 : i+1);
    Vector3 &p1 = *source->points[i];
    Vector3 &p2 = *source->points[j];
    Vector3 a = target - p1;
    Vector3 b = target - p2;
    double absa = a.norm();
    double absb = b.norm();
    if ((absa * absb * fabs(absa * absb + a.dot(b))) > coresize )
    {
      vel_infl -= (a.cross(b) * (absa + absb) / 4 / M_PI / (absa * absb * (absa * absb + a.dot(b))));
    }
  }
}

Vector3 doublet_3_0_vsaero_v(Vector3& target, Panel3* source)
{
  Vector3 vel = Vector3(0, 0, 0);
  Vector3 p_v = target - source->center;
  int l = source->points.size();

  for (int i = 0; i < l; i++){
    int j = (i == l-1? 0 : i+1);
    Vector3 &p1 = *source->points[i];
    Vector3 &p2 = *source->points[j];
    Vector3 a = target - p1;
    Vector3 b = target - p2;
    double absa = a.norm();
    double absb = b.norm();
    if ((absa * absb * fabs(absa * absb + a.dot(b))) > 0.000001 )
    {
      vel -= (a.cross(b) * (absa + absb) / 4 / M_PI / (absa * absb * (absa * absb + a.dot(b))));
    }
  }
  return vel;
}


void doublet_3_0_vsaero(Vector3 & target, Panel3* source, double& dip_infl)
{
  Vector3 p_v = target - source->center;
  double pn = p_v.dot(source->n);
  int l = source->points.size();
  double rnum;
  double dnom;


  for (int i = 0; i < l; i++){
    int j = (i == l-1? 0 : i+1);
    Vector3 &p1 = *source->points[i];
    Vector3 &p2 = *source->points[j];
    Vector3 s = p2 - p1;
    Vector3 a = target - p1;
    Vector3 b = target - p2;
    Vector3 h = a.cross(s);
    double bm = b.dot(source->m);
    double am = a.dot(source->m);
    double sl = s.dot(source->l);
    double sm = s.dot(source->m);
    double al = am * sl - a.dot(source->l) * sm;
    double pa = pow(pn, 2) * sl + al * am;
    double pb = pow(pn, 2) * sl + al * bm;
    double absa = a.norm();
    double absb = b.norm();
    double abss = s.norm();
    
    if (
      absa > coresize && 
      absb > coresize &&
      absa + absb - abss >= 0
    )
    {
    dnom = pa * pb + pow(pn, 2) * absa * absb * pow(sm, 2);
    rnum = sm * pn * (absb * pa - absa * pb);

    if (fabs(pn) < coresize)
    {
      double pn_sign = (pn >= 0 ? -1 : 1);
      double al_sign = (al < 0 ? 1 : -1);
      if (dnom < 0.){dip_infl += al_sign * pn_sign * M_PI;}
      else if (dnom == 0){dip_infl += al_sign * pn_sign * M_PI / 2;}
    }
    else
        {
        dip_infl += atan2(rnum, dnom);
        }
    }
  }
  dip_infl /= 4 * M_PI;
}

double doublet_3_0_sphere(Vector3& target, Panel3* source, double& dip_infl)
{
    Vector3 u1, u2, u3;
    
    u1 = source->center - target;
    if (u1.norm() < 0.000000000000000000001){
        // point lies on the center of the panel
        // this will return the value of the integral, with the panel wrapped around the center point
        // in n direction.
        dip_infl  = -0.5;
        return dip_infl;
    }
    u1.normalize();
    for (int i = 0; i < source->points.size(); i++){
        int j = (i == source->points.size()-1? 0 : i+1);
        u2 = *source->points[i] - target;
        u3 = *source->points[j] - target;
        u2.normalize();
        u3.normalize();
        dip_infl -= atan2((u1.cross(u2)).dot(u3), (1 + u2.dot(u3) + u3.dot(u1) + u1.dot(u2))) * 0.5 / M_PI;
    }
    return dip_infl;
}


void doublet_src_3_0_n0(Vector3& target, Panel3* source, double& dip_infl, double& mop_infl)
{
  Vector3 p_v = target - source->center;
  double absp = p_v.norm();
  double pn = p_v.dot(source->n);
  dip_infl = -pn * source->area / pow(absp, 3) / 4 / M_PI;
  mop_infl = source->area / absp / 4 / M_PI;
}

Vector3 doublet_src_3_0_vsaero_v(Vector3& target, Panel3* source, double mue, double sigma)
{
    Vector3 vel = Vector3(0, 0, 0);
    double dip_infl = 0;
    Vector3 p_v = target - source->center;
    double pn = p_v.dot(source->n);
    int l = source->points.size();
    double rnum;
    double dnom;


    for (int i = 0; i < l; i++){
        int j = (i == l-1? 0 : i+1);
        Vector3 &p1 = *source->points[i];
        Vector3 &p2 = *source->points[j];
        Vector3 s = p2 - p1;
        Vector3 a = target - p1;
        Vector3 b = target - p2;
        Vector3 h = a.cross(s);
        double bm = b.dot(source->m);
        double am = a.dot(source->m);
        double sl = s.dot(source->l);
        double sm = s.dot(source->m);
        double al = am * sl - a.dot(source->l) * sm;
        double pa = pow(pn, 2) * sl + al * am;
        double pb = pow(pn, 2) * sl + al * bm;
        double absa = a.norm();
        double absb = b.norm();
        double abss = s.norm();
        
        if (
        absa > coresize && 
        absb > coresize &&
        absa + absb - abss >= 0
        ){
        dnom = pa * pb + pow(pn, 2) * absa * absb * pow(sm, 2);
        rnum = sm * pn * (absb * pa - absa * pb);

        if (fabs(pn) < coresize){
            double pn_sign = (pn >= 0 ? -1 : 1);
            double al_sign = (al < 0 ? 1 : -1);
            if (dnom < 0.){dip_infl += al_sign * pn_sign * M_PI;}
            else if (dnom == 0){dip_infl += al_sign * pn_sign * M_PI / 2;}
        }
        else{
            dip_infl += atan2(rnum, dnom);
        }    if ((absa * absb * fabs(absa * absb + a.dot(b))) > coresize )
        vel -= mue * (a.cross(b) * (absa + absb) / (absa * absb * (absa * absb + a.dot(b))));
        vel -= sigma * log((absa + absb + abss) / (absa + absb - abss)) / abss * (sm * source->l - sl * source->m);
        }
    }
    vel -= sigma * source->n * dip_infl;
    return vel / 4 / M_PI;
}

Vector3 src_3_0_vsaero_v(Vector3& target, Panel3* source, double sigma)
{
    Vector3 vel = Vector3(0, 0, 0);
    double dip_infl = 0;
    Vector3 p_v = target - source->center;
    double pn = p_v.dot(source->n);
    int l = source->points.size();
    double rnum;
    double dnom;


    for (int i = 0; i < l; i++){
        int j = (i == l-1? 0 : i+1);
        Vector3 &p1 = *source->points[i];
        Vector3 &p2 = *source->points[j];
        Vector3 s = p2 - p1;
        Vector3 a = target - p1;
        Vector3 b = target - p2;
        Vector3 h = a.cross(s);
        double bm = b.dot(source->m);
        double am = a.dot(source->m);
        double sl = s.dot(source->l);
        double sm = s.dot(source->m);
        double al = am * sl - a.dot(source->l) * sm;
        double pa = pow(pn, 2) * sl + al * am;
        double pb = pow(pn, 2) * sl + al * bm;
        double absa = a.norm();
        double absb = b.norm();
        double abss = s.norm();
        
        if (
        absa > coresize && 
        absb > coresize &&
        absa + absb - abss >= 0
        ){
        dnom = pa * pb + pow(pn, 2) * absa * absb * pow(sm, 2);
        rnum = sm * pn * (absb * pa - absa * pb);

        if (fabs(pn) < coresize){
            double pn_sign = (pn >= 0 ? 1 : -1);
            double al_sign = (al < 0 ? 1 : -1);
            if (dnom < 0.){dip_infl += al_sign * pn_sign * M_PI;}
            else if (dnom == 0){dip_infl += al_sign * pn_sign * M_PI / 2;}
        }
        else{
            dip_infl += atan2(rnum, dnom);
        }    if ((absa * absb * fabs(absa * absb + a.dot(b))) > coresize )
        vel -= sigma * log((absa + absb + abss) / (absa + absb - abss)) / abss * (sm * source->l - sl * source->m);
        }
    }
    vel -= sigma * source->n * dip_infl;
    return vel / 4 / M_PI;
}


void doublet_src_3_0_vsaero(Vector3& target, Panel3* source, double& dip_infl, double& mop_infl)
{
    Vector3 p_v = target - source->center;
    double pn = p_v.dot(source->n);
    int l = source->points.size();
    double rnum;
    double dnom;


    for (int i = 0; i < l; i++){
        int j = (i == l-1? 0 : i+1);
        Vector3 &p1 = *source->points[i];
        Vector3 &p2 = *source->points[j];
        Vector3 s = p2 - p1;
        Vector3 a = target - p1;
        Vector3 b = target - p2;
        Vector3 h = a.cross(s);
        double bm = b.dot(source->m);
        double am = a.dot(source->m);
        double sl = s.dot(source->l);
        double sm = s.dot(source->m);
        double al = am * sl - a.dot(source->l) * sm;
        double pa = pow(pn, 2) * sl + al * am;
        double pb = pow(pn, 2) * sl + al * bm;
        double absa = a.norm();
        double absb = b.norm();
        double abss = s.norm();
        
        if (
        absa > coresize && 
        absb > coresize &&
        absa + absb - abss >= 0
        )
        {
        dnom = pa * pb + pow(pn, 2) * absa * absb * pow(sm, 2);
        rnum = sm * pn * (absb * pa - absa * pb);

        if (fabs(pn) < coresize)
        {
            double pn_sign = (pn >= 0 ? -1 : 1);
            double al_sign = (al < 0 ? 1 : -1);
            if (dnom < 0.){dip_infl += al_sign * pn_sign * M_PI;}
            else if (dnom == 0){dip_infl += al_sign * pn_sign * M_PI / 2;}
        }
        else
            {
            dip_infl += atan2(rnum, dnom); 
            }
        mop_infl -= al * log((absa + absb + abss) / (absa + absb - abss)) / abss;  
        }
    }
    mop_infl += pn * dip_infl;
    dip_infl /= 4. * M_PI;
    mop_infl /= 4. * M_PI;
}


Vector3 vortex_3_0_v(Vector3& target, Vector3& source_point_0, Vector3& source_point_1)
{
  Vector3 r0 = source_point_1 - source_point_0;
  Vector3 r1 = target - source_point_0;
  Vector3 r2 = target - source_point_1;
  Vector3 r_cross = r1.cross(r2);
  if (r_cross.norm() > 0.001){
    return r_cross * ((r0.dot(r1) / r1.norm() - r0.dot(r2) / r2.norm()) / r_cross.dot(r_cross) / 4 / M_PI);
  }
  return Vector3(0.,0.,0.);
}

Vector3 vortex_3_0_v(Vector3& target, Edge& e)
{
    return e.vorticity * vortex_3_0_v(target, *e.v1, *e.v2);
}



Vector3 vortex_3_0_half_infinity_v(Vector3& target, Vector3& source_direction, Vector3& source_point)
{
  Vector3 r1 = target - source_point;
  Vector3 cross = r1.cross(source_direction);
  if (cross.norm() > 0.00001){
    return cross / cross.dot(cross) * (1 + r1.dot(source_direction) / r1.norm()) / 4 / M_PI;
  }
  return Vector3(0.,0.,0.);
}

double doublet_3(Vector3& target, Vector3& source, Vector3& normal)
{
    Vector3 r = target - source;
    return -1. / 4 / M_PI / pow(r.dot(r), (3. / 2)) * r.dot(normal);
}

Vector3 doublet_3_v(Vector3& target, Vector3& source, Vector3& normal)
{
    Vector3 r = target - source;
    return ((normal.cross(r)).cross(r) + 2. * r * (r.dot(normal))) / 4. / M_PI / pow(r.dot(r), 5. / 2.);
}



/************************-2D ELEMENT INFLUENCE-*********************/

double doublet_2_0(Vector2& target, Panel2& source)
{
    double pn = (target - source.center).dot(source.n);
    double dx1 = (*source.points[0] - target).dot(source.t);
    double dx2 = (*source.points[1] - target).dot(source.t);
    if (pn * pn < coresize and (dx1 * dx1 + dx2 * dx2) <= (source.l * source.l))
    {
        return -0.5;
    }
    return - (atan2(pn , dx2) - atan2(pn , dx1)) / (2. * M_PI);  
}


double source_2_0(Vector2& target, Panel2& source)
{ 
    double pn = (target - source.center).dot(source.n);
    double dx1 = (*source.points[0] - target).dot(source.t);
    double dx2 = (*source.points[1] - target).dot(source.t);
    double a = (dx1 * log(dx1 * dx1 + pn * pn) - dx2 * log(dx2 * dx2 + pn * pn));
    double b = 2 * pn * (atan2(pn, dx2) -atan2(pn, dx1));
    return -(a + b) / (4. * M_PI);
}

Vector2 doublet_2_0_v(Vector2& target, Panel2& source)
{
    double pn = (target - source.center).dot(source.n);
    Vector2 r1 = target - *source.points[0];
    Vector2 r2 = target - *source.points[1];
    double dx1 = r1.dot(source.t);
    double dx2 = r2.dot(source.t);
    
    Vector2 u = (pn / r1.dot(r1) - pn / r2.dot(r2)) * source.t;
    Vector2 w = (dx1 / r1.dot(r1)- dx2 / r2.dot(r2)) * source.n;
    return - (w - u) / M_PI / 2;
}

Vector2 source_2_0_v(Vector2& target, Panel2& source)
{ 
    double pn = (target - source.center).dot(source.n);
    if (pn == 0)
    {
	return -0.5 * source.n;
    }
    double dx1 = (target - *source.points[0]).dot(source.t);
    double dx2 = (target - *source.points[1]).dot(source.t);
    
    Vector2 u = log((dx1 * dx1 + pn * pn) / (dx2 * dx2 + pn * pn) ) / 2  * source.t;
    Vector2 w = (atan2(pn, dx2) - atan2(pn, dx1)) * source.n;
    return (u + w) / M_PI / 2;
}

double vortex_2(Vector2& target, Vector2& source, Vector2& direction)
{
    direction.normalize();
    double pn = (target - source).dot(normal2(direction));
    double dx = (source - target).dot(direction);
    return - atan2(pn , dx) / (2. * M_PI);  
}

Vector2 vortex_2_v(Vector2& target, Vector2& source)
{
    Vector2 r = target - source;
    if (r.norm() < coresize)
    {
        r.normalize();
        r *= coresize;
    }
    return normal2(r / r.dot(r) / 2 / M_PI);
}


double doublet_2_1(Vector2& target, Panel2& source, bool left)
{
    /* return a linear distribution doublet distribution from 1..0 if left is true
     * otherwise the distribution is from 0...1
     * for a distribution from mue2 to mue1 this function has to be called twice 
     *          mue1 * doublet_2_1(,,True) + mue2 * doublet_2_1(,,false)
     * or call a constant doublet + this linear doublet
     *          mue1 * doublet_2_0() + (mue2 - mue1) * doublet_2_1(,,false)*/
    
    double pn = (target - source.center).dot(source.n);
    Vector2 p1_t = target - *source.points[0];
    Vector2 p2_t = target - *source.points[1];
    double dx1 = p1_t.dot(source.t);
    double dx2 = p2_t.dot(source.t);
    double r1_squared = p1_t.dot(p1_t);
    double r2_squared = p2_t.dot(p2_t);
    double phi1 = atan2(pn, dx1);
    double phi2 = atan2(pn, dx2);
    if (pn * pn < coresize)
    {
        if (r1_squared < (source.l * source.l) and
        r2_squared < (source.l * source.l))
        {
            if (left){return -0.5 * p2_t.norm() / source.l;}
            else {return -0.5 * p1_t.norm() / source.l;}
        }
        return 0.;
    }
    if (left)
    {
        return  -0.5 / M_PI / source.l * (dx1 * (phi2 - phi1) + pn / 2. * log(r2_squared / r1_squared));
    }
    return 0.5 / M_PI / source.l * (dx2 * (phi2 - phi1) + pn / 2. * log(r2_squared / r1_squared)); 
}

Vector2 doublet_2_1_v(Vector2& target, Panel2& source, bool left)
{
    double pn = (target - source.center).dot(source.n);
    Vector2 p1_t = target - *source.points[0];
    Vector2 p2_t = target - *source.points[1];
    double dx1 = p1_t.dot(source.t);
    double dx2 = p2_t.dot(source.t);
    double r1_squared = p1_t.dot(p1_t);
    double r2_squared = p2_t.dot(p2_t);
    double phi1 = atan2(pn, dx1);
    double phi2 = atan2(pn, dx2);
    
    if (left)
    {
        Vector2 u = ((phi2 - phi1) * source.l + pn / r2_squared) * source.t;
        Vector2 w = (0.5 * source.l * log(r2_squared / r1_squared)  - dx2 / r2_squared) * source.n;
        return -0.5 / M_PI * (u + w);
    }
    else
    {
        Vector2 u = ((phi2 - phi1) * source.l + pn / r1_squared) * source.t;
        Vector2 w = (0.5 * source.l * log(r2_squared / r1_squared)  - dx1 / r1_squared) * source.n;
        return 0.5 / M_PI * (u + w);
    }
    
}


double source_2_1(Vector2& target, Panel2& source, bool left)
{
    double pn = (target - source.center).dot(source.n);
    Vector2 p1_t = target - *source.points[0];
    Vector2 p2_t = target - *source.points[1];
    double dx1 = p1_t.dot(source.t);
    double dx2 = p2_t.dot(source.t);
    double r1_squared = p1_t.dot(p1_t);
    double r2_squared = p2_t.dot(p2_t);
    double phi1 = atan2(pn, dx1);
    double phi2 = atan2(pn, dx2);    
}


Vector2 source_2_1_v(Vector2& target, Panel2& source, bool left)
{

}


double source_2(Vector2& target, Vector2& source)
{
    return 1. / 2. / M_PI * log((target - source).norm());
}

Vector2 source_2_v(Vector2& target, Vector2& source)
{
    Vector2 r= target - source;
    return 1. / 2. / M_PI * r / r.dot(r);
}

double doublet_2(Vector2& target, Vector2& source, Vector2& direction)
{
    Vector2 r = target - source;
    direction.normalize();
    return - 1. / 2. / M_PI * direction.dot(r) / r.dot(r);
}


Vector2 doublet_2_v(Vector2& target, Vector2& source, Vector2& direction)
{
    direction.normalize();
    Vector2 r = (target - source);
    double r_dot_r = r.dot(r); r.normalize();
    Vector2 v = r.dot(direction) * r + r.dot(normal2(direction)) * normal2(r);
    return 1. / 2. / M_PI / r_dot_r * v;
}



