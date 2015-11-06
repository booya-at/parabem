#include "panel2.h"

void Panel2::append_point(PanelVector2* v)
{
    this->points.push_back(v);
    v->panels.push_back(this);
}

Panel2::Panel2(PanelVector2* p1, PanelVector2* p2)
{
    this->append_point(p1);
    this->append_point(p2);
    this->calc_geo();
}


vector< PanelVector2* > Panel2::get_points()
{
    return this->points;
}


void Panel2::calc_geo()
{
  this->center = (*this->points[0] + *this->points[1])/2;
  this->t = *this->points[0] - *this->points[1];
  this->l = this->t.norm();
  this->t = this->t / this->l;
  this->n = Vector2(this->t[1], -this->t[0]);
}

Panel2* Panel2::get_neigbour(int nr, Panel2* pan1)
{
  Panel2* pan = pan1->points[nr]->panels[0];
  if (pan == this && pan1 == this){
    return pan1->points[nr]->panels[1];
  }
  if (pan == pan1){
    pan = pan1->points[nr]->panels[1];
    if (pan == this){
      pan = pan1->points[1 - nr]->panels[0];
      if (pan = pan1){
	pan = pan1->points[1 - nr]->panels[1];
      }
    }
  }
  return pan;
}


void Panel2::compute_gradient()
{
//   compute the gradient of a quadric interpolation from the neighbours potential
//   C0 + C1*x + C2*x^2 = p(x)
//   C0 = p(0) = p0
//   C0 + C1*(l1) + C2*l1^2 = p(-l1) = p1
//   C0 + C1*l2 + C2*l2^2 = p(l2) = p2
//   C1 = (l2^2*(p0 - p1) + l1^2*(-p0 + p2))/(l1*l2*(l1 - l2))
//   if one corner of the panel is a wake-vertex the gradient the neigbours are choosen from
//   one side only.
  Panel2* pan1;
  Panel2* pan2;
  pan1 = this->get_neigbour(0, this);
  pan2 = this->get_neigbour(1, this);
  double l1 = -(pan1->l + this->l) / 2;
  double l2 = (pan2->l + this->l) / 2;
  if (this->points[0]->wake_vertex){
    pan1 = this->get_neigbour(1, pan2);
    l1 = pan2->l + (pan1->l + this->l) / 2;
  }
  else if (this->points[1]->wake_vertex){
    pan2 = this->get_neigbour(0, pan1);
    l2 = - pan1->l - (pan2->l + this->l) / 2;
  }
    
  double p0 = this->mue;
  double p1 = pan1->mue;
  double p2 = pan2->mue;
  this->velocity = -this->t * (l2 * l2 *(p0 - p1) + l1 * l1 *(p2 - p0)) / (l1 * l2 * (l1 - l2));
}

void Panel2::compute_pressure(Vector2& v_inf)
{
  this->cp = 1 - pow((this->velocity.norm() / v_inf.norm()),2);
}

Vector2 Panel2::force()
{
  return this->n * (this->cp * this->l);
}



PanelWakeVertex::PanelWakeVertex(PanelVector2* point)
{
    this->panel_vector = point;
    if (point == point->panels[0]->points[0])
    {
        this->upper_operating_panel = point->panels[0];
        this->lower_operating_panel = point->panels[1];
    }
    else
    {
        this->upper_operating_panel = point->panels[1];
        this->lower_operating_panel = point->panels[0];
    }
}


PointWakeVertex::PointWakeVertex(PanelVector2* point1, PanelVector2* point2)
{
    this->upper_point = point1;
    this->lower_point = point2;
}


