#include "panel3.h"
#include <Eigen/Cholesky>
#include <algorithm>

using namespace std;

typedef Eigen::Matrix<double, 5, 1> Vector5d;
typedef Eigen::Matrix<double, 5, 5> Matrix5d;

Panel3::~Panel3()
{
    for (PanelVector3* vec: this->points)
    {
        vec->owner -= 1;
        if (vec->owner == 0)
            delete vec;
    }
}

Panel3::Panel3(){};

Panel3::Panel3(std::vector< PanelVector3* > points)
{
    for (PanelVector3*& vec: points)
    {
        this->append_point(vec);
    }
    this->calc_geo();
}



void Panel3::append_point(PanelVector3* v)
{
  this->points.push_back(v);
  v->owner += 1;
  v->panels.push_back(this);
}

vector< PanelVector3*> Panel3::get_points()
{
    return this->points;
}

bool Panel3::is_symmetric()
{
    return this->_is_symmetric;
}

void Panel3::set_symmetric()
{
    this->_is_symmetric = true;
}


void Panel3::calc_geo()
{
    int l;
    this->area = 0.;
    this->center = Vector3(0.,0.,0.);
    this->n = Vector3(0.,0.,0.);
    
    l = this->points.size();
    
    for (PanelVector3* point: this->points){
        this->center += *point;
    }
    this->center /= double(l);
    
    for (int i=0; i < l; i++){
        int j = (i == l - 1 ? 0 : i+1);
        this->side_sum += (*this->points[j] - *this->points[i]).norm();
        Vector3 n_i = (this->center - *this->points[i]).cross(this->center - *this->points[j]);
        this->n -= n_i;
        this->area += n_i.norm();
    }
    this->area /= 2;
    this->n.normalize();
    this->l = (this->center - (*this->points[1] + *this->points[2]) / 2);
    this->l.normalize();
    this->m = this->n.cross(this->l);
}

bool Panel3::is_inside(Vector3& v)
{
    //     check if the point v lies inside the plane
    double coresize = 0.000000000000000000000001;
    if (v == this->center)
    {
        return 1;
    }
    int l = this->points.size();
    for (int i = 0; i < l; i++)
    {
        int j = (i == l-1 ? 0 : i+1);
        if ((*this->points[j]-v).cross(v-*this->points[i]).dot(this->n) <= coresize)
        {
            return 0;
        }
    }
    return 1;
}

void Panel3::compute_gradient(Vector3& v_inf)
{
    this->ls_gradient();
}


void Panel3::ls_gradient()
{
    //   computes the local velocity from linear leastsquare potential gradient
    
    int l = this->neighbours.size();
    bool quadratic = l > 1;
    Eigen::MatrixXd mat(l + 1, 3 + 2 * quadratic);
    Eigen::VectorXd rhs(l + 1);
    Eigen::VectorXd sol(l + 1);
    rhs.setZero();
    mat.setZero();
    for (int i = 0; i <= l; i++)
    {
        PanelVector3 proj (0,0,0);
        Vector3 loc_pos;
        if (i != l)
        {
            proj = project_to_panel(*this->neighbours[i]);
            loc_pos = this->local_coordinates(proj);
        }
        else
        {
            proj.potential = this->mue;
            loc_pos = Vector3(0, 0, 0);
        }
        rhs(i) = proj.potential;
        mat(i, 0) = 1;
        mat(i, 1) = loc_pos.x();
        mat(i, 2) = loc_pos.y();
        if (quadratic)
        {
            mat(i, 3) = loc_pos.x() * loc_pos.x();
            mat(i, 4) = loc_pos.y() * loc_pos.y();
        }
    }
    sol = mat.colPivHouseholderQr().solve(rhs);
    this->velocity = this->m * sol[1] + this->l * sol[2];
}

Vector3 Panel3::local_coordinates(Vector3 a)
{
    Vector3 point = a - this->center;
    Vector3 out (point.dot(this->m), point.dot(this->l), point.dot(this->n));
    return out;
}


void Panel3::set_neighbours()
{
        if (this->neighbours_set){return;}
    //   get all connected neighbours
    int l = this->points.size();
    for (int i = 0; i < l ; i++){
        int j = (i == l-1 ? 0 : i+1);
        if (! (this->points[i]->wake_vertex && this->points[j]->wake_vertex)){
            for (Panel3*& panel_k: this->points[i]->panels){
                for (Panel3*& panel_m: this->points[j]->panels){
                    if (panel_k == panel_m){
                        if (panel_k != this){
                            if (panel_k->n.dot(this->n) > 0.5){  // not add sharp edge panels
                                this->neighbours.push_back(panel_m);
                            }
                        }
                    }
                }
            }
        }
    }
    this->neighbours_set = true;
}


PanelVector3 Panel3::project_to_panel(Panel3& pan)
{
    //   get the vertices connecting both panels
    Vector3* edge [2];
    int input = 0;
    for (PanelVector3*& point_i: this->points)
    {
        for (PanelVector3*& point_j: pan.points)
        {
            if (point_i == point_j)
            {
                edge[input] = point_i;
                input++;
            }
            if (input == 2){break;}
        }
    if (input == 2){break;}
    }
    // computing shortest distance between center and neighbour center on the surface
    Vector3 &v0 = *edge[0];
    Vector3 &v1 = *edge[1];
    Vector3 v = v1 - v0;
    Vector3 c2 = pan.center;
    Vector3 c1 = this->center;
    // find shortest distance between center points and side 
    Vector3 proj_1 = v0 + v * (v.dot(c1 - v0) / v.dot(v));
    Vector3 proj_2 = v0 + v * (v.dot(c2 - v0) / v.dot(v));
    Vector3 midpoint = (proj_1 + proj_2) / 2;
    Vector3 dir = (midpoint - c1);
    double length = dir.norm() + (midpoint - c2).norm();
    dir.normalize();
    PanelVector3 out = c1 + (dir * length);
    out.potential = pan.get_mue();				// mue is the jump of the potential on the surface
    return out;
}

std::string Panel3::tostring(){
    std::ostringstream out;
    out << "Panel3(";
    for (int i = 0; i < this->points.size(); i++)
    {
        out << this->points[i];
        out << (i == this->points.size() - 1 ? ")" : ", ");
    }
    return out.str();
}


std::ostream& operator<<(std::ostream& output, Panel3& p){
    output <<p.tostring();
    return output;
}


void Panel3::compute_pressure(Vector3& v_inf)
{
    this->cp = 1 - pow((this->velocity.norm() / v_inf.norm()),2);
}

Vector3 Panel3::get_force()
{
    return this->n * (this->get_cp() * this->area);
}

vector< Edge > Panel3::get_edges()
{
    int i, j;
    vector<Edge> panel_edges;
    for (i = 0; i < this->points.size(); i++)
    {
        j = (i == this->points.size() - 1? 0: i + 1);
        Edge e(this->points[i], this->points[j]);
        panel_edges.push_back(e);
    }
    return panel_edges;
}


Panel3* Panel3::get_edge_neighbour(Edge e)
{
    if (e.p1 == this)
    {
        return e.p2;
    }
    return e.p1;
}


int Panel3::aligned(PanelVector3* a, PanelVector3* b)
{
  int ind_a, ind_b;
  bool found_a = false;
  bool found_b = false;
//  return 0 if edge is not part of panel
//  return 1 if aligned to panel
//  return -1 if not aligned to panel
  
  for (ind_a=0; ind_a < this->points.size(); ind_a++){
    // test if a in points and break if found
    if (this->points[ind_a] == a){
      found_a = true;
      break;
    }
  }
  if (! found_a){
    cout << "not_found_a" << endl;
    return 0;
  };
  for (ind_b=0; ind_b < this->points.size(); ind_b++){
    // test if b in points and break if found
    if (this->points[ind_b] == b){
      found_b = true;
      break;
    }
  }
  if (! found_b){
    cout << "not_found_b" << endl;
    return 0;
  }
  if (ind_a - ind_b == 1){
    return -1; // not aligned
  }
  else if (ind_a - ind_b == -1){
    return 1; // aligned
  }
  else if (ind_a == 0){
    return -(ind_b == this->points.size() - 1); // not aligned
  }
  else if (ind_b == 0){
    return (ind_a == this->points.size() - 1); // aligned
  }
  cout << "diagonal" << endl;
  return 0;  // diagonal
}

// ######################WAKE PANEL 3 ###################################

bool WakePanel3::is_symmetric()
{
    return (this->upper_operating_panel->is_symmetric() &&
            this->lower_operating_panel->is_symmetric());
}

Panel3* WakePanel3::get_lower_operating_panel()
{
    return this->lower_operating_panel;
}

Panel3* WakePanel3::get_upper_operating_panel()
{
    return this->upper_operating_panel;
}

void WakePanel3::set_lower_operating_panel(Panel3* pan)
{
    this->lower_operating_panel = pan;
}

void WakePanel3::set_upper_operating_panel(Panel3* pan)
{
    this->upper_operating_panel = pan;
}

// #######################SYMMETRIC PANEL 3##############################

SymmetricWakePanel3::SymmetricWakePanel3(WakePanel3* parent, Vector3 plane_n, Vector3 plane_p): 
    SymmetricPanel3(parent, plane_n, plane_p)
{
    this->parent = parent;
}

void SymmetricWakePanel3::calc_geo()
{
    SymmetricPanel3::calc_geo();
}

Panel3* SymmetricWakePanel3::get_lower_operating_panel()
{
    return this->parent->get_lower_operating_panel();
}

Panel3* SymmetricWakePanel3::get_upper_operating_panel()
{
    return this->parent->get_upper_operating_panel();
}

int SymmetricWakePanel3::get_nr()
{
    return SymmetricPanel3::get_nr();
}

Vector3 SymmetricPanel3::get_velocity()
{
    return this->parent->get_velocity() - this->parent->get_velocity().dot(plane_n) * plane_n * 2;
}


int SymmetricPanel3::get_nr()
{
    return this->parent->get_nr();
}

std::vector< PanelVector3* > SymmetricWakePanel3::get_points()
{
    return SymmetricPanel3::get_points();
}

vector<PanelVector3*> SymmetricPanel3::get_points()
{   // call in constructor
    vector<PanelVector3*> points;
    for(PanelVector3* parent_point: this->parent->points){
        PanelVector3* mirrored_point = this->mirror_point(parent_point);
        mirrored_point->panels.push_back(this);
        mirrored_point->owner ++;
        points.push_back(mirrored_point);
    }
    reverse(points.begin(), points.end());
    this->points = points;
    return Panel3::get_points();
}

SymmetricPanel3::SymmetricPanel3(Panel3* parent, Vector3 plane_n, Vector3 plane_p){
    // creating relation
    this->parent = parent;
    // set the symmetry plane
    this->plane_n = plane_n;
    this->plane_p = plane_p;
}


PanelVector3* SymmetricPanel3::mirror_point(PanelVector3* point)
{
//      if point lies on the symmetry plane return the input pointer
//      else: return a new mirrored panelvector
    PanelVector3* mirrored_point = new PanelVector3(*point);
    *mirrored_point = *point - (plane_n.dot(*point - plane_p)) * plane_n * 2;
    if ((*mirrored_point - *point).norm() < 0.00000000000001)
    {
        delete mirrored_point;
        return point;
    }
    return mirrored_point;
}

void SymmetricPanel3::calc_geo()
{
    this->get_points();
    Panel3::calc_geo();
}


void SymmetricPanel3::set_neighbours()
{
    return;
}


Edge::Edge(PanelVector3* v1, PanelVector3* v2){
  this->v1 = v1;
  this->v2 = v2;
  this->setup_Edge();
}

void Edge::setup_Edge()
{
    // Überprüfen, ob v1 und v2 null sind
    if (!this->v1 || !this->v2) {
        std::cout << "Edge::setup_Edge(): v1 or v2 is null" << std::endl;
        return;
    }
    if (this->v1->panels.size() == 0 || this->v2->panels.size() == 0){
        return;
    }

//   assuming that a edge has two connected panels
  vector<Panel3*> connected_panels;
  for (Panel3*& panel_i: this->v1->panels){
    for (Panel3*& panel_j: this->v2->panels){
      if (panel_i == panel_j){
	connected_panels.push_back(panel_i);
      }
    }
  }
// the trailing edge has 3 connected panels (p1, p2, wp)
// how to order these panels?
  if (connected_panels.size() != 2){
    return;
  }

//    ordering the panels (first Panel3 has [..., this->v1, this->v2, ...] order) (is aligned)
  bool first = false;
  this->p1 = connected_panels[0];
  this->p2 = connected_panels[1];
  int l = this->p1->points.size();
  for (int i = 0; i < l; i++){
    int j = (i == l ? 0 : i + 1);
    if (this->p1->points[i] == this->v1 && this->p1->points[j] == this->v2){
      return;
    }      
  }
  this->p1 = connected_panels[1];
  this->p2 = connected_panels[0];
  return;
}

bool Edge::is_wake_edge()
{
    return (v1->wake_vertex && v2->wake_vertex);
}


double Edge::lambda_cut(Vector3 &p, Vector3 &n){
  return (p-*this->v1).dot(n) / (*this->v2 - *this->v1).dot(n);
}

Vector3 Edge::lambda_pos(double lambda){
  return *this->v1 + lambda * (*this->v2 - *this->v1);
}

Vector3 Edge::center(){
  return (*this->v1 + *this->v2) / 2;
}

Vector3 Edge::t(){
  return this->diff_vector() / this->l();
}

Vector3 Edge::diff_vector(){
  return (*this->v2 - *this->v1);
}

double Edge::l(){
   return this->diff_vector().norm();
}
