#include "case2.h"
#include <Eigen/LU>
#include <set>

/************************-2D CASE BASE CLASS-*********************/


/** construction of the case with a list of panels.*/
Case2::Case2(vector<Panel2*> panels)
{
    int i = 0;
    for(Panel2*& panel: panels){
        panel->nr = this->panels.size();
        this->panels.push_back(panel);
        for (PanelVector2*& point: panel->points)
        {
            if (point->nr == -1)
            {
                point->nr = i;
                i++;
            }
        }
    }
    for(PanelVector2*& point: this->get_all_points())
    {
        if (point->wake_vertex)
        {
            PanelWakeVertex w(point);
            this->wake_vertices.push_back(w);
        }
    }
}


/** a set of PanelVectors which is ordered by the operator function */
struct sort_set
{
    bool operator() (const PanelVector2* lhs, const PanelVector2* rhs) const
    {
        return rhs->nr > lhs->nr;
    }
};


/** return all panelvectors ordered by their number.*/
vector< PanelVector2* > Case2::get_all_points()
{
    set<PanelVector2*, sort_set> point_set;
    for(Panel2*& panel: this->panels)
    {
        for(PanelVector2*& point: panel->points)
        {
            point_set.insert(point);
        }
    }
    return vector<PanelVector2*>(point_set.begin(), point_set.end());
}

/** sum all panel-forces and set the case-force and case-moment*/
void Case2::sum_forces()
{
    this->moment = 0;
    this->force.setZero();
    for (int i = 0; i < this->mat_size; i++){
        force += this->panels[i]->force();
        Vector2 distance = (this->panels[i]->center - this->mom_ref_point);
        this->moment += this->panels[i]->force().y() * distance.x();
        this->moment -= this->panels[i]->force().x() * distance.y();
    }
    this->ca = force.dot(Vector2(-this->v_inf.y(), this->v_inf.x()));
    this->center_of_pressure = this->mom_ref_point + Vector2(this->moment / this->force.y(), 0.);
    this->force = force;
}


/************************************************/
void Case2::run(){}					/** solving the boundary integral equation*/
double Case2::off_body_potential(Vector2&){}		/** return the potenial in the region of interest*/
Vector2 Case2::off_body_velocity(Vector2&){}		/** return the velocity in the region of interest*/
/************************************************/

/** return a flow path (list of Vectors)*/
vector< Vector2 > Case2::flow_path(Vector2 start, double inc, int count)
{
    PanelVector2 v1(start);
    vector<Vector2> path;
    path.push_back(v1);
    for (int i = 0; i < count; i++)
    {
        this->off_body_velocity(v1);
        v1.velocity.normalize();
        v1 += v1.velocity * inc;
        path.push_back(v1);
    }
    return path;
}


/************************-DIRICHLET DOUPLET CASE-*********************/
DirichletDoublet0Case2::DirichletDoublet0Case2(vector<Panel2*> panels):
    Case2(panels){}


/** solving the boundary integral equation and set the panel-flow-results (cp, vel, pot) */
void DirichletDoublet0Case2::run()
{

    //   cout << "CALCULATE GEOMETRY" << endl;
    this->mat_size = this->panels.size();
    for (int i = 0; i < this->mat_size; i++)
    {
        this->panels[i]->calc_geo();
    }
    //   cout << "NUMBER OF PANELS = " << this->mat_size << endl;
    this->matrix.resize(this->mat_size, this->mat_size);
    this->matrix.setZero();
    this->rhs.resize(this->mat_size);
    this->rhs.setZero();
    this->result.resize(this->mat_size);
    this->result.setZero();


    //   cout << "FILLING MATRIX & RHS" << endl;

    for (Panel2*& panel_i: this->panels)
    {
        this->rhs(panel_i->nr) -= this->freeflow_influence(*panel_i);
        for (Panel2*& panel_j: this->panels)
        {
            double influence_j_i = this->doublet_influence(*panel_i, *panel_j);
            this->matrix(panel_i->nr, panel_j->nr) += influence_j_i;
            this->rhs(panel_i->nr) -= influence_j_i * this->freeflow_influence(*panel_j);
        };
        for (PanelWakeVertex wake_point: this->wake_vertices)
        {
            double influence = this->wake_influence(*panel_i, *wake_point.panel_vector);
            int nu = wake_point.upper_operating_panel->nr;
            int nl = wake_point.lower_operating_panel->nr;
            matrix(panel_i->nr, nu) -= influence;
            matrix(panel_i->nr, nl) += influence;
        }
    }

//   cout << "SOLVE SYSTEM -> POTENTIAL" << endl;

    this->result = this->matrix.lu().solve(this->rhs);
    for (int i = 0; i < this->mat_size; i++){
        this->panels[i]->mue = this->result(i);  // here mue = mue - phi_inf
    }

//   cout << "COMPUTE VELOCITY" << endl;
    for (Panel2*& panel: this->panels)
    {
        panel->compute_gradient();
        panel->velocity += this->v_inf - panel->n * this->v_inf.dot(panel->n);
        panel->potential = -panel->mue - this->freeflow_influence(*panel);
        panel->compute_pressure(this->v_inf);
    }
    this->sum_forces();
}


//* return the influence of the wake */
double DirichletDoublet0Case2::wake_influence(Panel2& target, Vector2& source)
{
    return vortex_2(target.center, source, this->v_inf);
}


double DirichletDoublet0Case2::off_body_potential(Vector2& vec)
{
    double pot = 0;;
    for (int i = 0; i < this->mat_size; i++){
        pot += -this->panels[i]->potential * doublet_2_0(vec, *this->panels[i]);
    }
    for (PanelWakeVertex wake_point: this->wake_vertices)
    {
        pot += vortex_2(vec, *wake_point.panel_vector) *
                (wake_point.upper_operating_panel->potential -
                 wake_point.lower_operating_panel->potential);
    }
    pot += this->freeflow_influence(vec);
    return pot;
}

Vector2 DirichletDoublet0Case2::off_body_velocity(Vector2& vec)
{
    Vector2 velocity(0, 0);
    Vector2 velocity_influence;
    for (Panel2*& panel: this->panels)
    {
        velocity_influence = doublet_2_0_v(vec, *panel);
        velocity -= panel->potential * velocity_influence;  //negative because the pot is used
    }
    for (PanelWakeVertex wake_point: this->wake_vertices)
    {
        velocity_influence = vortex_2_v(vec, *wake_point.panel_vector);
        velocity += velocity_influence * (
            wake_point.upper_operating_panel->potential -
            wake_point.lower_operating_panel->potential);
    }

    return velocity + this->v_inf;
}


double DirichletDoublet0Case2::doublet_influence(Vector2& target, Panel2& source){
    return doublet_2_0(target, source);
}

double DirichletDoublet0Case2::doublet_influence(Panel2& target, Panel2& source){
    return doublet_2_0(target.center, source);
}

double DirichletDoublet0Case2::freeflow_influence(Vector2& vec){
    return vec.dot(this->v_inf);
}

double DirichletDoublet0Case2::freeflow_influence(Panel2& pan)
{
    return pan.center.dot(this->v_inf);
}




/************************-NEUMANN DOUPLET CASE-*********************/

NeumannDoublet0Case2::NeumannDoublet0Case2(vector<Panel2*> panels):
    Case2(panels){}

double NeumannDoublet0Case2::surface_influence(Panel2& target, Panel2& source)
{
    return doublet_2_0_v(target.center, source).dot(target.n);
}

double NeumannDoublet0Case2::wake_influence(Panel2& target, Vector2& source)
{
    return vortex_2_v(target.center, source).dot(target.n);
}


void NeumannDoublet0Case2::run()
{
//   cout << "CALCULATE GEOMETRY" << endl;
    this->mat_size = this->panels.size();
    for (int i = 0; i < this->mat_size; i++){
        this->panels[i]->calc_geo();
    }
//   cout << "NUMBER OF PANELS = " << this->mat_size << endl;
    this->matrix.resize(this->mat_size, this->mat_size);
    this->matrix.setZero();
    this->rhs.resize(this->mat_size);
    this->rhs.setZero();
    this->result.resize(this->mat_size);
    this->result.setZero();


//   cout << "FILLING MATRIX & RHS" << endl;

    for (Panel2*& panel_i: this->panels)
    {
        this->rhs(panel_i->nr) -= this->v_inf.dot(panel_i->n);
        for (Panel2* panel_j: this->panels)
        {
            double influence_j_i = this->surface_influence(*panel_i, *panel_j);
            this->matrix(panel_i->nr, panel_j->nr) += influence_j_i;
        }
        for (PanelWakeVertex wake_point: this->wake_vertices)
        {
            double influence = this->wake_influence(*panel_i, *wake_point.panel_vector);
            int nu = wake_point.upper_operating_panel->nr;
            int nl = wake_point.lower_operating_panel->nr;
            matrix(panel_i->nr, nu) -= influence;
            matrix(panel_i->nr, nl) += influence;
        }
    }

//   cout << "SOLVE SYSTEM -> POTENTIAL" << endl;
//   using a more accurate matrix solver because
//   the neumann bc has some problems with lu-decomposition
    this->result = this->matrix.fullPivHouseholderQr().solve(this->rhs);
    for (int i = 0; i < this->mat_size; i++){
        this->panels[i]->mue = this->result(i);
    }

//   cout << "COMPUTE VELOCITY" << endl;
    for (int i=0; i < this->mat_size; i++){

    //  velocity = half the duplet gradient over the surface (LSA->11.38)...
        this->panels[i]->velocity = Vector2(0, 0);
        this->panels[i]->compute_gradient();
        this->panels[i]->velocity *= 0.5;

    //  ... + the sum of all panel influences acting on the point of interest (LSA->11.37)
        this->panels[i]->velocity += this->off_body_velocity(this->panels[i]->center);
        this->panels[i]->compute_pressure(this->v_inf);
    }
    this->sum_forces();
}

Vector2 NeumannDoublet0Case2::off_body_velocity(Vector2& vec)
{
    Vector2 velocity(0, 0);
    for (Panel2*& panel: this->panels)
    {
        velocity += panel->mue * doublet_2_0_v(vec, *panel);
    }

    for (PanelWakeVertex wake_point: this->wake_vertices)
    {
        velocity -= vortex_2_v(vec, *wake_point.panel_vector) *
                    (wake_point.upper_operating_panel->mue -
                     wake_point.lower_operating_panel->mue);
    }
    return velocity + this->v_inf;
}

double NeumannDoublet0Case2::off_body_potential(Vector2& vec)
{
    double pot = 0;;
    for (Panel2*& panel: this->panels)
    {
        pot += panel->mue * doublet_2_0(vec, *panel);
    }
    for (PanelWakeVertex wake_point: this->wake_vertices)
    {
        pot -= vortex_2(vec, *wake_point.panel_vector) *
                    (wake_point.upper_operating_panel->mue -
                     wake_point.lower_operating_panel->mue);
    }
    pot += vec.dot(this->v_inf);
    return pot;
}


/************************-DIRICHLET DOUPLET + SOURCE-*********************/

DirichletDoublet0Source0Case2::DirichletDoublet0Source0Case2(vector<Panel2*> panels):
    DirichletDoublet0Case2(panels) {}

double DirichletDoublet0Source0Case2::src_influence(Vector2& target, Panel2& source)
{
    return source_2_0(target, source);
}


double DirichletDoublet0Source0Case2::src_influence(Panel2& target, Panel2& source)
{
    return source_2_0(target.center, source);
}


void DirichletDoublet0Source0Case2::run()
{

//   cout << "CALCULATE GEOMETRY" << endl;
    this->mat_size = this->panels.size();
    for (int i = 0; i < this->mat_size; i++){
        this->panels[i]->calc_geo();
        this->panels[i]->sigma = this->panels[i]->n.dot(this->v_inf);
    }
//   cout << "NUMBER OF PANELS = " << this->mat_size << endl;
    this->matrix.resize(this->mat_size, this->mat_size);
    this->matrix.setZero();
    this->rhs.resize(this->mat_size);
    this->rhs.setZero();
    this->result.resize(this->mat_size);
    this->result.setZero();

//   cout << "FILLING MATRIX & RHS" << endl;

    for (Panel2*& panel_i: this->panels)
    {
        for (Panel2*& panel_j: this->panels)
        {
            double dub_influence_j_i = this->doublet_influence(*panel_i, *panel_j);
            double src_influence_j_i = this->src_influence(*panel_i, *panel_j);
            this->matrix(panel_i->nr, panel_j->nr) += dub_influence_j_i;
            this->rhs(panel_i->nr) += src_influence_j_i * panel_j->sigma;
        }
        for (PanelWakeVertex wake_point: this->wake_vertices)
        {
            double influence = this->wake_influence(*panel_i, *wake_point.panel_vector);
            int nu = wake_point.upper_operating_panel->nr;
            int nl = wake_point.lower_operating_panel->nr;
            matrix(panel_i->nr, nu) -= influence;
            matrix(panel_i->nr, nl) += influence;
        }
    }

//   cout << "SOLVE SYSTEM -> POTENTIAL" << endl;

    this->result = this->matrix.lu().solve(this->rhs);
    for (int i = 0; i < this->mat_size; i++){
        this->panels[i]->mue = this->result(i);
    }

//   cout << "COMPUTE VELOCITY" << endl;
    for (int i=0; i < this->mat_size; i++){
        this->panels[i]->compute_gradient();
        this->panels[i]->velocity += this->v_inf;
        this->panels[i]->velocity -= this->panels[i]->n * this->panels[i]->sigma;
        this->panels[i]->compute_pressure(this->v_inf);
    }
    this->sum_forces();
}


double DirichletDoublet0Source0Case2::off_body_potential(Vector2& vec)
{
    double pot = 0;
    for (Panel2*& panel: this->panels)
    {
        pot += panel->mue * doublet_2_0(vec, *panel);
        pot -= panel->sigma * source_2_0(vec, *panel);
    }
    for (PanelWakeVertex wake_point: this->wake_vertices)
    {
        pot -= vortex_2(vec, *wake_point.panel_vector) *
                    (wake_point.upper_operating_panel->mue -
                     wake_point.lower_operating_panel->mue);
    }
    pot += this->freeflow_influence(vec);
    return pot;
}

Vector2 DirichletDoublet0Source0Case2::off_body_velocity(Vector2& vec)
{
    Vector2 velocity = this->v_inf;
    //   adding velocity influence from the wing
    for (Panel2*& panel: this->panels)
    {
    //     adding the duplet contribution
        velocity += doublet_2_0_v(vec, *panel) * panel->mue;
    //     adding the source contribution
        velocity -= source_2_0_v(vec, *panel) * panel->sigma;
    }
    for (PanelWakeVertex wake_point: this->wake_vertices)
    {
        velocity -= vortex_2_v(vec, *wake_point.panel_vector) *
                    (wake_point.upper_operating_panel->mue -
                     wake_point.lower_operating_panel->mue);
    }
    return velocity;
}

/************************Linear-DIRICHLET DOUPLET-*********************/
// problem: there is a need for 2 values at the trailing edge.

DirichletDoublet1Case2::DirichletDoublet1Case2(vector< Panel2* > panels): Case2(panels)
{
    // creating the linear wake vertices
    for (PanelWakeVertex wake_point: this->wake_vertices)
    {
        PanelVector2* lower_wake_point = new PanelVector2(*wake_point.panel_vector);
        lower_wake_point->nr = this->get_all_points().size();
        if (wake_point.lower_operating_panel->points[1]->wake_vertex)
        {
            wake_point.lower_operating_panel->points[1] = lower_wake_point;
        }
        else
        {
            wake_point.lower_operating_panel->points[0] = lower_wake_point;
        }
        this->linear_wake_vertices.push_back(
            PointWakeVertex(wake_point.panel_vector, lower_wake_point)); //upper first
    }
}

vector< PanelVector2 > DirichletDoublet1Case2::collocation_points()
{
    PanelVector2 coll_point;
    vector<PanelVector2> collocation_points;
    int i = 0;
    for (Panel2*& panel: this->panels)
    {
        coll_point = panel->center + 0.533 * (panel->center - *panel->points[1]);
        coll_point.nr = i++;
        collocation_points.push_back(coll_point);
        coll_point = panel->center + 0.533 * (panel->center - *panel->points[0]);
        coll_point.nr = i++;
        collocation_points.push_back(coll_point);

    }
    return collocation_points;
}

double DirichletDoublet1Case2::surface_influence(Vector2& target, PanelVector2* source)
{
    double infl = 0;

    for (Panel2*& panel: source->panels)
    {
        if (source->wake_vertex)
        {
//             only one side of the panel will be handled, the other is corresponding to the other panel of the wake
            if (panel->points[0] != source and panel->points[1] != source)
            {
                continue;
            }
        }
        infl += doublet_2_1(target, *panel, panel->points[0] == source);
    }
    return infl;
}

double DirichletDoublet1Case2::wake_influence(Vector2& target, Vector2& source)
{
    return vortex_2(target, source);
}


double DirichletDoublet1Case2::freeflow_influence(Vector2& target)
{
    return target.dot(this->v_inf);
}


Vector2 DirichletDoublet1Case2::surface_velocity_influence(Vector2& target, PanelVector2* source)
{
    Vector2 infl(0, 0);
    for (Panel2*& panel: source->panels)
    {   if (source->wake_vertex)
            {
            // only one side of the panel will be handled, the other is corresponding to the other panel of the wake
                if (panel->points[0] != source and panel->points[1] != source)
                {
                    continue;
                }
            }
        infl += doublet_2_1_v(target, *panel, panel->points[0] == source);
    }
    return infl;
}


void DirichletDoublet1Case2::run()
{
    this->mat_size = this->get_all_points().size();
    int row_size = this->collocation_points().size();

    for (Panel2*& panel: this->panels){
        panel->calc_geo();
    }

    this->matrix.resize(row_size, this->mat_size);
    this->matrix.setZero();
    this->rhs.resize(row_size);
    this->rhs.setZero();
    this->result.resize(this->mat_size);
    this->result.setZero();
    for (PanelVector2& coll_point: this->collocation_points())
    {
        this->rhs(coll_point.nr) -= this->freeflow_influence(coll_point);
        for (PanelVector2*& panel_point: this->get_all_points())
        {
            double influence_j_i = this->surface_influence(coll_point, panel_point);
            this->matrix(coll_point.nr, panel_point->nr) = influence_j_i;
        };
        for (PointWakeVertex wake_point: this->linear_wake_vertices)
        {
            double influence = this->wake_influence(coll_point, *wake_point.upper_point );
            int nu = wake_point.upper_point->nr;
            int nl = wake_point.lower_point->nr;
            matrix(coll_point.nr, nu) -= influence;
            matrix(coll_point.nr, nl) += influence;
        }
    }
    this->result = this->matrix.colPivHouseholderQr().solve(this->rhs);   //LeastSquareSolve
    for (PanelVector2*& point: this->get_all_points())
    {
        point->potential = this->result(point->nr); // + this->freeflow_influence(*point);  // here mue = mue - phi_inf
    }
    for (Panel2*& panel: this->panels)
    {
        panel->velocity = (panel->points[0]->potential - panel->points[1]->potential) / panel->l * panel->t;
        panel->potential = (panel->points[0]->potential + panel->points[1]->potential) / 2.;
        panel->mue =  panel->potential + this->freeflow_influence(panel->center);
        panel->compute_pressure(this->v_inf);
    }
//     this->sum_forces();

}

double DirichletDoublet1Case2::off_body_potential(Vector2& target)
{
    double target_pot = 0;
    for (PanelVector2*& point: this->get_all_points())
    {
        target_pot += point->potential * this->surface_influence(target, point);
    }
    for (PointWakeVertex wake_point: this->linear_wake_vertices)
    {
        target_pot -= (wake_point.upper_point->potential - wake_point.lower_point->potential) *
                        vortex_2(target, *wake_point.upper_point);
    }
    target_pot += this->freeflow_influence(target);
    return target_pot;
}

Vector2 DirichletDoublet1Case2::off_body_velocity(Vector2& target)
{
    Vector2 target_vel(0, 0);
    for (PanelVector2*& point: this->get_all_points())
    {
        target_vel += point->potential *  this->surface_velocity_influence(target, point) / 2;
    }
    for (PointWakeVertex wake_point: this->linear_wake_vertices)
    {
        target_vel -= (wake_point.upper_point->potential - wake_point.lower_point->potential) *
                        vortex_2_v(target, *wake_point.upper_point) / 2;
    }
    target_vel += this->v_inf;
    return target_vel;
}



/************************Neumann CONSTANT SOURCE-*********************/

NeumannSource0Case2::NeumannSource0Case2(vector<Panel2*> panels):
    Case2(panels){};
    
void NeumannSource0Case2::run()
{

    cout << "CALCULATE GEOMETRY" << endl;
    this->mat_size = this->panels.size();
    for (int i = 0; i < this->mat_size; i++){
        this->panels[i]->calc_geo();
    }
//   cout << "NUMBER OF PANELS = " << this->mat_size << endl;
    this->matrix.resize(this->mat_size, this->mat_size);
    this->matrix.setZero();
    this->rhs.resize(this->mat_size);
    this->rhs.setZero();
    this->result.resize(this->mat_size);
    this->result.setZero();

//   cout << "FILLING MATRIX & RHS" << endl;

    for (Panel2*& panel_i: this->panels)
    {
        for (Panel2*& panel_j: this->panels)
        {
            this->matrix(panel_i->nr, panel_j->nr) += this->surface_influence(*panel_i, *panel_j);
            if (panel_i->nr == panel_j->nr)
            {
                this->matrix(panel_i->nr, panel_j->nr) *= -1;
            }
        }
        this->rhs[panel_i->nr] += this->freeflow_influence(*panel_i);
    }

//   cout << "SOLVE SYSTEM -> POTENTIAL" << endl;

    this->result = this->matrix.fullPivHouseholderQr().solve(this->rhs);
    for (int i = 0; i < this->mat_size; i++){
        this->panels[i]->sigma = this->result(i);
    }
    for (Panel2*& panel: this->panels)
    {
        panel->velocity = this->off_body_velocity(panel->center);
        panel->velocity += 2 * source_2_0_v(panel->center, *panel) * panel->sigma;
        panel->compute_pressure(this->v_inf);
    }
}


double NeumannSource0Case2::surface_influence(Panel2& target, Panel2& source)
{
    return source_2_0_v(target.center, source).dot(target.n);
}

double NeumannSource0Case2::freeflow_influence(Panel2& target)
{
    return this->v_inf.dot(target.n);
}

double NeumannSource0Case2::wake_influence(Panel2& target, Vector2& source)
{
    return vortex_2_v(target.center, source).dot(target.n);
}


double NeumannSource0Case2::off_body_potential(Vector2& target)
{
    double pot = 0;
    for (Panel2*& panel: this->panels)
    {
        pot -= panel->sigma * source_2_0(target, *panel);
    }
    pot += this->v_inf.dot(target);
    return pot;
}

Vector2 NeumannSource0Case2::off_body_velocity(Vector2& target)
{
    Vector2 velocity(0., 0.);
    for (Panel2*& panel: this->panels)
    {
        velocity -= panel->sigma * source_2_0_v(target, *panel);
    }
    velocity += this->v_inf;
    return velocity;
}


