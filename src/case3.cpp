#include "case3.h"

using namespace std;

AeroCoef3::AeroCoef3(Case3 c, Vector3 v_inf)
{
    this->v_inf = v_inf;
    this->cop = c.center_of_pressure;
    this->force = c.force;
    this->cL = c.cL;
    this->cD = c.cD;
    this->cS = c.cS;
    this->cR = c.cM.x();
    this->cP = c.cM.y();
    this->cY = c.cM.z();
}

double AeroCoef3::alpha()
{
    return atan(this->v_inf.z() / this->v_inf.x());
}


vector< double > AeroCoef3::as_vector()
{
    vector<double> cooeficients = 
    {
        this->alpha(),
        this->v_inf.x(), this->v_inf.y(),
        this->v_inf.z(), this->force.x(),
        this->force.y(), this->force.z(),
        this->cop.x(), this->cop.y(),
        this->cop.z(), this->cL,
        this->cD, this->cS, this->cR,
        this->cP, this->cY
    };
    return cooeficients;
}

vector< string > AeroCoef3::labels()
{
    vector<string> labels = 
    {
        "alpha", "v_inf.x" "v_inf.y" "v_inf.z",
        "force.x", "force.y", "force.z",
        "cop.x", "cop.y", "cop.z",
        "cL", "cD", "cS", "cR", "cP",
        "cY"
    };
    return labels;
}


vector<vector<double>> Polar3::as_matrix()
{
    vector<vector<double>> mat;
    for(int i = 0; i < this->size(); i++)
    {
        mat.push_back(this->at(i).as_vector());
    }
    return mat;
}

vector<string> Polar3::labels()
{
    return this->at(0).labels();
}


CoordSys::CoordSys(Vector3 v_inf, Vector3 lift_ref)
{
  this->drag_direction = v_inf;
  this->lift_direction = v_inf.cross(lift_ref).cross(v_inf);
  this->side_lift_direction = this->lift_direction.cross(this->drag_direction);
  this->drag_direction.normalize();
  this->lift_direction.normalize();
  this->side_lift_direction.normalize();
};


Case3::Case3(vector<Panel3*> panels, vector<PanelVector3*> trailing_edge)
{
    for(Panel3*& panel: panels){
        this->append_panel(panel);
    }
    for(PanelVector3*& trailing_point: trailing_edge){
        this->append_wakepoint(trailing_point);
    }
    for(Panel3*& panel: this->get_all_panels()){
        panel->set_neighbours();
    }
}

DirichletDoublet0Case3::DirichletDoublet0Case3(vector<Panel3*> panels, vector<PanelVector3*> trailing_edge):
    Case3(panels, trailing_edge){}

DirichletDoublet0Source0Case3::DirichletDoublet0Source0Case3(vector<Panel3*> panels, vector<PanelVector3*> trailing_edge):
    DirichletDoublet0Case3(panels, trailing_edge){}

void Case3::append_panel(Panel3* p)
{
    p->set_nr(this->panels.size());
    p->calc_geo();
    this->panels.push_back(p);
    for (PanelVector3*& point: p->points){
        if (point->nr == -1){
            point->nr = this->vert_size;
            this->vert_size ++;
            this->vertices.push_back(point);
        }
    }
    if (p->is_symmetric()){
        SymmetricPanel3* sym_pan = new SymmetricPanel3(p, this->symmetric_plane_n, this->symmetric_plane_p);
        sym_pan->calc_geo();
        this->sym_panels.push_back(sym_pan);
    }
}


void Case3::append_wakepoint(PanelVector3* v)
{
  this->trailing_edge.push_back(v);
  v->wake_vertex = 1;
}

void Case3::append_panel(WakePanel3* p)
{
  this->wake_panels.push_back(p);
}

vector<Panel3*> Case3::get_all_panels()
{
    vector<Panel3*> all;
    for (Panel3* panel: this->panels){
        all.push_back(panel);
    }
    for (Panel3* panel: this->sym_panels){
        all.push_back(panel);
    }
    return all;
}

vector<WakePanel3*> Case3::get_all_wake_panels()
{
    vector<WakePanel3*> all;
    for (WakePanel3* panel: this->wake_panels){
        all.push_back(panel);
    }
    for (WakePanel3* panel: this->sym_wake_panels){
        all.push_back(panel);
    }
    return all;
}


void Case3::write_to_verts()
{
    for (PanelVector3*& point: this->vertices)
    {
        point->reset_properties();
        for (Panel3*& panel: point->panels)
        {
        point->potential += panel->get_potential() / point->owner;
        point->velocity +=  panel->get_velocity() / point->owner;
        point->cp += panel->get_cp() / point->owner;
        }
    }
}


Vector3 Case3::mirror(Vector3& vec)
{
    return vec - vec.dot(this->symmetric_plane_n) * this->symmetric_plane_n * 2;
}


void Case3::sum_forces(Vector3 vinf_)
{
    int i, j, k;
    double sign, l;
    double cD_factor = 1;
    Vector3 induced_velocity, a, center, locale_force;
    CoordSys sys(vinf_, Vector3(0, 0, 1));
    this->cD = 0;
    this->cL = 0;
    this->cS = 0;
    this->cM.setZero();
    this->force.setZero();
    
    for (Panel3* panel: this->get_all_panels()){
        this->force += panel->get_force();
        this->cM += (panel->center - this->mom_ref_point).cross(panel->get_force());
    }
    a = cM.cross(force);
    a.normalize();
    l = 0;
    if (force.norm() > 0){l = cM.norm() / force.norm();}
    this->center_of_pressure = this->mom_ref_point - a * l;
    if (this->drag_calc == "trefftz" || this->drag_calc == "trailing_edge")
    {
        vector<Edge> wake_cut;
        this->force.setZero();
        
        if (this->drag_calc == "trefftz")
        {
            wake_cut = this->trefftz_cut();
        }
        
        if (this->drag_calc == "trailing_edge")
        {
            // first get the trailing edges
            wake_cut = this->get_trailing_edge();
            cD_factor = 2;
            
        }
        for (Edge edge: wake_cut)
        {
            // now sum over all the edges add to cd and ca
            // this assumes that the plane is parallel to v_inf
            center = edge.center();
            induced_velocity = this->off_body_wake_velocity_vec(center);
            locale_force = edge.vorticity * edge.diff_vector().cross(sys.drag_direction) * 2 / this->v_inf.norm();
            for (Edge te: this->get_trailing_edge())
            {   
                induced_velocity += vortex_3_0_v(center, te);
            }
            locale_force += (cD_factor * edge.vorticity *
                             edge.diff_vector().cross(induced_velocity) / 
                             pow(this->v_inf.norm(), 2)).dot(sys.drag_direction) * sys.drag_direction;
            if (edge.is_symmetric){
                locale_force += this->mirror(locale_force);
            }
            this->force += locale_force;
        }

    }
    this->cL = this->force.dot(sys.lift_direction) / this->A_ref;
    this->cD = this->force.dot(sys.drag_direction) / this->A_ref;
    this->cS = this->force.dot(sys.side_lift_direction) / this->A_ref;
    this->cM /= this->A_ref;
}


void Case3::create_wake(double length, int count, Vector3 direction){
    if (direction == Vector3(0,0,0)){
        direction = this->v_inf;
    }
    direction.normalize();
    this->wake_panels.clear();
    this->wake_streams.clear();
    this->first_wake_row.clear();
    this->sym_wake_panels.clear();
    int edge_len = this->trailing_edge.size();
    for (PanelVector3*& trailing_point: this->trailing_edge){
        vector<PanelVector3*> stream_line;
        for (int k = 0; k < count; k++){
            if (k == 0){
            stream_line.push_back(trailing_point);
        }
        else{
            PanelVector3* point = new PanelVector3(*trailing_point + direction * (length / count * k));
            stream_line.push_back(point);
        }
        }
    this->wake_streams.push_back(stream_line);
    }
    for (int i = 0; i < edge_len -1; i++)
    {
        int j = i+1;
        Edge e(this->trailing_edge[i], this->trailing_edge[j]);
        
        for (int k = 0; k < count - 1; k++)
        {
            int l = k+1;
            WakePanel3* w = new WakePanel3();
            w->append_point(this->wake_streams[i][k]);
            w->append_point(this->wake_streams[j][k]);
            w->append_point(this->wake_streams[j][l]);
            w->append_point(this->wake_streams[i][l]);
            w->set_upper_operating_panel(e.p2);
            w->set_lower_operating_panel(e.p1);
            this->wake_panels.push_back(w);
            if (k == 0){
                this->first_wake_row.push_back(w);
            }
            w->calc_geo();
            if (e.p1->is_symmetric() and e.p2->is_symmetric()){
                w->set_symmetric();
                SymmetricWakePanel3* sym_wake_pan = new SymmetricWakePanel3(w, this->symmetric_plane_n, this->symmetric_plane_p);
                sym_wake_pan->calc_geo();
                this->sym_wake_panels.push_back(sym_wake_pan);
            }
        }
    }
}

void Case3::relax_wake(int iterations, double smoothening)
{
    if (this->trailing_edge.size() == 0){return;}
    Vector3 p1;
    Vector3 p2;
    Vector3 c;
    Vector3 t;
    Vector3 v;
    double l;
    int i, j, k, iter;

    vector<Vector3> wake_offset;
    vector<Vector3> copy_wake_offset;
    wake_offset.resize(this->wake_streams.size());
    for (iter = 0; iter < iterations; iter++){
        for (i = 0; i < this->wake_streams[0].size() - 1; i++){
            for (j = 0; j < this->wake_streams.size(); j++){
                p1 = *this->wake_streams[j][i];
                p2 = *this->wake_streams[j][i + 1];
                l = (p2 - p1).norm();
                c = (p1 + p2) / 2;
                v = this->off_body_velocity_vec(p1) + this->off_body_velocity_vec(p2); v.normalize();
                t = p1 + v * l - p2;
                wake_offset[j] = t / 4;
            }
            // SMOOTHING THE WAKE
                copy_wake_offset = wake_offset;
            if (this->sym_wake_panels.size() == 0){
                for (j = 1; j < wake_offset.size() - 1; j++){
                    copy_wake_offset[j] += (wake_offset[j-1] +  wake_offset[j+1]) * smoothening / 2;
                    copy_wake_offset[j] *= 1 / (1 + smoothening);
                };
                for (j = 1; j < wake_offset.size() - 1; j++){
                    wake_offset[j] = copy_wake_offset[j];
                }
            }
            for (j = 0; j <this->wake_streams.size(); j++){
                for (k = i + 1; k < this->wake_streams[0].size(); k++){
                    *this->wake_streams[j][k] += wake_offset[j];
                }
            }
            for (WakePanel3* w: this->get_all_wake_panels()){
            w->calc_geo();
            }
        }
    }
}

vector< Edge > Case3::trefftz_cut()
{
    int i, j, k;
    vector<PanelVector3*> cut_line;
    vector<Edge> segments;
    bool found_cut = false;

    for (vector<PanelVector3*> stream: this->wake_streams){
        found_cut = false;
        for (j = 1; j < stream.size(); j++){
            Edge e(stream[j-1], stream[j]);
            double lambda = e.lambda_cut(this->trefftz_cut_pos, this->v_inf);
            if (0 <= lambda && lambda <=1){
                PanelVector3* cut_pos = new PanelVector3(e.lambda_pos(lambda));
                cut_line.push_back(cut_pos);
                found_cut = true;
                break;
            }
        }
        if (not found_cut){
            cout << "trefftz plane not aligned with wake" << endl;
            break;
        }
    }
//   creating the edges
    for (k = 1; k < cut_line.size(); k++){
        Edge e(cut_line[k - 1], cut_line[k]);
        e.vorticity = this->first_wake_row[k - 1]->get_mue();
        e.is_symmetric = this->first_wake_row[k - 1]->is_symmetric();
        segments.push_back(e);
    }
    return segments;
}

vector< Edge > Case3::get_trailing_edge()
{
    vector<Edge> segments;
    int i;
    for (i=0; i < this->trailing_edge.size() - 1; i++)
    {
        Edge e(this->trailing_edge[i], this->trailing_edge[i + 1]);
        e.vorticity = this->first_wake_row[i]->get_mue();
        e.is_symmetric = this->first_wake_row[i]->is_symmetric();
        segments.push_back(e);
    }
    return segments;
}

void Case3::calc_geo()
{
    for (Panel3* panel: this->panels){
        panel->calc_geo();
        panel->set_neighbours();
    }
}

void Case3::run()
{
    vector<Vector3> vinf_range;
    this->polars(vinf_range);
}
Polar3 Case3::polars(vector<Vector3> v_inf){}
void Case3::off_body_potential(PanelVector3& point){}
void Case3::off_body_velocity(PanelVector3& point){}
void Case3::off_body_wake_velocity(PanelVector3& point){}


Vector3 Case3::off_body_velocity_vec(Vector3& point)
{
    PanelVector3 temp_point(point);
    this->off_body_velocity(temp_point);
    return temp_point.velocity;
}

Vector3 Case3::off_body_wake_velocity_vec(Vector3& point)
{
    PanelVector3 temp_point(point);
    this->off_body_wake_velocity(temp_point);
    return temp_point.velocity;
}


vector<Vector3> Case3::flow_path(Vector3 start, double inc, int count)
{
    PanelVector3 v1(start);
    vector<Vector3> path;
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

vector< Vector3 > Case3::body_flow_path(Panel3* start, int num_traverse)
{
//     start on the center of the edge
//     calculate edge intersecting point
    vector<Vector3> flow_path;
    PanelVector3 current_point = start->center;
    flow_path.push_back(current_point);
    Panel3* current_panel = start;
    Vector3 cut_pos, plane_n, vel_dir;
    double lambda;
    Panel3* neighbour;
    for (int i = 0; i < num_traverse; i++)
    {
        vel_dir = current_panel->get_velocity();
        vel_dir.normalize();
        plane_n = current_panel->n.cross(vel_dir);
        // checking intersection on all Panel edges
        for (Edge e: current_panel->get_edges())
        {
            lambda = e.lambda_cut(current_point, plane_n);
            if (lambda >= 0 && lambda <= 1)     //--> Intersection inside the Edge
            {
                cut_pos = e.lambda_pos(lambda);
                if ((cut_pos - current_point).dot(vel_dir) > 0.000001)
                {
                    // check if neigbour is available
                    flow_path.push_back(cut_pos);
                    if (e.is_wake_edge())
                    {
                        return flow_path;
                    }
                    neighbour = current_panel->get_edge_neighbour(e);
                    if (neighbour)
                    {
                        current_point = cut_pos;
                        current_panel = neighbour;
                        break;
                                            
                    }
                    else // no neighbour is available
                    {
                        return flow_path;
                    }
                }
            }
        }
    }
    return flow_path;
}



double Case3::get_volume()
{
    double volume = 0.;
    for (Panel3* panel: this->get_all_panels())
    {
        volume += panel->center.dot(panel->n) * panel->area / 3;
    }
    return volume;
}

double Case3::get_surface_area()
{
    double area = 0.;
    for (Panel3* panel: this->get_all_panels())
    {
        area += panel->area;
    }
    return area;
}

double Case3::get_projected_area()
{
    double proj_area = 0;
    for (Panel3*& panel: this->get_all_panels())
    {
        proj_area += (panel->n.z() > 0) * panel->n.z() * panel->area;
    }
    return proj_area;
}


/************************-DIRICHLET DOUBLET0 CASE-************************/
Polar3 DirichletDoublet0Case3::polars(vector<Vector3> vinf_range)
{
    vinf_range.push_back(this->v_inf); // adding v_inf for visualization
    this->mat_size = panels.size();
    cout << "CALCULATE GEOMETRY" << endl;
    for (Panel3* panel: this->get_all_panels()){
        panel->reset_properties();
        panel->calc_geo();
        panel->set_neighbours();
    }
    cout << "NUMBER OF PANELS = " << this->mat_size << endl;
    this->matrix.resize(this->mat_size, this->mat_size);
    this->matrix.setZero();
    this->rhs.resize(this->mat_size, vinf_range.size());
    this->rhs.setZero();
    this->result.resize(this->mat_size, vinf_range.size());
    this->result.setZero();
    Polar3 polars;

    cout << "FILLING MATRIX & RHS" << endl;
    #pragma omp parallel for
    for (int i = 0; i < this->mat_size; i++){
        Panel3* panel_i = this->panels[i];
        for (Panel3* panel_j: this->get_all_panels()){
            this->surface_influence(panel_i, panel_j, vinf_range);
        }
        for (WakePanel3* panel_k: this->get_all_wake_panels()){
            this->wake_influence(panel_i, panel_k, vinf_range);
        }
        for (int h = 0; h < vinf_range.size(); h++)
        {
            this->rhs(panel_i->get_nr(), h) -= this->freeflow_influence(panel_i, vinf_range[h]);
        }
    }

    cout << "SOLVE SYSTEM -> POTENTIAL" << endl;
    Eigen::setNbThreads(omp_get_num_procs());
    this->result = this->matrix.lu().solve(this->rhs);
    cout << "COMPUTE VELOCITY & SUM FORCES" << endl;
    for (int h = 0; h < vinf_range.size() - (vinf_range.size() != 1); h++)
    {
        Vector3 vinf_ = vinf_range[h];
        for (Panel3*& panel_i: this->panels)
        {
            panel_i->reset_properties();
            panel_i->set_mue(this->result(panel_i->get_nr(), h));
        }
        for (WakePanel3*& wake_panel: this->wake_panels){
            double p_u = wake_panel->get_upper_operating_panel()->get_mue();
            double p_l = wake_panel->get_lower_operating_panel()->get_mue();
            p_u += this->freeflow_influence(wake_panel->get_upper_operating_panel(), vinf_);
            p_l += this->freeflow_influence(wake_panel->get_lower_operating_panel(), vinf_);
            wake_panel->set_mue(-p_l + p_u);
        }
        for (Panel3*& panel_i: this->panels)
        {
            panel_i->compute_gradient(vinf_);
            panel_i->add_velocity(vinf_);                                                //freestream contribution
            panel_i->add_velocity(- panel_i->n * vinf_.dot(panel_i->n));
            panel_i->set_potential(panel_i->get_mue() + this->freeflow_influence(panel_i, vinf_));          //doublet contribution
            panel_i->compute_pressure(vinf_);
        }
        this->sum_forces(vinf_);
        polars.push_back(AeroCoef3(*this, vinf_));
    }
    this->write_to_verts();
    return polars;
}



void DirichletDoublet0Case3::surface_influence(Panel3* pan_i, Panel3* pan_j, vector<Vector3> vinf_range)
{
    double influence = 0;
    if ((pan_i->center - pan_j->center).norm() > (this->farfield * pan_j->side_sum / 4)){
        doublet_3_0_n0(pan_i->center, pan_j, influence); 
    }
    else{
        doublet_3_0_vsaero(pan_i->center, pan_j, influence);
    }
    this->matrix(pan_i->get_nr(), pan_j->get_nr()) += influence;
    for (int h = 0; h < vinf_range.size(); h++)
    {
        this->rhs(pan_i->get_nr(), h) -= influence * this->freeflow_influence(pan_j, vinf_range[h]);
    }
}


void DirichletDoublet0Case3::wake_influence(Panel3* pan_i, WakePanel3* pan_k, vector<Vector3> vinf_range)
{
    double influence = 0;
    if ((pan_i->center -pan_k->center).norm() > (this->farfield * pan_k->side_sum / 4))
    {
        doublet_3_0_n0(pan_i->center, pan_k, influence); 
    }
    else
    {
        doublet_3_0_sphere(pan_i->center, pan_k, influence);
    }
    this->matrix(pan_i->get_nr(), pan_k->get_upper_operating_panel()->get_nr()) += influence;
    this->matrix(pan_i->get_nr(), pan_k->get_lower_operating_panel()->get_nr()) -= influence;
    
    for (int h = 0; h < vinf_range.size(); h++)
    {
        this->rhs(pan_i->get_nr(), h) -= influence * this->freeflow_influence(pan_k->get_upper_operating_panel(), vinf_range[h]);
        this->rhs(pan_i->get_nr(), h) += influence * this->freeflow_influence(pan_k->get_lower_operating_panel(), vinf_range[h]);

    }
}

void DirichletDoublet0Case3::off_body_potential(PanelVector3& vec){
    // adding the influence of the wing
    for (Panel3*& panel: this->get_all_panels()){
        double influence = 0;
        if ((vec - panel->center).norm() > (this->farfield * panel->side_sum / 4)){
            doublet_3_0_n0(vec, panel,  influence);
        }
        else{
            doublet_3_0_sphere(vec, panel, influence);
        }
        vec.potential += panel->get_potential() * influence;
    }
  
    // adding influence of the wake
    for (WakePanel3*& panel: this->get_all_wake_panels()){
        double dip = 0;
        if ((vec - panel->center).norm() > (this->farfield * panel->side_sum / 4)){
            doublet_3_0_n0(vec, panel, dip);
        }
        else{
        doublet_3_0_sphere(vec, panel, dip);
        }
        vec.potential += panel->get_mue() * dip;
    }

    // adding influence of the freeflow
    vec.potential += this->freeflow_influence(vec);
}


void DirichletDoublet0Case3::off_body_velocity(PanelVector3& point){
    // adding velocity influence from the wing
    for (Panel3*& panel: this->get_all_panels()){
        Vector3 velocity(0,0,0);
        doublet_3_0_vsaero_v(point, panel, velocity);
        point.velocity += velocity * panel->get_potential();
    }
  
    // adding velocity influence from the wake
    for (WakePanel3*& panel: this->get_all_wake_panels()){
        Vector3 velocity(0,0,0);
        doublet_3_0_vsaero_v(point, panel, velocity); 
        point.velocity += velocity * panel->get_mue();
    }
    point.velocity += this->v_inf;
}


void DirichletDoublet0Case3::off_body_wake_velocity(PanelVector3& point)
{
    // adding velocity influence from the wake
    for (WakePanel3*& panel: this->get_all_wake_panels()){
        Vector3 velocity(0,0,0);
        doublet_3_0_vsaero_v(point, panel, velocity); 
        point.velocity += velocity * panel->get_mue();
    }
}


double DirichletDoublet0Case3::freeflow_influence(Vector3 vec)
{
    return this->v_inf.dot(vec);
}

double DirichletDoublet0Case3::freeflow_influence(Panel3* pan_i, Vector3 vinf_)
{
    return vinf_.dot(pan_i->center);
}

double DirichletDoublet0Case3::freeflow_influence(Panel3* pan_i)
{
    return this->freeflow_influence(pan_i, this->v_inf);
}



/************************-DIRICHLET-DOUBLE0-SOURCE0-************************/

Polar3 DirichletDoublet0Source0Case3::polars(vector<Vector3> vinf_range)
{
    vinf_range.push_back(this->v_inf); // adding v_inf for visualization
    this->mat_size = panels.size();
    cout << "CALCULATE GEOMETRY" << endl;
    for (Panel3* panel: this->get_all_panels()){
        panel->reset_properties();
        panel->set_neighbours();
    }
    
    cout << "NUMBER OF PANELS = " << this->mat_size << endl;
    this->matrix.resize(this->mat_size, this->mat_size);
    this->matrix.setZero();
    this->rhs.resize(this->mat_size, vinf_range.size());
    this->rhs.setZero();
    this->result.resize(this->mat_size, vinf_range.size());
    this->result.setZero();
    Polar3 polars;
    
    cout << "FILLING MATRIX & RHS" << endl;
    #pragma omp parallel for
    for (int i = 0; i < this->mat_size; i++){
        Panel3* panel_i = this->panels[i];
        for (Panel3* panel_j: this->get_all_panels()){
            this->surface_influence(panel_i, panel_j, vinf_range);
        }
        for (WakePanel3* panel_k: this->get_all_wake_panels()){
            this->wake_influence(panel_i, panel_k, vinf_range);
        }
    }
    
    cout << "SOLVE SYSTEM -> POTENTIAL" << endl;
    Eigen::setNbThreads(omp_get_num_procs());
    this->result = this->matrix.lu().solve(this->rhs);

    cout << "COMPUTE VELOCITY" << endl;
    cout << "SUM FORCES" << endl;

    for (int h = 0; h < vinf_range.size() - (vinf_range.size() != 1); h++)
    {
        Vector3 vinf_ = vinf_range[h];
        for (Panel3*& panel: this->panels)
        {
            panel->set_mue(this->result(panel->get_nr(), h));
            panel->set_sigma(vinf_.dot(panel->n));
        }
        for (WakePanel3* wake_panel: this->get_all_wake_panels()){
            double p_u = wake_panel->get_upper_operating_panel()->get_mue();
            double p_l = wake_panel->get_lower_operating_panel()->get_mue();
            p_u += this->freeflow_influence(wake_panel->get_upper_operating_panel(), vinf_);
            p_l += this->freeflow_influence(wake_panel->get_lower_operating_panel(), vinf_);
            wake_panel->set_mue(-p_l + p_u);
        }
        for (int i = 0; i < this->mat_size; i++)
        {
            this->panels[i]->compute_gradient(vinf_);
            this->panels[i]->add_velocity(vinf_);
            this->panels[i]->add_velocity(- this->panels[i]->n * this->panels[i]->get_sigma());
            this->panels[i]->set_potential(this->panels[i]->get_mue() + this->freeflow_influence(this->panels[i], vinf_));
            this->panels[i]->compute_pressure(vinf_);
        }
        this->sum_forces(vinf_);
        polars.push_back(AeroCoef3(*this, vinf_));
    }
    this->write_to_verts();
    return polars;
}

void DirichletDoublet0Source0Case3::surface_influence(Panel3* pan_i, Panel3* pan_j, vector<Vector3> vinf_range)
{
    double dip = 0;
    double src = 0;
    if ((pan_i->center - pan_j->center).norm() > (this->farfield * pan_j->side_sum / 4))
    {
        doublet_src_3_0_n0(pan_i->center, pan_j, dip, src); 
    }
    else
    {
        doublet_src_3_0_vsaero(pan_i->center, pan_j, dip, src);
    }
    this->matrix(pan_i->get_nr(), pan_j->get_nr()) += dip;
    for (int h = 0; h < vinf_range.size(); h++)
    {
        this->rhs(pan_i->get_nr(), h) += pan_j->n.dot(vinf_range[h]) * src;
    }
}

void DirichletDoublet0Source0Case3::off_body_velocity(PanelVector3& vec){

    //   adding velocity influence from the wing
    for (Panel3* panel: this->get_all_panels())
    {
    //     adding the doublet and source contribution
        vec.velocity += doublet_src_3_0_vsaero_v(vec, panel, panel->get_mue(), panel->get_sigma());
    }
    
    // adding velocity influence from the wake
    for (WakePanel3* wake_panel: this->get_all_wake_panels())
    {
        Vector3 velocity(0,0,0);
        vec.velocity += doublet_3_0_vsaero_v(vec, wake_panel) * wake_panel->get_mue();
    }
    vec.velocity += this->v_inf;
}

void DirichletDoublet0Source0Case3::off_body_potential(PanelVector3& vec)
{
    //   adding the influence of the wing
    for (Panel3* panel: this->get_all_panels())
    {
        double dip = 0;
        double src = 0;
        if ((vec - panel->center).norm() > (this->farfield * panel->side_sum / 4)){
            doublet_src_3_0_n0(vec, panel, dip, src);
        }
        else{
            doublet_src_3_0_vsaero(vec, panel, dip, src);
        }
        vec.potential += panel->get_mue() * dip;
        vec.potential -= panel->get_sigma() * src;
    }

    //   adding influence of the wake
    for (WakePanel3* panel: this->get_all_wake_panels())
    {
        double dip = 0;
        if ((vec - panel->center).norm() > (this->farfield * panel->side_sum / 4)){
            doublet_3_0_n0(vec, panel, dip);
        }
        else{
            doublet_3_0_sphere(vec, panel, dip);
        }
        vec.potential += panel->get_potential() * dip;
    }
    //    adding influence of the freeflow
    vec.potential += this->freeflow_influence(vec);
}
