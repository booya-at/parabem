#include "ppm_wrap.h"
#include <math.h>

using namespace std;

// ###########################-PY_CASE_WRAPPER-##################################

void py_case::set_wake_points(Case3& c, boost::python::list& points)
{
    for (int i = 0; i < boost::python::len(points); i++){
        c.append_wakepoint(boost::python::extract<PanelVector3*>(points[i]));
    }
}

Polar3 py_case::polar3(Case3& c, boost::python::list& vinf_range)
{
    vector<Vector3> vector_vinf_range;
    for (int i = 0; i < boost::python::len(vinf_range); i++){
        vector_vinf_range.push_back(boost::python::extract<Vector3>(vinf_range[i]));
    }
    return c.polars(vector_vinf_range);
}


// ###########################-PY_EIGEN_WRAPPER-##################################

boost::python::list py_eigen::mat_to_list(Eigen::MatrixXd mat){
    boost::python::list l;
    for (int i = 0; i < mat.rows() ; i++){
        boost::python::list li;
        for (int j = 0; j < mat.cols(); j++){
            li.append<double>(mat(i, j));
        }
        l.append<boost::python::list>(li);
    }
    return l;
}


boost::python::list py_eigen::vec_to_list(Eigen::VectorXd vec){
    boost::python::list l;
    for (int i = 0; i < vec.rows(); i++){
        l.append<double>(vec[i]);
    }
    return l;
}



// ############################-PY_ELEMENT_INFLUENCE-#############################

double py_element_influence::wrap_doublet_3_0_vsaero(Vector3& target, Panel3* source){
    double a = 0.;
    doublet_3_0_vsaero(target, source, a);
    return a;
}

double py_element_influence::wrap_doublet_3_0_sphere(Vector3& target, Panel3* source){
    double a = 0.;
    doublet_3_0_sphere(target, source, a);
    return a;
}

double py_element_influence::wrap_doublet_3_0_n0(Vector3& target, Panel3* source){
    double a = 0.;
    doublet_3_0_n0(target, source, a);
    return a;
}

Vector3 py_element_influence::wrap_doublet_3_0_vsaero_v(Vector3& target, Panel3* source)
{
    Vector3 v (0,0,0);
    doublet_3_0_vsaero_v(target, source, v);
    return v;
}

boost::python::tuple py_element_influence::wrap_doublet_src_3_0_n0(Vector3& target, Panel3* source)
{
    double a = 0.;
    double b = 0.;
    doublet_src_3_0_n0(target, source, a, b);
    return boost::python::make_tuple(a, b);
}


boost::python::tuple py_element_influence::wrap_doublet_src_3_0_vsaero(Vector3& target, Panel3* source)
{
    double a = 0.;
    double b = 0.;
    doublet_src_3_0_vsaero(target, source, a, b);
    return boost::python::make_tuple(a, b);
}

Vector3 py_element_influence::wrap_vortex_3_0_edge_v(Vector3& target, Edge& e)
{
    return vortex_3_0_v(target, e);
}

Vector3 py_element_influence::wrap_vortex_3_0_v(Vector3& target, Vector3& p1, Vector3& p2)
{
    return vortex_3_0_v(target, p1, p2);
}

Vector3 py_element_influence::wrap_source_3_0_vsaero_v(Vector3& target, Panel3* source)
{
    return src_3_0_vsaero_v(target, source, 1.);
}



// #################################### PY_Vector2 ##########################################


boost::shared_ptr< Vector2 > py_vector_2::init_from_Vector3(Vector3 v)
{
    Vector2* a = new Vector2(v.x(), v.y());
    return boost::shared_ptr<Vector2>(a);
}

boost::shared_ptr< Vector2 > py_vector_2::init_from_Vector2(Vector2 v)
{
    Vector2* a = new Vector2(v.x(), v.y());
    return boost::shared_ptr<Vector2>(a);
}



double py_vector_2::getitem(Vector2& v, int i)
{
    if (i < 2){
        return v(i);
    }
    else{
        throw out_of_range("list index out of range");
    }
}

void py_vector_2::setitem(Vector2& v, int i, double value)
{
    if (0 <= i && i <= 1){
        v(i) = value;
    }
    else{
        throw out_of_range("list index out of range");
    }
}

double py_vector_2::get_x(Vector2& v)
{
    double out;
    out = v.x();
    return out;
}

double py_vector_2::get_y(Vector2& v)
{
    double out;
    out = v.y();
    return out;
}

void py_vector_2::set_x(Vector2& v, double x)
{
    v.x() = x; 		//copy
}

void py_vector_2::set_y(Vector2& v, double y)
{
    v.y() = y;		//copy
}


string py_vector_2::repr(Vector2& v)
{
    ostringstream out_str;
    out_str << "Eigen::Vector2(" << v.x() << ", " << v.y() << ")";
    return out_str.str();
}

double py_vector_2::dot(Vector2& v1, Vector2& v2){return v1.dot(v2);}
double py_vector_2::abs(Vector2& v){return v.norm();}
Vector2 py_vector_2::add(Vector2& v1, Vector2& v2){return v1 + v2;}
Vector2 py_vector_2::sub(Vector2& v1, Vector2& v2){return v1 - v2;}
Vector2 py_vector_2::mul(Vector2& v, double d){return v * d;}
Vector2 py_vector_2::div(Vector2& v, double d){return v * (1/d);}
Vector2 py_vector_2::neg(Vector2& v){return -v;}
Vector2 py_vector_2::pos(Vector2& v){return v;}


// #################################### PY_Vector3 ##########################################

boost::shared_ptr< Vector3 > py_vector_3::init_from_Vector2(Vector2 v)
{
    Vector3* a = new Vector3(v.x(), v.y(), 0);
    return boost::shared_ptr<Vector3>(a);
}

boost::shared_ptr< Vector3 > py_vector_3::init_from_Vector3(Vector3 v)
{
    Vector3* a = new Vector3(v.x(), v.y(), v.z());
    return boost::shared_ptr<Vector3>(a);
}


double py_vector_3::getitem(Vector3& v, int i)
{
    if (0 <= i && i <= 2){
        return v(i);
    }
    else{
        throw out_of_range("list index out of range");
    }
}

void py_vector_3::setitem(Vector3& v, int i, double value)
{
    if (0 <= i && i <= 2){
        v(i) = value;
    }
    else{
        throw out_of_range("list index out of range");
    }
}

double py_vector_3::get_x(Vector3& v)
{
    double out;
    out = v.x();
    return out;
}

double py_vector_3::get_y(Vector3& v)
{
    double out;
    out = v.y();
    return out;
}

double py_vector_3::get_z(Vector3& v)
{
    double out;
    out = v.z();
    return out;
}


void py_vector_3::set_x(Vector3& v, double x)
{
    v.x() = x; 		//copy
}

void py_vector_3::set_y(Vector3& v, double y)
{
    v.y() = y;		//copy
}

void py_vector_3::set_z(Vector3& v, double z)
{
    v.z() = z;		//copy
}

string py_vector_3::repr(Vector3& v)
{
    ostringstream out_str;
    out_str << "Eigen::Vector3(" << v.x() << ", " << v.y() << ", " << v.z() << ")";
    return out_str.str();
}

double py_vector_3::abs(Vector3& v)
{
    return v.norm();
}


Vector3 py_vector_3::add(Vector3& v1, Vector3& v2){return v1 + v2;}
Vector3 py_vector_3::sub(Vector3& v1, Vector3& v2){return v1 - v2;}
Vector3 py_vector_3::mul(Vector3& v, double d){return v * d;}
Vector3 py_vector_3::div(Vector3& v, double d){return v * (1/d);}
Vector3 py_vector_3::neg(Vector3& v){return -v;}
Vector3 py_vector_3::pos(Vector3& v){return v;}
