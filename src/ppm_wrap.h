/**
    Purpose: wrapper functions to get a (cleaner) python module
    @author L Lechner
*/


#ifndef ppm_wrap_H
#define ppm_wrap_H

#include <boost/python.hpp>
#include <boost/python/errors.hpp>
#include <iostream>
#include <vector>

#include "vector.h"
#include "string.h"
#include "panel3.h"
#include "panel2.h"
#include "case3.h"
#include "case2.h"
#include "element_influence.h"


template <typename T>
boost::python::list to_python_list(vector<T> vec)
{
    boost::python::list l;
    for (T value: vec)
        {
            l.append<T>(value);
        }
    return l;
}

template <typename T>
boost::python::list mat_to_python_list(vector<vector<T>> mat)
{
    boost::python::list l;
    for (vector<T> vec: mat){ 
        boost::python::list li;
        for (T value: vec)
        {
            li.append<T>(value);
        }
        l.append<boost::python::list>(li); 
    }
    return l;
}

class py_eigen{
public:
    static boost::python::list mat_to_list(Eigen::MatrixXd mat);
    static boost::python::list vec_to_list(Eigen::VectorXd vec);
};

class py_panel{
public:
    template <typename t_panel, typename t_vector> static boost::shared_ptr<t_panel> init_Panel(boost::python::list& points);
    template <typename t_panel, typename t_vector> static boost::python::list panel_vectors(t_panel& p);
};

class py_polar{
public:
    template <typename t_polar> 
    static boost::python::list as_matrix(t_polar& polar)
        {return mat_to_python_list(polar.as_matrix());};
    template <typename t_polar> 
    static boost::python::list labels(t_polar& polar)
        {return to_python_list(polar.labels());};
};

class py_case{
public:
    template <typename t_case> static boost::shared_ptr<t_case> init_from_pans(
        boost::python::list& panels);
    template <typename t_case> static boost::shared_ptr<t_case> init_from_pans_and_te(
        boost::python::list& panels,
        boost::python::list& points);
    template <typename t_case> static boost::shared_ptr<t_case> init_from_pans_2(
        boost::python::list& panels);
    template <typename t_case, typename t_panel> static boost::python::list case_panels(t_case& c);
    template <typename t_case, typename t_panel> static boost::python::list case_wake_panels(t_case& c);
    template <typename t_case, typename t_panel> static void set_panels(t_case& c, boost::python::list& panels);
    static void set_wake_points(Case3& c, boost::python::list& points);
    static Polar3 polar3(Case3& c, boost::python::list& vinf_range);
};


class py_element_influence{
public:
    static double wrap_doublet_3_0_vsaero(Vector3& target, Panel3* source);
    static double wrap_doublet_3_0_sphere(Vector3& target, Panel3* source);
    static double wrap_doublet_3_0_n0(Vector3& target, Panel3* source);
    static boost::python::tuple wrap_doublet_src_3_0_vsaero(Vector3& target, Panel3* source);
    static boost::python::tuple wrap_doublet_src_3_0_n0(Vector3& target, Panel3* source);
    static Vector3 wrap_doublet_3_0_vsaero_v(Vector3& target, Panel3* source);
    static Vector3 wrap_source_3_0_vsaero_v(Vector3& target, Panel3*source);
    static Vector3 wrap_vortex_3_0_v(Vector3& target, Vector3& p1, Vector3& p2);
    static Vector3 wrap_vortex_3_0_edge_v(Vector3& target, Edge& e);
};


class py_vector_2{
public:
    template <class TVECTOR2>  static boost::shared_ptr<TVECTOR2> init_from_double(double x, double y);
    template <class TVECTOR2>  static boost::shared_ptr<TVECTOR2> init_from_list(boost::python::list l);
    static boost::shared_ptr<Vector2> init_from_Vector3(Vector3);
    static boost::shared_ptr<Vector2> init_from_Vector2(Vector2);
    static double getitem(Vector2& v, int i);
    static void setitem(Vector2& v, int i, double value);
    static double get_x(Vector2& v);
    static double get_y(Vector2& v);
    static void set_x(Vector2& v, double x);
    static void set_y(Vector2& v, double y);
    static string repr(Vector2& v); 
    static double dot(Vector2& v1, Vector2& v2);
    static Vector2 add(Vector2& v1, Vector2& v2);
    static Vector2 sub(Vector2& v1, Vector2& v2);
    static Vector2 mul(Vector2& v, double d);
    static Vector2 div(Vector2& v, double d);
    static double abs(Vector2& v);
    static Vector2 neg(Vector2& v);
    static Vector2 pos(Vector2& v);
};


class py_vector_3{
public:
    template <class TVECTOR3>  static boost::shared_ptr<TVECTOR3> init_from_double(double x, double y, double z);
    template <class TVECTOR3>  static boost::shared_ptr<TVECTOR3> init_from_list(boost::python::list l);
    static boost::shared_ptr<Vector3> init_from_Vector2(Vector2);
    static boost::shared_ptr<Vector3> init_from_Vector3(Vector3);
    static double getitem(Vector3& v, int i);
    static void setitem(Vector3& v, int i, double value);
    static double get_x(Vector3& v);
    static double get_y(Vector3& v);
    static double get_z(Vector3& v);
    static void set_x(Vector3& v, double x);
    static void set_y(Vector3& v, double y);
    static void set_z(Vector3& v, double z);
    static string repr(Vector3& v); 
    static Vector3 add(Vector3& v1, Vector3& v2);
    static Vector3 sub(Vector3& v1, Vector3& v2);
    static Vector3 mul(Vector3& v, double d);
    static Vector3 div(Vector3& v, double d);
    static Vector3 cross(Vector3& v1, Vector3& v2){return v1.cross(v2);}
    static double dot(Vector3& v1, Vector3& v2){return v1.dot(v2);}
    static double abs(Vector3& v);
    static Vector3 neg(Vector3& v);
    static Vector3 pos(Vector3& v);
};



template<typename t_panel, typename t_vector>
boost::shared_ptr<t_panel> py_panel::init_Panel(boost::python::list& points){
    t_panel* p = new t_panel();
    for (int i = 0; i < boost::python::len(points); i++){
        p->append_point(boost::python::extract<t_vector*>(points[i]));
    }
    p->calc_geo();
    return boost::shared_ptr<t_panel>(p);
}


template<typename t_panel, typename t_vector>
boost::python::list py_panel::panel_vectors(t_panel& p){
    boost::python::list l;
    for (t_vector* point: p.points){
        l.append<t_vector>(*point);
    }
    return l;
}

template <typename t_case>
boost::shared_ptr<t_case> py_case::init_from_pans_and_te(
    boost::python::list& panels, 
    boost::python::list& points)
{
    vector<Panel3*> pans;
    vector<PanelVector3*> trailing_edge;
    for (int i = 0; i < boost::python::len(panels); i++){
        // add if wakePanel ...
        pans.push_back(boost::python::extract<Panel3*>(panels[i]));
    }
    for (int i = 0; i < boost::python::len(points); i++){
        trailing_edge.push_back(boost::python::extract<PanelVector3*>(points[i]));
    }
    t_case* c =  new t_case(pans, trailing_edge);
    return boost::shared_ptr<t_case>(c);
}

template <typename t_case>
boost::shared_ptr<t_case> py_case::init_from_pans(boost::python::list& panels)
{
    boost::python::list l;
    return py_case::init_from_pans_and_te<t_case>(panels, l);
}

template <typename t_case>
boost::shared_ptr<t_case> py_case::init_from_pans_2(boost::python::list& panels)
{   
    vector<Panel2*> pans;
    for (int i = 0; i < boost::python::len(panels); i++){
        // add if wakePanel ...
        pans.push_back(boost::python::extract<Panel2*>(panels[i]));
    }
    t_case* c =  new t_case(pans);
    return boost::shared_ptr<t_case>(c);
}

template<typename t_case, typename t_panel>
void py_case::set_panels(t_case& c, boost::python::list& panels){
    c.panels.resize(0);
    for (int i = 0; i < boost::python::len(panels); i++){
        // add if wakePanel ...
        c.append_panel(boost::python::extract<t_panel*>(panels[i]));
    }
}

template <typename t_case, typename t_panel>
boost::python::list py_case::case_panels(t_case& c){
    boost::python::list l;
    vector<t_panel*> panels = c.panels;
    for (t_panel* panel: panels){
        l.append<t_panel>(*panel);
    }
    return l;
}

template <typename t_case, typename t_panel>
boost::python::list py_case::case_wake_panels(t_case& c){
    boost::python::list l;
    for (t_panel* pan: c.get_all_wake_panels()){
            l.append<t_panel>(*pan);
    }
    return l;
}


template <class TVECTOR2>
boost::shared_ptr<TVECTOR2> py_vector_2::init_from_list(boost::python::list l){
    TVECTOR2* v = new TVECTOR2;
    v->x() = boost::python::extract<double>(l[0]);
    v->y() = boost::python::extract<double>(l[1]);
    return boost::shared_ptr<TVECTOR2>(v);
}

template <class TVECTOR2>
boost::shared_ptr<TVECTOR2> py_vector_2::init_from_double(double x, double y){
    TVECTOR2* v = new TVECTOR2;
    v->x() = x;
    v->y() = y;
    return boost::shared_ptr<TVECTOR2>(v);
}


template <class TVECTOR3>
boost::shared_ptr<TVECTOR3> py_vector_3::init_from_list(boost::python::list l){
    TVECTOR3* v = new TVECTOR3;
    v->x() = boost::python::extract<double>(l[0]);
    v->y() = boost::python::extract<double>(l[1]);
    v->z() = boost::python::extract<double>(l[2]);
    return boost::shared_ptr<TVECTOR3>(v);
}

template <class TVECTOR3>
boost::shared_ptr<TVECTOR3> py_vector_3::init_from_double(double x, double y, double z){
    TVECTOR3* v = new TVECTOR3;
    v->x() = x;
    v->y() = y;
    v->z() = z;
    return boost::shared_ptr<TVECTOR3>(v);
}

#endif // vector_H