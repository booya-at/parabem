#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <iostream>
#include "string.h"
#include "Eigen/Core"

#include <vector>
#include <tuple>

#include "vector.h"
#include "panel3.h"
#include "panel2.h"
#include "case3.h"
#include "case2.h"
#include "element_influence.h"
#include "lifting_line.h"


using std::cout;
using std::endl;

namespace py = pybind11;

void init_eigen(py::module &m);

template<typename Vector>
void from_python_list(Vector& instance, vector<double> v){
    if (v.size() != instance.size())
        throw py::index_error("number of elements given, not match Vector size");
    new (&instance) Vector();
    int i = 0;
    for (double vi: v){
        instance[i] = vi;
        i++;
    }
}

double wrap_doublet_3_0_vsaero(Vector3& target, Panel3* source){
    double output=0; doublet_3_0_vsaero(target, source, output); return output;}

double wrap_doublet_3_0_sphere(Vector3& target, Panel3* source){
    double output=0; doublet_3_0_sphere(target, source, output); return output;}

double wrap_doublet_3_0_n0(Vector3& target, Panel3* source){
    double output=0; doublet_3_0_n0(target, source, output); return output;}

std::tuple<double, double> wrap_doublet_src_3_0_vsaero(Vector3& target, Panel3* source){
    double src_infl = 0; double dip_infl = 0;
    doublet_src_3_0_vsaero(target, source, dip_infl, src_infl);
    return std::make_tuple(dip_infl, src_infl);}

std::tuple<double, double> wrap_doublet_src_3_0_n0(Vector3& target, Panel3* source){
    double src_infl = 0; double dip_infl = 0;
    doublet_src_3_0_n0(target, source, dip_infl, src_infl);
    return std::make_tuple(dip_infl, src_infl);}

Vector3 wrap_doublet_3_0_vsaero_v(Vector3& target, Panel3* source){
    Vector3 output=Vector3(0,0,0); doublet_3_0_vsaero_v(target, source, output); return output;}

Vector3 wrap_source_3_0_vsaero_v(Vector3& target, Panel3*source){
    return src_3_0_vsaero_v(target, source);}

Vector3 wrap_vortex_3_0_v(Vector3& target, Vector3& p1, Vector3& p2){
    return vortex_3_0_v(target, p1, p2);}

Vector3 wrap_vortex_3_0_edge_v(Vector3& target, Edge& e){
    return vortex_3_0_v(target, e);}


PYBIND11_PLUGIN(_ppm) {
    py::module m("_ppm", "pybind11 example plugin");
    init_eigen(m);
    
    py::class_<PanelVector2>(m, "PanelVector2", py::base<Vector2>())
        .def(py::init<double, double>())
        .def("__init__", &from_python_list<PanelVector2>)
        .def_readwrite("wake_vertex", &PanelVector2::wake_vertex, "if true this point shades a wake")
        .def_readwrite("potential", &PanelVector2::potential, "potential at this point (not all methodes set this value)")
        .def_readwrite("velocity", &PanelVector2::velocity, "velocity at this point (not all methodes set this value)")
        .def_readonly("nr", &PanelVector2::nr, "nr of this point");
        
    py::class_<Panel2>(m, "Panel2")
        //docs: "Panel2 represents a straigth line which is used to approximate 2d geometry\n\
        //            constructor: Panel2([p1, p2])\n\n\
        //            p1, p2 -> PanelVector2 which have to be stored somewhere in python (list, variable,...)"
        .def(py::init<PanelVector2*, PanelVector2*>())
        .def(py::init<vector<PanelVector2*>>())
        .def_readonly("l", &Panel2::l, "length of panel")
        .def_readonly("t", &Panel2::t, "panel direction")
        .def_readonly("n", &Panel2::n, "normal direction")
        .def_readonly("center", &Panel2::center, "panel center")
        .def_property_readonly("points", &Panel2::get_points)
        .def_readwrite("potential", &Panel2::potential, "value is set from case")
        .def_readwrite("mue", &Panel2::mue, "value is set from case")
        .def_readwrite("sigma", &Panel2::sigma, "value is set from case")
        .def_readonly("cp", &Panel2::cp, "value is set from case")
        .def_readonly("velocity", &Panel2::velocity, "value is set from case");
        
    py::class_<Case2>(m, "Case2")
        .def_readonly("panels", &Case2::panels, "panels of case")
        .def_property_readonly("points", &Case2::get_all_points, "all points")
        .def_readwrite("v_inf", &Case2::v_inf, "direction of parallell flow")
        .def_readwrite("A_ref", &Case2::A_ref, "reference Area")
        .def_readwrite("mom_ref_point", &Case2::mom_ref_point, "point to sum moment")
        .def_readonly("matrix", &Case2::matrix, "resulting matrix")
        .def_readonly("result", &Case2::result, "solution of the linear sysem")
        .def_readonly("rhs", &Case2::rhs, "the right hand side of the equation")
        .def_readonly("mat_size", &Case2::mat_size, "size of the matrix")
        .def_readonly("force", &Case2::force, "sum of all panel-forces")
        .def_readonly("center_of_pressure", &Case2::center_of_pressure, "center of pressure")
        .def_readonly("cl", &Case2::ca, "force normal to vinf")
        .def_readonly("cm", &Case2::moment, "resulting moment")
        .def("run", &Case2::run, "run the case")
        .def("off_body_velocity", &Case2::off_body_velocity, "velocity field")
        .def("off_body_potential", &Case2::off_body_potential, "potetnial field");
        
    py::class_<DirichletDoublet0Case2>(m, "DirichletDoublet0Case2", py::base<Case2>())
        .def(py::init<vector<Panel2*>>());
        
    py::class_<NeumannDoublet0Case2>(m, "NeumannDoublet0Case2", py::base<Case2>())
        .def(py::init<vector<Panel2*>>());
        
    py::class_<NeumannSource0Case2>(m, "NeumannSource0Case2", py::base<Case2>())
        .def(py::init<vector<Panel2*>>());
    
    py::class_<DirichletDoublet0Source0Case2>(m, "DirichletDoublet0Source0Case2", py::base<DirichletDoublet0Case2>())
        .def(py::init<vector<Panel2*>>());
        
    py::class_<DirichletDoublet1Case2>(m, "DirichletDoublet1Case2", py::base<Case2>())
        .def(py::init<vector<Panel2*>>());

    py::class_<PanelVector3>(m, "PanelVector3", py::base<Vector3>())
    //           "Panel3 represents a straigth line which is used to approximate 2d geometry\n\
    //           constructor: Panel2([p1, p2, ...])\n\n\
    //           p1, p2... -> PanelVector2 which have to be stored somewhere in python (list, variable,...)")
        .def(py::init<double, double, double>())
        .def("__init__", &from_python_list<Vector3>)
        .def_readonly("nr", &PanelVector3::nr)
        .def_readonly("potential", &PanelVector3::potential)
        .def_readonly("cp", &PanelVector3::cp)
        .def_readonly("vorticity", &PanelVector3::vorticity)
        .def_readonly("velocity", &PanelVector3::velocity);

    py::class_<Panel3>(m, "Panel3")
        .def(py::init<vector<PanelVector3*>>())
        .def(py::init<>())
        .def("append_point", &Panel3::append_point)
        .def_readonly("center", &Panel3::center)
        .def_readonly("m", &Panel3::m)
        .def_readonly("n", &Panel3::n)
        .def_readonly("l", &Panel3::l)
        .def_readonly("area", &Panel3::area)
        .def_property("potential", &Panel3::get_potential, &Panel3::set_potential)
        .def_property("sigma", &Panel3::get_sigma, &Panel3::set_sigma)
        .def_property_readonly("mue", &Panel3::get_mue)
        .def_property_readonly("velocity", &Panel3::get_velocity)
        .def_property_readonly("cp", &Panel3::get_cp)
        .def_readonly("points", &Panel3::points)
        .def("get_points", &Panel3::get_points, "", py::return_value_policy::copy)
        .def("set_symmetric", &Panel3::set_symmetric);

    py::class_<Edge>(m, "Edge")
        .def_readonly("v1", &Edge::v1)
        .def_readonly("v2", &Edge::v2)
        .def_readonly("vorticity", &Edge::vorticity)
        .def_property_readonly("center", &Edge::center);

// workaround because subclassing messed things up (double frees, ...)
    py::class_<WakePanel3>(m, "WakePanel3")
        .def("get_points", [](WakePanel3& pan){return pan.get_points();}, "", py::return_value_policy::copy)
        .def_property_readonly("n", [](WakePanel3& pan){return pan.n;})
        .def_property("potential", [](WakePanel3& pan){return pan.get_potential();}, 
                                   [](WakePanel3& pan, double pot){pan.set_potential(pot);})
        .def_property_readonly("upper_operating_Panel", &WakePanel3::get_upper_operating_panel)
        .def_property_readonly("lower_operating_Panel", &WakePanel3::get_lower_operating_panel);
    
    py::class_<AeroCoef3>(m, "AeroCoef3")
        .def_property_readonly("alpha", &AeroCoef3::alpha)
        .def_readonly("v_inf", &AeroCoef3::v_inf)
        .def_readonly("force", &AeroCoef3::force)
        .def_readonly("cop", &AeroCoef3::cop)
        .def_readonly("cL", &AeroCoef3::cL)
        .def_readonly("cD", &AeroCoef3::cD)
        .def_readonly("cS", &AeroCoef3::cS)
        .def_readonly("cR", &AeroCoef3::cR)
        .def_readonly("cP", &AeroCoef3::cP)
        .def_readonly("cY", &AeroCoef3::cY);

    py::class_<Polar3>(m, "Polar3")
        .def_readonly("values", &Polar3::values)
        .def_property_readonly("labels", &Polar3::labels)
        .def_property_readonly("as_matrix", &Polar3::as_matrix);

    py::class_<Case3>(m, "Case3")
        .def_readonly("panels", &Case3::panels)
        .def_property_readonly("trailing_edge", &Case3::get_trailing_edge)
        .def_readwrite("v_inf", &Case3::v_inf)
        .def_readwrite("farfield", &Case3::farfield)
        .def_readonly("mat_size", &Case3::mat_size)
        .def_readonly("matrix", &Case3::matrix)
        .def_readonly("result", &Case3::result)
        .def_readonly("rhs", &Case3::rhs)
        .def_readwrite("drag_calc", &Case3::drag_calc)
        .def_readwrite("trefftz_cut_pos", &Case3::trefftz_cut_pos)
        .def("trefftz_cut", &Case3::trefftz_cut)
        .def("create_wake", (void (Case3::*)(double, int, Vector3&)) &Case3::create_wake, "create wake",
             py::arg("length") = 100., py::arg("count") = 10, py::arg("direction") = Vector3(0, 0, 0))
        .def("relax_wake", &Case3::relax_wake, "relax the wake", py::arg("iterations")=1, py::arg("smoothening")=1)
        .def_readwrite("symmetric_plane_n", &Case3::symmetric_plane_n)
        .def_readwrite("symmetric_plane_p", &Case3::symmetric_plane_p)
        .def_readonly("vertices", &Case3::vertices)
        .def_readwrite("A_ref", &Case3::A_ref)
        .def_readwrite("mom_ref_point", &Case3::mom_ref_point)
        .def_readwrite("lift_ref", &Case3::lift_ref)
        .def_readonly("force", &Case3::force)
        .def_readonly("center_of_pressure", &Case3::center_of_pressure)
        .def_readonly("cL", &Case3::cL, "scalar lift cooeficient: cL = L * 2 / (rho * u**2 * A_ref)")
        .def_readonly("cD", &Case3::cD, "scalar drag cooeficient: cD = D * 2 / (rho * u**2 * A_ref)")        
        .def_readonly("cS", &Case3::cS, "scalar side lift cooeficient: cS = S * 2 / (rho * u**2 * A_ref)")
        .def_readonly("cM", &Case3::cM, "Vector3 moment cooeficients: cM = M * 2 / (rho * u**2 * A_ref * 1)")
        .def_readonly("wake_panels", &Case3::wake_panels)
        .def("get_surface_area", &Case3::get_surface_area)
        .def("get_volume", &Case3::get_volume)
        .def("get_projected_area", &Case3::get_projected_area)
        .def("off_body_velocity", &Case3::off_body_velocity)
        .def("off_body_potential",&Case3::off_body_potential)
        .def("flow_path", &Case3::flow_path)
        .def("body_flow_path", &Case3::body_flow_path)
        .def("run", &Case3::run)
        .def("polars", &Case3::polars);
        
    py::class_<DirichletDoublet0Case3>(m, "DirichletDoublet0Case3", py::base<Case3>())
        .def(py::init<vector<Panel3*>>())
        .def(py::init<vector<Panel3*>, vector<PanelVector3*>>());
        
    py::class_<DirichletDoublet0Source0Case3>(m, "DirichletDoublet0Source0Case3", py::base<DirichletDoublet0Case3>())
        .def(py::init<vector<Panel3*>>())
        .def(py::init<vector<Panel3*>, vector<PanelVector3*>>());
        
    py::class_<LineSegment>(m, "LineSegment")
        .def(py::init<Vector3&, Vector3&, Vector3&>())
        .def_readonly("best_gamma", &LineSegment::best_gamma)
        .def_readonly("gamma", &LineSegment::gamma)
        .def_property_readonly("ca", &LineSegment::ca)
        .def_readonly("n", &LineSegment::n)
        .def_readonly("t", &LineSegment::t)
        .def_property_readonly("b", &LineSegment::b)
        .def_readonly("mids", &LineSegment::mid) 
        .def_readonly("v1", &LineSegment::v1)
        .def_readonly("v2", &LineSegment::v2)
        .def_readonly("lift_factor", &LineSegment::lift_factor);

    py::class_<LiftingLine>(m, "LiftingLine")
        .def(py::init<vector<Vector3*>>())
        .def_readwrite("v_inf", &LiftingLine::v_inf)
        .def_property_readonly("segments", &LiftingLine::get_segments)
        .def_property_readonly("points", &LiftingLine::get_points)
        .def("solve_for_best_gamma", &LiftingLine::solve_for_best_gamma);

    /************************-3D SINGULARITY ELEMENTS-************************/
    m.def("doublet_3_0_vsaero", &wrap_doublet_3_0_vsaero);
    m.def("doublet_3_0_sphere", &wrap_doublet_3_0_sphere);
    m.def("doublet_3_0_n0", &wrap_doublet_3_0_n0);
    m.def("doublet_src_3_0_vsaero", &wrap_doublet_src_3_0_vsaero);
    m.def("doublet_src_3_0_n0", &wrap_doublet_src_3_0_n0);
    m.def("doublet_3", &doublet_3);
    m.def("doublet_3_0_vsaero_v", &wrap_doublet_3_0_vsaero_v);
    m.def("doublet_src_3_0_vsaero_v", &doublet_src_3_0_vsaero_v);
    m.def("src_3_0_vsaero_v", &wrap_source_3_0_vsaero_v);
    m.def("vortex_3_0_v", &wrap_vortex_3_0_v);
    m.def("vortex_3_0_v", &wrap_vortex_3_0_edge_v);
    m.def("vortex_3_0_half_infinity_v", &vortex_3_0_half_infinity_v);
    m.def("doublet_3_v", &doublet_3_v);

    /************************-2D SINGULARITY ELEMENTS-************************/
    m.def("source_2", &source_2);
    m.def("doublet_2", &doublet_2);
    m.def("vortex_2", &vortex_2);
    m.def("source_2_0", &source_2_0);
    m.def("doublet_2_0", &doublet_2_0);
    m.def("doublet_2_1", &doublet_2_1);
    m.def("source_2_v", &source_2_v);
    m.def("doublet_2_v", &doublet_2_v);
    m.def("vortex_2_v", &vortex_2_v);
    m.def("source_2_0_v", &source_2_0_v);
    m.def("doublet_2_0_v", &doublet_2_0_v);
    m.def("doublet_2_1_v", &doublet_2_1_v);

    return m.ptr();
};

