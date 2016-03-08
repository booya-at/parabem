#include <boost/python.hpp>
#include <boost/python/overloads.hpp>
#include <boost/python/stl_iterator.hpp>
#include <vector>

#include "vector.h"
#include "panel3.h"
#include "panel2.h"
#include "case3.h"
#include "case2.h"
#include "ppm_wrap.h"
#include "element_influence.h"
#include "lifting_line.h"

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(Case3_relax_wake_overloads, Case3::relax_wake, 0, 2);
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(Case3_create_wake_overloads, Case3::create_wake, 0, 3);
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(Case2_create_wake_overloads, Case2::create_wake, 0, 3);

using namespace boost::python;
BOOST_PYTHON_MODULE(_ppm){
    docstring_options local_docstring_options(true, false, false);
    
    class_<vector<Vector3>>("VectorList", no_init)
        .add_property("values", &to_python_list<Vector3>);
        
    class_<vector<Edge>>("EdgeList", no_init)
        .add_property("values", &to_python_list<Edge>);
        
    class_<vector<PanelVector3*>>("PanelVectorList", no_init)
        .add_property("values", &to_python_list<PanelVector3*>);
        
    class_<vector<PanelVector2*>>("PanelVector2List", no_init)
        .add_property("values", &to_python_list<PanelVector2*>);
    
    class_<Eigen::VectorXd>("VectorXd", no_init)
        .add_property("values", &py_eigen::vec_to_list);
    
    class_<Eigen::MatrixXd>("MatrixXd", no_init)
        .add_property("values", &py_eigen::mat_to_list);
        
    class_<vector<AeroCoef3>>("AeroCoef3List", no_init)
        .add_property("values", &to_python_list<AeroCoef3>);
        
    class_<AeroCoef3>("AeroCoef3", no_init)
        .add_property("alpha", &AeroCoef3::alpha)
        .add_property("v_inf", &AeroCoef3::v_inf)
        .add_property("force", &AeroCoef3::force)
        .add_property("cop", &AeroCoef3::cop)
        .add_property("cL", &AeroCoef3::cL)
        .add_property("cD", &AeroCoef3::cD)
        .add_property("cS", &AeroCoef3::cS)
        .add_property("cR", &AeroCoef3::cR)
        .add_property("cP", &AeroCoef3::cP)
        .add_property("cY", &AeroCoef3::cY);
        
    class_<Polar3, bases<vector<AeroCoef3>>>("Polar3")
        .add_property("labels", &py_polar::labels<Polar3>)
        .add_property("as_matrix", &py_polar::as_matrix<Polar3>);

    class_<Vector2>("Vector2", "Vector2 is the ppm version of a Eigen::Vector2d\n\
		    constructor: Vector2(x, y)")
        .def("__init__", make_constructor(&py_vector_2::init_from_double<Vector2>))
        .def("__init__", make_constructor(&py_vector_2::init_from_list<Vector2>))
        .def("__init__", make_constructor(&py_vector_2::init_from_Vector3))
        .def("__init__", make_constructor(&py_vector_2::init_from_Vector2))
        .def("__getitem__", &py_vector_2::getitem)
        .def("__setitem__", &py_vector_2::setitem)
        .add_property("x", &py_vector_2::get_x, &py_vector_2::set_x)
        .add_property("y", &py_vector_2::get_y, &py_vector_2::set_y)
        .def("norm", &Vector2::norm)
        .add_property("normal", &normal2)
        .def("dot", &py_vector_2::dot)
        .def("__abs__", &py_vector_2::abs)
        .def("__len__", &Vector2::size)
        .def("__repr__", &py_vector_2::repr)
        .def("__add__", &py_vector_2::add)
        .def("__sub__", &py_vector_2::sub)
        .def("__neg__", &py_vector_2::neg)
        .def("__pos__", &py_vector_2::pos)
        .def("__mul__", &py_vector_2::mul)
        .def("__rmul__", &py_vector_2::mul)
        .def("__truediv__", &py_vector_2::div)
        ;
        
    class_<PanelVector2, bases<Vector2>>("PanelVector2")
        .def("__init__", make_constructor(&py_vector_2::init_from_double<PanelVector2>))
        .def("__init__", make_constructor(&py_vector_2::init_from_list<PanelVector2>))
        .def_readwrite("wake_vertex", &PanelVector2::wake_vertex, "if true this point shades a wake")
        .def_readwrite("potential", &PanelVector2::potential, "potential at this point (not all methodes set this value)")
        .def_readwrite("velocity", &PanelVector2::velocity, "velocity at this point (not all methodes set this value)")
        .def_readonly("nr", &PanelVector2::nr, "nr of this point");
        
    class_<Panel2>("Panel2", "Panel2 represents a straigth line which is used to approximate 2d geometry\n\
		    constructor: Panel2([p1, p2])\n\n\
		    p1, p2 -> PanelVector2 which have to be stored somewhere in python (list, variable,...)")
        .def("__init__", make_constructor(&py_panel::init_Panel<Panel2, PanelVector2>))
        .add_property("l", &Panel2::l, "length of panel")
        .add_property("t", &Panel2::t, "panel direction")
        .add_property("n", &Panel2::n, "normal direction")
        .add_property("center", &Panel2::center, "panel center")
        .add_property("points", &py_panel::panel_vectors<Panel2, PanelVector2>)
        .def_readwrite("potential", &Panel2::potential, "value is set from case")
        .def_readwrite("mue", &Panel2::mue, "value is set from case")
        .def_readwrite("sigma", &Panel2::sigma, "value is set from case")
        .add_property("cp", &Panel2::cp, "value is set from case")
        .add_property("velocity", &Panel2::velocity, "value is set from case");

    class_<Case2>("Case2", no_init)
        .add_property("panels", &py_case::case_panels<Case2, Panel2>, "panels of case")
        .add_property("points", &Case2::get_all_points, "all points")
        .def_readwrite("v_inf", &Case2::v_inf, "direction of parallell flow")
        .def_readwrite("A_ref", &Case2::A_ref, "reference Area")
        .def_readwrite("mom_ref_point", &Case2::mom_ref_point, "point to sum moment")
        .add_property("matrix", &Case2::matrix, "resulting matrix")
        .add_property("result", &Case2::result, "solution of the linear sysem")
        .add_property("rhs", &Case2::rhs, "the right hand side of the equation")
        .add_property("mat_size", &Case2::mat_size, "size of the matrix")
        .add_property("force", &Case2::force, "sum of all panel-forces")
        .add_property("center_of_pressure", &Case2::center_of_pressure, "center of pressure")
        .add_property("cl", &Case2::ca, "force normal to vinf")
        .add_property("cm", &Case2::moment, "resulting moment")
        .def("run", &Case2::run, "run the case")
        .def("off_body_velocity", &Case2::off_body_velocity, "velocity field")
        .def("off_body_potential", &Case2::off_body_potential, "potetnial field");
        
    class_<DirichletDoublet0Case2, bases<Case2>>("DirichletDoublet0Case2", no_init)
        .def("__init__", make_constructor(&py_case::init_from_pans_2<DirichletDoublet0Case2>));
        
    class_<NeumannDoublet0Case2, bases<Case2>>("NeumannDoublet0Case2", no_init)
        .def("__init__", make_constructor(&py_case::init_from_pans_2<NeumannDoublet0Case2>));
        
    class_<NeumannSource0Case2, bases<Case2>>("NeumannSource0Case2", no_init)
        .def("__init__", make_constructor(&py_case::init_from_pans_2<NeumannSource0Case2>));
        
    class_<DirichletDoublet0Source0Case2, bases<DirichletDoublet0Case2>>("DirichletDoublet0Source0Case2", no_init)
        .def("__init__", make_constructor(&py_case::init_from_pans_2<DirichletDoublet0Source0Case2>));
    
    class_<DirichletDoublet1Case2, bases<Case2>>("DirichletDoublet1Case2", no_init)
        .def("__init__", make_constructor(&py_case::init_from_pans_2<DirichletDoublet1Case2>));
    
    class_<Vector3>("Vector3")
        .def("__init__", make_constructor(&py_vector_3::init_from_double<Vector3>))
        .def("__init__", make_constructor(&py_vector_3::init_from_list<Vector3>))
        .def("__init__", make_constructor(&py_vector_3::init_from_Vector2))
        .def("__init__", make_constructor(&py_vector_3::init_from_Vector3))
        .def("__getitem__", &py_vector_3::getitem)
        .def("__setitem__", &py_vector_3::setitem)
        .add_property("x", &py_vector_3::get_x, &py_vector_3::set_x)
        .add_property("y", &py_vector_3::get_y, &py_vector_3::set_y)
        .add_property("z", &py_vector_3::get_z, &py_vector_3::set_z)
        .def("norm", &Vector3::norm)
        .def("__abs__", &py_vector_3::abs)
        .def("__len__", &Vector3::size)
        .def("__repr__", &py_vector_3::repr)
        .def("__add__", &py_vector_3::add)
        .def("__sub__", &py_vector_3::sub)
        .def("__mul__", &py_vector_3::mul)
        .def("__rmul__", &py_vector_3::mul)
        .def("__truediv__", &py_vector_3::div)
        .def("__neg__", &py_vector_3::neg)
        .def("__pos__", &py_vector_3::pos)
        .def("cross", &py_vector_3::cross)
        .def("dot", &py_vector_3::dot);
        
    class_<PanelVector3, bases<Vector3>>("PanelVector3", 
		"Panel3 represents a straigth line which is used to approximate 2d geometry\n\
		constructor: Panel2([p1, p2, ...])\n\n\
		p1, p2... -> PanelVector2 which have to be stored somewhere in python (list, variable,...)")
        .def("__init__", make_constructor(&py_vector_3::init_from_double<PanelVector3>))
        .def("__init__", make_constructor(&py_vector_3::init_from_list<PanelVector3>))
        .add_property("nr", &PanelVector3::nr)
        .def_readonly("potential", &PanelVector3::potential)
        .def_readonly("cp", &PanelVector3::cp)
        .def_readonly("vorticity", &PanelVector3::vorticity)
        .def_readonly("velocity", &PanelVector3::velocity);
    
    class_<Panel3>("Panel3", no_init)
        .def("__init__", make_constructor(&py_panel::init_Panel<Panel3, PanelVector3>))
        .def("append_point", &Panel3::append_point)
        .add_property("center", &Panel3::center)
        .add_property("m", &Panel3::m)
        .add_property("n", &Panel3::n)
        .add_property("l", &Panel3::l)
        .add_property("area", &Panel3::area)
        .add_property("potential", &Panel3::get_potential, &Panel3::set_potential)
        .add_property("sigma", &Panel3::get_sigma, &Panel3::set_sigma)
        .add_property("mue", &Panel3::get_mue)
        .add_property("velocity", &Panel3::get_velocity)
        .add_property("cp", &Panel3::get_cp)
        .add_property("points", &py_panel::panel_vectors<Panel3, PanelVector3>)
        .def("set_symmetric", &Panel3::set_symmetric);
        
    class_<Edge>("Edge")
        .def_readonly<>("v1", &Edge::v1)
        .def_readonly<>("v2", &Edge::v2)
        .def_readonly<>("vorticity", &Edge::vorticity)
        .add_property("center", &Edge::center);
        
    class_<WakePanel3, bases<Panel3>>("WakePanel3");
//         .add_property<>("upper_operating_Panel3", &WakePanel3::get_upper_operating_panel)
//         .add_property<>("lower_operating_Panel3", &WakePanel3::get_lower_operating_panel);

    class_<Case3>("Case3", no_init)
        .add_property("panels", &py_case::case_panels<Case3, Panel3>)
        .add_property("trailing_edge", &Case3::get_trailing_edge)
        .def_readwrite("v_inf", &Case3::v_inf)
        .def_readwrite("farfield", &Case3::farfield)
        .add_property("mat_size", &Case3::mat_size)
        .add_property("matrix", &Case3::matrix)
        .add_property("result", &Case3::result)
        .add_property("rhs", &Case3::rhs)
        .def_readwrite("drag_calc", &Case3::drag_calc)
        .def_readwrite("trefftz_cut_pos", &Case3::trefftz_cut_pos)
        .def("trefftz_cut", &Case3::trefftz_cut)
        .def("create_wake", &Case3::create_wake, Case3_create_wake_overloads(args("length", "count", "direction")))
        .def("relax_wake", &Case3::relax_wake,  Case3_relax_wake_overloads(args("iterations", "smoothing")))
        .def_readwrite("symmetric_plane_n", &Case3::symmetric_plane_n)
        .def_readwrite("symmetric_plane_p", &Case3::symmetric_plane_p)
        .add_property("vertices", &Case3::vertices)
        .def_readwrite("A_ref", &Case3::A_ref)
        .def_readwrite("mom_ref_point", &Case3::mom_ref_point)
        .def_readwrite("lift_ref", &Case3::lift_ref)
        .add_property("force", &Case3::force)
        .add_property("center_of_pressure", &Case3::center_of_pressure)
        .add_property("cL", &Case3::cL, "scalar lift cooeficient: cL = L * 2 / (rho * u**2 * A_ref)")
        .add_property("cD", &Case3::cD, "scalar drag cooeficient: cD = D * 2 / (rho * u**2 * A_ref)")        
        .add_property("cS", &Case3::cS, "scalar side lift cooeficient: cS = S * 2 / (rho * u**2 * A_ref)")
        .add_property("cM", &Case3::cM, "Vector3 moment cooeficients: cM = M * 2 / (rho * u**2 * A_ref * 1)")
        .add_property("wake_panels", &py_case::case_wake_panels< Case3, WakePanel3>)
        .def("get_surface_area", &Case3::get_surface_area)
        .def("get_volume", &Case3::get_volume)
        .def("get_projected_area", &Case3::get_projected_area)
        .def("off_body_velocity", &Case3::off_body_velocity)
        .def("off_body_potential",&Case3::off_body_potential)
        .def("flow_path", &Case3::flow_path)
        .def("body_flow_path", &Case3::body_flow_path)
        .def("run", &Case3::run)
        .def("polars", &py_case::polar3);
        
    class_<DirichletDoublet0Case3, bases<Case3>>("DirichletDoublet0Case3", no_init)
        .def("__init__", make_constructor(&py_case::init_from_pans_and_te<DirichletDoublet0Case3>))
        .def("__init__", make_constructor(&py_case::init_from_pans<DirichletDoublet0Case3>));
        
    class_<DirichletDoublet0Source0Case3, bases<DirichletDoublet0Case3>>("DirichletDoublet0Source0Case3", no_init)
        .def("__init__", make_constructor(&py_case::init_from_pans_and_te<DirichletDoublet0Source0Case3>))
        .def("__init__", make_constructor(&py_case::init_from_pans<DirichletDoublet0Source0Case3>));

    class_<LineSegment>("LineSegment", init<Vector3&, Vector3&, Vector3&>())
        .add_property("best_gamma", &LineSegment::best_gamma)
        .def_readonly("gamma", &LineSegment::gamma)
        .def_readonly("ca", &LineSegment::ca)
        .def_readonly("n", &LineSegment::n)
        .def_readonly("t", &LineSegment::t)
        .def_readonly("b", &LineSegment::b)
        .def_readonly("mids", &LineSegment::mid) 
        .add_property("v1", &LineSegment::py_v1)
        .add_property("v2", &LineSegment::py_v2)
        .def_readonly("lift_factor", &LineSegment::lift_factor);
       
    class_<vector<LineSegment*>>("VectorList")
        .add_property<>("values", &to_python_list<LineSegment*>);

    class_<LiftingLine>("LiftingLine")
        .def("append_point", &LiftingLine::append_point)
        .def("initialize", &LiftingLine::initialize)
        .def_readonly("segments", &LiftingLine::segments)
        .def("best_gamma", &LiftingLine::best_gamma);
    
    /************************-3D SINGULARITY ELEMENTS-************************/
    def("doublet_3_0_vsaero", &py_element_influence::wrap_doublet_3_0_vsaero);
    def("doublet_3_0_sphere", &py_element_influence::wrap_doublet_3_0_sphere);
    def("doublet_3_0_n0", &py_element_influence::wrap_doublet_3_0_n0);
    def("doublet_src_3_0_vsaero", &py_element_influence::wrap_doublet_src_3_0_vsaero);
    def("doublet_src_3_0_n0", &py_element_influence::wrap_doublet_src_3_0_n0);
    def("doublet_3", &doublet_3);
    def("doublet_3_0_vsaero_v", &py_element_influence::wrap_doublet_3_0_vsaero_v);
    def("doublet_src_3_0_vsaero_v", &doublet_src_3_0_vsaero_v);
    def("src_3_0_vsaero_v", &py_element_influence::wrap_source_3_0_vsaero_v);
    def("vortex_3_0_v", &py_element_influence::wrap_vortex_3_0_v);
    def("vortex_3_0_v", &py_element_influence::wrap_vortex_3_0_edge_v);
    def("vortex_3_0_half_infinity_v", &vortex_3_0_half_infinity_v);
    def("doublet_3_v", &doublet_3_v);

    /************************-2D SINGULARITY ELEMENTS-************************/
    def("source_2", &source_2);
    def("doublet_2", &doublet_2);
    def("vortex_2", &vortex_2);
    def("source_2_0", &source_2_0);
    def("doublet_2_0", &doublet_2_0);
    def("doublet_2_1", &doublet_2_1);
    def("source_2_v", &source_2_v);
    def("doublet_2_v", &doublet_2_v);
    def("vortex_2_v", &vortex_2_v);
    def("source_2_0_v", &source_2_0_v);
    def("doublet_2_0_v", &doublet_2_0_v);
    def("doublet_2_1_v", &doublet_2_1_v);
}