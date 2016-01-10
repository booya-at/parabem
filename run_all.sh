#!/usr/bin/bash
# argument is the version number of the python interpreter
BASEDIR=$(dirname $0)
cd $BASEDIR/ppm_tests

cd tests

$1 ../test_it.py test_airfoil.py
$1 ../test_it.py test_case.py
$1 ../test_it.py test_doublet.py
$1 ../test_it.py test_doublet_src.py
$1 ../test_it.py test_linear_2d_doublet.py
$1 ../test_it.py test_mesh.py
$1 ../test_it.py test_parallel.py
$1 ../test_it.py test_parallel_2d.py

cd ../plots

$1 ../test_it.py 2d_elements.py
$1 ../test_it.py alpha_cL_2d.py
$1 ../test_it.py alpha_cL_3d.py
$1 ../test_it.py cluster_test.py
$1 ../test_it.py cylinder_2d.py
$1 ../test_it.py far_field_error.py
$1 ../test_it.py far_field_error_src.py
$1 ../test_it.py joukowsky_surface.py
$1 ../test_it.py lifting_line.py
$1 ../test_it.py panel_doublet.py
$1 ../test_it.py panel_doublet_far.py
$1 ../test_it.py panel_doublet_src.py
$1 ../test_it.py panel_src.py
$1 ../test_it.py panel_src_far.py
$1 ../test_it.py sphere_error.py
$1 ../test_it.py sphere_error_new.py
$1 ../test_it.py vortex.py

cd ../vtk

$1 ../test_it.py vtk_airfoil.py
$1 ../test_it.py vtk_airfoil_joukowsky.py
$1 ../test_it.py vtk_airfoil_linear_dirichlet.py
$1 ../test_it.py vtk_airfoil_neumann.py
$1 ../test_it.py vtk_airfoil_vandevooren.py
$1 ../test_it.py vtk_coordinates_not_aligned.py
$1 ../test_it.py vtk_cylinder_linear_dirichlet.py
$1 ../test_it.py vtk_element_vortex.py
$1 ../test_it.py vtk_matrix_plot.py
$1 ../test_it.py vtk_pan3d_obj.py
$1 ../test_it.py vtk_panel_doublet.py
$1 ../test_it.py vtk_panel_source.py
$1 ../test_it.py vtk_parallel_flow.py
$1 ../test_it.py vtk_symmetric.py
$1 ../test_it.py vtk_vortex_2.py

# cd ../openglider

# $1 ../test_it.py alpha_cL_cD.py
# $1 ../test_it.py shape_cl_cd.py
# $1 ../test_it.py lifting_line_opt.py
# $1 ../test_it.py vtk_glider.py
# $1 ../test_it.py wake_models.py