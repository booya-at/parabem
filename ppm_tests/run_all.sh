#!/usr/bin/bash
# argument is the version number of the python interpreter
cd tests

python$1 ../test_it.py test_airfoil.py
python$1 ../test_it.py test_case.py
python$1 ../test_it.py test_doublet.py
python$1 ../test_it.py test_doublet_src.py
python$1 ../test_it.py test_linear_2d_doublet.py
python$1 ../test_it.py test_mesh.py
python$1 ../test_it.py test_parallel.py
python$1 ../test_it.py test_parallel_2d.py

cd ../plots

python$1 ../test_it.py 2d_elements.py
python$1 ../test_it.py alpha_cL_2d.py
python$1 ../test_it.py alpha_cL_3d.py
python$1 ../test_it.py cluster_test.py
python$1 ../test_it.py cylinder_2d.py
python$1 ../test_it.py far_field_error.py
python$1 ../test_it.py far_field_error_src.py
python$1 ../test_it.py joukowsky_surface.py
python$1 ../test_it.py lifting_line.py
python$1 ../test_it.py panel_doublet.py
python$1 ../test_it.py panel_doublet_far.py
python$1 ../test_it.py panel_doublet_src.py
python$1 ../test_it.py panel_src.py
python$1 ../test_it.py panel_src_far.py
python$1 ../test_it.py sphere_error.py
python$1 ../test_it.py vortex.py

cd ../openglider

python$1 ../test_it.py alpha_cL_cD.py
python$1 ../test_it.py shape_cl_cd.py
python$1 ../test_it.py lifting_line_opt.py
python$1 ../test_it.py vtk_glider.py
python$1 ../test_it.py wake_models.py

cd ../vtk

python$1 ../test_it.py vtk_airfoil.py
python$1 ../test_it.py vtk_airfoil_joukowsky.py
python$1 ../test_it.py vtk_airfoil_linear_dirichlet.py
python$1 ../test_it.py vtk_airfoil_neumann.py
python$1 ../test_it.py vtk_airfoil_vandevooren.py
python$1 ../test_it.py vtk_coordinates_not_aligned.py
python$1 ../test_it.py vtk_cylinder_linear_dirichlet.py
python$1 ../test_it.py vtk_element_vortex.py
python$1 ../test_it.py vtk_matrix_plot.py
python$1 ../test_it.py vtk_pan3d_obj.py
python$1 ../test_it.py vtk_panel_doublet.py
python$1 ../test_it.py vtk_panel_source.py
python$1 ../test_it.py vtk_parallel_flow.py
python$1 ../test_it.py vtk_symmetric.py
python$1 ../test_it.py vtk_vortex_2.py
