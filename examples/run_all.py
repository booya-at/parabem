import subprocess as sub
import os

directory = os.path.dirname(__file__)
print(directory)


def run_example(sub_dir, file_name):
    print(f"\n\nRunning {sub_dir}/{file_name}\n")
    sub.run(["python", os.path.join(directory, sub_dir, file_name)],
            cwd = os.path.join(directory, sub_dir),
            check=False)

run_example("tests", "test_airfoil.py")
run_example("tests", "test_case.py")
run_example("tests", "test_doublet.py")
run_example("tests", "test_doublet_src.py")
run_example("tests", "test_linear_2d_doublet.py")
run_example("tests", "test_mesh.py")
run_example("tests", "test_parallel.py")
run_example("tests", "test_parallel_2d.py")

run_example("plots", "2d_elements.py")
run_example("plots", "alpha_cL_2d.py")
run_example("plots", "alpha_cL_3d.py")
run_example("plots", "cluster_test.py")
run_example("plots", "cylinder_2d.py")
run_example("plots", "far_field_error.py")
run_example("plots", "far_field_error_src.py")
run_example("plots", "joukowsky_surface.py")
run_example("plots", "lifting_line.py")
run_example("plots", "panel_doublet.py")
run_example("plots", "panel_doublet_far.py")
run_example("plots", "panel_doublet_src.py")
run_example("plots", "panel_src.py")
run_example("plots", "panel_src_far.py")
run_example("plots", "sphere_error.py")
run_example("plots", "sphere_error_new.py")
run_example("plots", "vortex.py")

run_example("vtk", "vtk_airfoil.py")
run_example("vtk", "vtk_airfoil_joukowsky.py")
run_example("vtk", "vtk_airfoil_linear_dirichlet.py")
run_example("vtk", "vtk_airfoil_neumann.py")
run_example("vtk", "vtk_airfoil_vandevooren.py")
run_example("vtk", "vtk_coordinates_not_aligned.py")
run_example("vtk", "vtk_cylinder_linear_dirichlet.py")
run_example("vtk", "vtk_element_vortex.py")
run_example("vtk", "vtk_pan3d_obj.py")
run_example("vtk", "vtk_panel_doublet.py")
run_example("vtk", "vtk_panel_source.py")
run_example("vtk", "vtk_parallel_flow.py")
run_example("vtk", "vtk_symmetric.py")
run_example("vtk", "vtk_vortex_2.py")

run_example("openglider", "alpha_cL_cD.py")
run_example("openglider", "shape_cl_cd.py")
run_example("openglider", "lifting_line_opt.py")
run_example("openglider", "vtk_glider.py")
run_example("openglider", "wake_models.py")
run_example("openglider", "crash.py")
