import ppm
from ppm import pan3d
from ppm.vtk_export import CaseToVTK
from ppm.mesh import mesh_object
import numpy as np

mesh = mesh_object.from_OBJ("../mesh/sphere_low_tri.obj")

case = pan3d.DirichletDoublet0Case3(mesh.panels)
case.v_inf = ppm.Vector3(1, 0, 0.)
case.farfield = 100
case.run()

lin = np.linspace(-0.5, 0.5, 5)
grid = [[-2, k, j] for j in lin for k in lin]

vtk_writer = CaseToVTK(case, "results/vtk_test_case")
vtk_writer.write_panels(data_type="point")
# vtk_writer.write_wake_panels()
vtk_writer.write_body_stream(mesh.panels, 50)
vtk_writer.write_field([-2, 2, 100], [-2, 2, 100], [-1, 1, 3])
