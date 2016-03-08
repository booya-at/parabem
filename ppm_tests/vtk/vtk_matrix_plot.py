import numpy

import ppm
from ppm.mesh import mesh_object
from ppm import pan3d
from ppm.vtk_export import VtkWriter
from ppm.utils import check_path


mesh = mesh_object.from_OBJ("../mesh/sphere_tri.obj")

case = pan3d.DirichletDoublet0Case3(mesh.panels, mesh.trailing_edges)
case.drag_calc = "on_body"
case.v_inf = ppm.Vector3(0, 1, 0)
case.create_wake(3, 3)
case.run()

n = case.mat_size

space = numpy.linspace(-2, 2, n).tolist()
index_pos = [[x, y, 0] for y in space for x in space]

mat = case.matrix.values
mat_flat = []
for row in mat:
    for value in row:
        if abs(value) < 1:
            mat_flat.append(value)
        else:
            mat_flat.append(0.)


writer = VtkWriter()
with open(check_path("results/matrix.vtk"), "w") as _file:
    writer.structed_grid(_file, "matrix", [n, n, 1])
    writer.points(_file, index_pos)
    writer.data(_file, mat_flat, name="matrix", _type="SCALARS", data_type="POINT_DATA")
