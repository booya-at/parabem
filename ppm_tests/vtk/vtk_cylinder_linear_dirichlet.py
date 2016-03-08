import numpy as np
import ppm
from ppm.pan2d import DirichletDoublet1Case2 as Case
from ppm.utils import check_path
from ppm.vtk_export import VtkWriter

# geometry
numpoints = 100
phi = np.linspace(0, 2 * np.pi, numpoints + 1)
x = np.cos(phi)[:-1]
y = np.sin(phi)[:-1]
xy = np.transpose(np.array([x, y]))

# mapping the geometry
vector = [ppm.PanelVector2(*i) for i in xy]
vector += [vector[0]]               # important for calculating the gradients
panels = [ppm.Panel2([vec, vector[i+1]]) for i, vec in enumerate(vector[:-1])]
vector[0].wake_vertex = True
# setting up the case
case = Case(panels)
case.v_inf = ppm.Vector2(1, 0.5)
case.run()

nx = 100
ny = 100
space_x = np.linspace(-2, 2, nx)
space_y = np.linspace(-2, 2, ny)
grid = [ppm.Vector2(x, y) for y in space_y for x in space_x]

velocity = list(map(case.off_body_velocity, grid))
pot = list(map(case.off_body_potential, grid))

with open(check_path("results/cylinder_2d_linear/field.vtk"), "w") as _file:
    writer = VtkWriter()
    writer.structed_grid(_file, "airfoil", [nx, ny, 1])
    writer.points(_file, grid)
    writer.data(_file, velocity, name="velocity", _type="VECTORS", data_type="POINT_DATA")
    writer.data(_file, pot, name="pot", _type="SCALARS", data_type="POINT_DATA")

with open(check_path("results/cylinder_2d_linear/shape.vtk"), "w") as _file:
    writer = VtkWriter()
    writer.unstructed_grid(_file, "airfoil")
    writer.points(_file, [[i[0], i[1], 0]for i in vector])
    writer.lines(_file, [range(len(vector))])