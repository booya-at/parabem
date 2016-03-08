# -*- coding: utf-8 -*-
import numpy as np

import ppm
from ppm.vtk_export import VtkWriter
from ppm import PanelVector2, Vector2, Panel2
from ppm.airfoil import Airfoil
from ppm.pan2d import NeumannSource0Case2 as Case
from ppm.utils import check_path

airfoil = Airfoil.joukowsky(m=-0.1 +0.1j)

airfoil.numpoints = 50
alpha = np.deg2rad(10)

# panelmethode
case = Case(airfoil.panels)
case.v_inf = Vector2(np.cos(alpha), np.sin(alpha))
case.run()
print(np.array(case.matrix.values))
nx = 300
ny = 300

space_x = np.linspace(-1, 2, nx)
space_y = np.linspace(-0.2, 0.2, ny)
vec = lambda x: ppm.Vector2(x[0], x[1])
vec3 = lambda x: [x[0], x[1], 0]
grid = [ppm.Vector2(x, y) for y in space_y for x in space_x]

velocity = list(map(vec3, map(case.off_body_velocity, grid)))
vel1 = [(i[0]**2 + i[1]**2)**(0.5) for i in velocity]
pot = list(map(case.off_body_potential, grid))

with open(check_path("results/neumann/field.vtk"), "w") as _file:
    writer = VtkWriter()
    writer.structed_grid(_file, "airfoil", [nx, ny, 1])
    writer.points(_file, grid)
    writer.data(_file, velocity, name="velocity", _type="VECTORS", data_type="POINT_DATA")
    writer.data(_file, pot, name="pot", _type="SCALARS", data_type="POINT_DATA")
    writer.data(_file, vel1, name="vel", _type="SCALARS", data_type="POINT_DATA")

with open(check_path("results/neumann/airfoil.vtk"), "w") as _file:
    writer = VtkWriter()
    writer.unstructed_grid(_file, "airfoil")
    writer.points(_file, [list(i.points[0]) + [0] for i in case.panels])
    writer.lines(_file, [range(len(case.panels))])
