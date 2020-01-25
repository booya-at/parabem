# -*- coding: utf-8 -*-
import numpy as np

import parabem
from parabem.airfoil import Airfoil
from parabem.utils import check_path
from parabem.vtk_export import VtkWriter
from parabem.pan2d import DirichletDoublet1Case2 as Case

# geometry
airfoil = Airfoil.trefftz_kutta(-0.1 + 0.0j, np.deg2rad(30))
airfoil.numpoints = 50
points = [parabem.PanelVector2(*i) for i in airfoil.coordinates[:-1]]
points[0].wake_vertex = True
panels = [parabem.Panel2([coord, points[i+1]]) for i, coord in enumerate(points[:-1])]
panels.append(parabem.Panel2([points[-1], points[0]]))

# panelmethode
case = Case(panels)
case.v_inf = parabem.Vector2(1, 0.3)
case.run()
# plt.plot(np.array(case.matrix.values).T)
# # plt.show()

nx = 200
ny = 200
space_x = np.linspace(-0.2, 1.5, nx)
space_y = np.linspace(-0.5, 0.5, ny)
grid = [parabem.Vector2(x, y) for y in space_y for x in space_x]

velocity = list(map(case.off_body_velocity, grid))
pot = list(map(case.off_body_potential, grid))

with open(check_path("results/airfoil_2d_linear/field.vtk"), "w") as _file:
    writer = VtkWriter()
    writer.structed_grid(_file, "airfoil", [nx, ny, 1])
    writer.points(_file, grid)
    writer.data(_file, velocity, name="velocity", _type="VECTORS", data_type="POINT_DATA")
    writer.data(_file, pot, name="pot", _type="SCALARS", data_type="POINT_DATA")

with open(check_path("results/airfoil_2d_linear/shape.vtk"), "w") as _file:
    writer = VtkWriter()
    writer.unstructed_grid(_file, "airfoil")
    writer.points(_file, [[i[0], i[1], 0]for i in airfoil.coordinates])
    writer.lines(_file, [range(len(airfoil.coordinates))])