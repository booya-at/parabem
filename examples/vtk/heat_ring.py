from __future__ import division
import numpy as np
import paraBEM
from paraBEM.pan2d import doublet_2_0, doublet_2_0_v, source_2_0, source_2_0_v
from paraBEM.vtk_export import VtkWriter
from paraBEM.utils import check_path

R = 1

points = [
    [-2, -2], [2, -2], [2, 2], [-2, 2],
    [-1, -1], [1, -1], [1, 1], [-1, 1]
]

points = [paraBEM.PanelVector2(*point) for point in points]
panels = [
    [1, 0], [2, 1], [3, 2], [0, 3],
    [4, 5], [5, 6], [6, 7], [7, 4]
]
panels = [paraBEM.Panel2([points[i] for i in panel]) for panel in panels]

T_boundary = [0, 0, 0, 0, 1, 1, 1, 1]
mat = np.zeros([8, 8])
rhs = np.zeros([8])
for i, panel_i in enumerate(panels):
    rhs[i] = T_boundary[i]
    for j, panel_j in enumerate(panels):
        mat[i, j] = - source_2_0(panel_i.center, panel_j)
        mat[i, j] += doublet_2_0(panel_i.center, panel_j) * R

sol = np.linalg.solve(mat, rhs)


nx = 300
ny = 300
x_grid = np.linspace(-3, 3, nx)
y_grid = np.linspace(-3, 3, ny)
grid = [paraBEM.Vector2(x, y) for y in y_grid for x in x_grid]
t_list = []
for point in grid:
    t = 0
    for i, panel_i in enumerate(panels):
        t -= source_2_0(point, panel_i) * sol[i]
        t += doublet_2_0(point, panel_i) * sol[i] * R
    t_list.append(t)

q_list = []
for point in grid:
    q = paraBEM.Vector2(0, 0)
    for i, panel_i in enumerate(panels):
        q -= source_2_0_v(point, panel_i) * sol[i]
        q += doublet_2_0_v(point, panel_i) * sol[i] * R
    q_list.append(q)


writer = VtkWriter()
with open(check_path("results/heat_test.vtk"), "w") as _file:
    writer.structed_grid(_file, "element_2", [nx, ny, 1])
    writer.points(_file, grid)
    writer.data(_file, t_list, name="temperature", _type="SCALARS", data_type="POINT_DATA")
    writer.data(_file, q_list, name="q", _type="VECTORS", data_type="POINT_DATA")
