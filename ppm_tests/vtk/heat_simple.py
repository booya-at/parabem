from __future__ import division
import numpy as np
import ppm
from ppm.pan2d import doublet_2_0, doublet_2_0_v, source_2_0, source_2_0_v
from ppm.vtk_export import VtkWriter
from ppm.utils import check_path

p1 = ppm.PanelVector2(-1, -10)
p2 = ppm.PanelVector2(1, -10)
p3 = ppm.PanelVector2(-1, 10)
p4 = ppm.PanelVector2(1, 10)

pan1 = ppm.Panel2([p1, p3])
pan2 = ppm.Panel2([p4, p2])

mat = np.zeros([2, 2])
rhs = np.zeros([2])

#    |    |
# col|+  +|
#    |    |
# x -1    1
# T  0    1
# R  1    1
# l    1
#   panel1: temp-formulation
T1 = 0
T2 = 1
R1 = 0
R2 = 0
l = 1
mat[0, 0] = source_2_0(pan1.center, pan1)
mat[0, 1] = source_2_0(pan1.center, pan2)
mat[0, 0] -= doublet_2_0(pan1.center, pan1) * l * R1
mat[0, 1] -= doublet_2_0(pan1.center, pan2) * l * R2
rhs[0] += T1

#   panel3: temp-formulation
mat[1, 0] = source_2_0(pan2.center, pan1)
mat[1, 1] = source_2_0(pan2.center, pan2)
mat[1, 0] -= doublet_2_0(pan2.center, pan1) * l * R1
mat[1, 1] -= doublet_2_0(pan2.center, pan2) * l * R2
rhs[1] += T2

sol = np.linalg.solve(mat, rhs)
print(mat)
print(rhs)
print(sol)


nx = 300
ny = 300
x_grid = np.linspace(-3, 3, nx)
y_grid = np.linspace(-3, 3, ny)


grid = [ppm.Vector2(x, y) for y in y_grid for x in x_grid]
temp_list = []
for point in grid:
    t = 0
    t -= doublet_2_0(point, pan1) * sol[0] * R1 * l
    t -= doublet_2_0(point, pan2) * sol[1] * R2 * l
    t += source_2_0(point, pan1) * sol[0]
    t += source_2_0(point, pan2) * sol[1]
    temp_list.append(t)

q_list = []
for point in grid:
    q = ppm.Vector2(0, 0)
    q -= doublet_2_0_v(point, pan1) * sol[0] * R1 * l
    q -= doublet_2_0_v(point, pan2) * sol[1] * R2 * l
    q += source_2_0_v(point, pan1) * sol[0]
    q += source_2_0_v(point, pan2) * sol[1]
    q_list.append(q)

writer = VtkWriter()
with open(check_path("results/heat_test.vtk"), "w") as _file:
    writer.structed_grid(_file, "element_2", [nx, ny, 1])
    writer.points(_file, grid)
    writer.data(_file, temp_list, name="temperature", _type="SCALARS", data_type="POINT_DATA")
    writer.data(_file, q_list, name="q", _type="VECTORS", data_type="POINT_DATA")
