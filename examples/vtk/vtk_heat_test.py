from __future__ import division
import numpy as np
import parabem
from parabem.pan2d import doublet_2_0, doublet_2_0_v, source_2_0, source_2_0_v
from parabem.vtk_export import VtkWriter
from parabem.utils import check_path

p1 = parabem.PanelVector2(-1, -10)
p2 = parabem.PanelVector2(0, -10)
p3 = parabem.PanelVector2(1, -10)
p4 = parabem.PanelVector2(-1, 10)
p5 = parabem.PanelVector2(0, 10)
p6 = parabem.PanelVector2(1, 10)

pan1 = parabem.Panel2([p4, p1])
pan2 = parabem.Panel2([p2, p5])
pan3 = parabem.Panel2([p3, p6])


mat = np.zeros([3, 3])
rhs = np.zeros([3])
panels = [pan1, pan2, pan3]

#    |   |   |
# col|+ +|  +|
#    |   |   |
# x -1   0   1
# T  0   ?   1
# l    1   2
#   panel1: temp-formulation
T1 = -10
T2 = 10
l1 = 1
l2 = 2
mat[0, 0] = source_2_0(panels[0].center, panels[0])
mat[0, 1] = source_2_0(panels[0].center, panels[1]) * (1 - l2 / l1)
mat[0, 2] = source_2_0(panels[0].center, panels[2])
rhs[0] += doublet_2_0(panels[0].center, panels[0]) * T1
rhs[0] += doublet_2_0(panels[0].center, panels[2]) * T2
rhs[0] += T1

#   panel2: velocity formulation
mat[1, 0] = source_2_0_v(panels[1].center, panels[0]).dot(panels[1].n)
mat[1, 1] = source_2_0_v(panels[1].center, panels[1]).dot(panels[1].n) * (1 - l2 / l1) - 1
mat[1, 2] = source_2_0_v(panels[1].center, panels[2]).dot(panels[1].n)
rhs[1] += doublet_2_0_v(panels[1].center, panels[0]).dot(panels[1].n) * T1
rhs[1] += doublet_2_0_v(panels[1].center, panels[2]).dot(panels[1].n) * T2

#   panel3: temp-formulation
mat[2, 0] = source_2_0(panels[2].center, panels[0])
mat[2, 1] = source_2_0(panels[2].center, panels[1]) * (1 - l2 / l1)
mat[2, 2] = source_2_0(panels[2].center, panels[2])
rhs[2] += doublet_2_0(panels[2].center, panels[0]) * T1
rhs[2] += doublet_2_0(panels[2].center, panels[2]) * T2
rhs[2] += T2

sol = np.linalg.solve(mat, rhs)
print(mat)
print(rhs)
print(sol)


nx = 300
ny = 300
x_grid = np.linspace(-3, 3, nx)
y_grid = np.linspace(-3, 3, ny)


grid = [parabem.Vector2(x, y) for y in y_grid for x in x_grid]
t_list = []
for point in grid:
    t = 0
    t -= doublet_2_0(point, pan1) * T1
    t -= doublet_2_0(point, pan3) * T2
    t += source_2_0(point, pan1) * sol[0]
    t += source_2_0(point, pan2) * sol[1] * (1 - l2 / l1)
    t += source_2_0(point, pan3) * sol[2]

    t_list.append(t)

q_list = []
for point in grid:
    q = parabem.Vector2(0, 0)
    q -= doublet_2_0_v(point, pan1) * T1
    q -= doublet_2_0_v(point, pan3) * T2
    q += source_2_0_v(point, pan1) * sol[0]
    q += source_2_0_v(point, pan2) * sol[1] * (1 - l2 / l1)
    q += source_2_0_v(point, pan3) * sol[2]
    if point.x > -1 and point.x < 0:
        q *= l1
    if point.x > 0 and point.x < 1:
        q *= l2
    q_list.append(q)


writer = VtkWriter()
with open(check_path("results/heat_test.vtk"), "w") as _file:
    writer.structed_grid(_file, "element_2", [nx, ny, 1])
    writer.points(_file, grid)
    writer.data(_file, t_list, name="temperature", _type="SCALARS", data_type="POINT_DATA")
    writer.data(_file, q_list, name="q", _type="VECTORS", data_type="POINT_DATA")
