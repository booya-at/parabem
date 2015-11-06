import matplotlib.pyplot as plt
import numpy as np

from openglider.jsonify import load
from openglider.airfoil import Profile2D
from openglider.glider.in_out.export_3d import ppm_Panels

import ppm
from ppm.pan3d import DirichletDoublet0Source0Case3 as Case
from ppm.vtk_export import CaseToVTK
from ppm.utils import vinf_deg_range3

n_x = 80

with open("glider/referenz_schirm_berg.json") as _file:
    glider_2d = load(_file)["data"]
    glider_2d.set_const_cell_dist()
    glider = glider_2d.get_glider_3d()

_, panels, trailing_edge = ppm_Panels(
    glider,
    midribs=7,
    profile_numpoints=n_x,
    distribution=Profile2D.cos_2_distribution,
    num_average=0,
    symmetric=False)

v_inf = [8, 0, 1]

case = Case(panels, trailing_edge)
case.mom_ref_point = ppm.Vector3(1.25, 0, 0)
case.A_ref = glider_2d.flat_area
case.farfield = 5
case.drag_calc = "on_body"

case.v_inf = ppm.Vector3(v_inf)
case.create_wake(length=10000, count=4)    # length, count
polars = case.polars(vinf_deg_range3(case.v_inf, -5, 15, 50))


vtk_writer = CaseToVTK(case, "results/vtk_glider_case")
vtk_writer.write_panels(data_type="cell")
vtk_writer.write_wake_panels()
vtk_writer.write_body_stream(panels, 100)

arr = []
for i in range(len(case.panels) / n_x):
    arr.append([[i.center.x, i.cp] for i in panels[i * n_x: (i + 1) * n_x]])
for i in arr:
    plt.plot(*zip(*i))
plt.show()

arr = []
for i in range(len(case.panels) / n_x):
    arr.append([[i.center.y, i.cp] for i in panels[i * n_x: (i + 1) * n_x]])
for i in zip(*arr):
    plt.plot(*zip(*i))
plt.show()

p = []
alpha = []
for i in polars.values:
    alpha.append(i.alpha)
    p.append((i.cL, i.cD, i.cP))
plt.plot(np.array(p), alpha)
plt.grid()
plt.show()