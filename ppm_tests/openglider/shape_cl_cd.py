from copy import deepcopy
import numpy as np
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt

from openglider.glider import glider_2d
from openglider.jsonify import dump, load
from openglider.glider.in_out.export_3d import ppm_Panels
from openglider.airfoil import Profile2D

import ppm
from ppm.pan3d import DirichletDoublet0Source0Case3 as Case
from ppm.vtk_export import CaseToVTK
from ppm.utils import check_path, vinf_deg_range3

count = 0
#   load the glider
with open("glider/referenz_schirm_berg.json") as _file:
    glider2d = load(_file)["data"]
    area = glider2d.flat_area
    aspect_ratio = glider2d.aspect_ratio


def glider_set_controlpoint(glider2d, x):
    c = deepcopy(glider2d.front.controlpoints)
    c[0, 0] = x
    glider2d.front.controlpoints = list(c)
    c = deepcopy(glider2d.back.controlpoints)
    c[0, 0] = x
    glider2d.back.controlpoints = list(c)
    glider2d.set_aspect_ratio(aspect_ratio, "span")
    # glider2d.set_flat_area(area, "aspect_ratio")
    glider2d.set_const_cell_dist()

def shape_plot(glider2d):
    shape = glider2d.shape()
    line = list(shape.front) + list(shape.back[::-1]) + [shape.front[0]]
    x, y = zip(*line)
    middle = sum(y) / len(y)
    y = [i - middle for i in y]
    return x, y

def min_func(x):
    global count
    count += 1
    _glider2d = deepcopy(glider2d)
    glider_set_controlpoint(_glider2d, x)
    glider3d = _glider2d.get_glider_3d()
    panels = ppm_Panels(glider3d,
                        midribs=0,
                        profile_numpoints=40,
                        symmetric=True,
                        distribution=Profile2D.nose_cos_distribution(0.2),
                        num_average=0)
    case = Case(panels[1], panels[2])
    case.farfield = 5
    case.drag_calc = "trefftz"
    case.A_ref = 23

    case.create_wake(length=10000, count=5)
    case.v_inf = ppm.Vector3(8, 0, 1)
    case.trefftz_cut_pos = case.v_inf * 50
    alpha = vinf_deg_range3(case.v_inf, 5, 9, 20)
    polars = case.polars(alpha)
    vtk_writer = CaseToVTK(case, "results/vtk_opt", suffix=str(count))
    vtk_writer.write_panels(data_type="point")
    cL = []
    cD = []
    for i in polars.values:
        cL.append(i.cL)
        cD.append(i.cD)
    return np.interp(0.6, cL, cD)


y_positions = np.linspace(-6, 10, 30)
cD = [min_func(i) for i in y_positions]
# plt.show()
plt.figure(figsize=(10,4))
plt.plot(y_positions, cD, marker="x")
plt.xlabel('y-Position der Kontrollpunkte [m]', fontsize=15)
plt.ylabel('induzierter Widerstandsbeiwert $c_{Wi}$', fontsize=15)
plt.grid(True)
plt.savefig(check_path('results/vtk_opt/induced_drag.png'),  bbox_inches='tight')
# # plt.show()
plt.close()


best_ind = list(cD).index(min(cD))
best_y = y_positions[best_ind]

glider_set_controlpoint(glider2d, best_y)
with open(check_path("results/glider/optimized.json"), "w") as _file:
    dump(glider2d, _file)

plt.figure(figsize=(10, 5))
plt.axes().set_aspect("equal", "datalim")
plt.grid(True)
plt.plot(*shape_plot(glider2d), color="0", label="optimal", linewidth=2)

glider_set_controlpoint(glider2d, y_positions[0])
with open(check_path("results/glider/min.json"), "w") as _file:
    dump(glider2d, _file)
plt.plot(*shape_plot(glider2d), color="r", label="untere Grenze")


glider_set_controlpoint(glider2d, y_positions[-1])
with open(check_path("results/glider/max.json"), "w") as _file:
    dump(glider2d, _file)

plt.plot(*shape_plot(glider2d), color="b", label="obere Grenze")
plt.legend()
plt.ylim((-2, 3))
plt.savefig(check_path('results/vtk_opt/shapes.png'), bbox_inches='tight')
