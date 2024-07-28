import os
from copy import deepcopy
import numpy as np
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt

# from openglider.glider.parametric import ParametricGlider
from openglider.jsonify import dump, load
from openglider.glider.in_out.export_3d import parabem_Panels
from openglider.utils.distribution import Distribution

import parabem
from parabem.pan3d import DirichletDoublet0Source0Case3 as Case
from parabem.vtk_export import CaseToVTK
from parabem.utils import check_path, v_inf_deg_range3

directory = os.path.dirname(__file__)

count = 0
with open(os.path.join(directory, "glider", "referenz_schirm_berg.json"), "r") as _file:
    glider2d = load(_file)["data"]
    area = glider2d.shape.area
    aspect_ratio = glider2d.shape.aspect_ratio

def results_path(*args):
    return check_path(os.path.join(directory, "..", "results", *args))

def glider_set_controlpoint(glider2d, x):
    c = deepcopy(glider2d.shape.front_curve.controlpoints)
    c[0, 0] = x
    glider2d.shape.front_curve.controlpoints = list(c)
    c = deepcopy(glider2d.shape.back_curve.controlpoints)
    c[0, 0] = x
    glider2d.shape.back_curve.controlpoints = list(c)
    glider2d.shape.set_aspect_ratio(aspect_ratio, "span")
    # glider2d.set_flat_area(area, "aspect_ratio")
    glider2d.shape.set_const_cell_dist()

def shape_plot(glider2d):
    shape = glider2d.shape.get_shape()
    line = list(shape.front) + list(shape.back[::-1]) + [shape.front[0]]
    x, y = zip(*line)
    middle = sum(y) / len(y)
    y = [i - middle for i in y]
    return x, y

def min_func(x):
    global count
    count += 1
    print("\niteration: " + str(count))
    _glider2d = deepcopy(glider2d)
    glider_set_controlpoint(_glider2d, x)
    glider3d = _glider2d.get_glider_3d()
    panels = parabem_Panels(glider3d,
                        midribs=0,
                        profile_numpoints=40,
                        symmetric=True,
                        num_average=0)
    case = Case(panels[1], panels[2])
    case.farfield = 5
    case.drag_calc = "trefftz"
    case.A_ref = 23

    case.create_wake(length=10000, count=5)
    case.v_inf = parabem.Vector3(8, 0, 1)
    case.trefftz_cut_pos = case.v_inf * 50
    alpha = v_inf_deg_range3(case.v_inf, 5, 9, 20)
    polars = case.polars(alpha)
    vtk_writer = CaseToVTK(case, "results/vtk_opt", suffix=str(count))
    vtk_writer.write_panels(data_type="point")
    cL = []
    cD = []
    for i in polars.values:
        cL.append(i.cL)
        cD.append(i.cD)
    return np.interp(0.6, cL, cD)


y_positions = np.linspace(-6, 10, 20)
cD = [min_func(i) for i in y_positions]
# plt.show()
plt.figure(figsize=(10, 4))
plt.plot(y_positions, cD, marker="x")
plt.xlabel('y-Position der Kontrollpunkte [m]', fontsize=15)
plt.ylabel('induzierter Widerstandsbeiwert $c_{Wi}$', fontsize=15)
plt.grid(True)
plt.savefig(results_path("vtk_opt", "induced_drag.png"),  bbox_inches='tight')
# # plt.show()
plt.close()
plt.figure(figsize=(10, 4))
plt.plot(y_positions, 0.6 / np.array(cD), marker="x")
plt.xlabel('y-Position der Kontrollpunkte [m]', fontsize=15)
plt.ylabel('Gleitzahl $\\epsilon$', fontsize=15)
plt.grid(True)
plt.savefig(results_path("vtk_opt", "glide.png"),  bbox_inches='tight')
plt.close()


best_ind = list(cD).index(min(cD))
best_y = y_positions[best_ind]

glider_set_controlpoint(glider2d, best_y)

with open(results_path("glider", "optimized.json"), "w") as _file:
    dump(glider2d, _file)

plt.figure(figsize=(10, 5))
plt.axes().set_aspect("equal", "datalim")
plt.grid(True)
plt.plot(*shape_plot(glider2d), color="0", label="optimal", linewidth=2)

glider_set_controlpoint(glider2d, y_positions[0])
with open(results_path("glider", "min.json"), "w") as _file:
    dump(glider2d, _file)
plt.plot(*shape_plot(glider2d), color="r", label="untere Grenze")


glider_set_controlpoint(glider2d, y_positions[-1])
with open(results_path("glider", "max.json"), "w") as _file:
    dump(glider2d, _file)

plt.plot(*shape_plot(glider2d), color="b", label="obere Grenze")
plt.legend()
plt.ylim((-2, 3))
plt.savefig(results_path("vtk_opt", "shapes.png"), bbox_inches='tight')