from copy import deepcopy
import numpy as np
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt

# from openglider.glider.parametric import ParametricGlider
from openglider.jsonify import dump, load
from openglider.glider.in_out.export_3d import paraBEM_Panels
from openglider.utils.distribution import Distribution

import paraBEM
from paraBEM.pan3d import DirichletDoublet0Source0Case3 as Case
from paraBEM.vtk_export import CaseToVTK
from paraBEM.utils import check_path, v_inf_deg_range3

count = 0
#   load the glider
with open("glider/referenz_schirm_berg.json") as _file:
    glider2d = load(_file)["data"]
    area = glider2d.shape.area
    aspect_ratio = glider2d.shape.aspect_ratio


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
    panels = paraBEM_Panels(glider3d,
                        midribs=0,
                        profile_numpoints=40,
                        symmetric=True,
                        distribution=Distribution.nose_cos_distribution(0.2),
                        num_average=0)
    case = Case(panels[1], panels[2])
    case.farfield = 5
    case.drag_calc = "trefftz"
    case.A_ref = 23

    case.create_wake(length=10000, count=5)
    case.v_inf = paraBEM.Vector3(8, 0, 1)
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
plt.savefig(check_path('results/vtk_opt/induced_drag.png'),  bbox_inches='tight')
# # plt.show()
plt.close()
plt.figure(figsize=(10, 4))
plt.plot(y_positions, 0.6 / np.array(cD), marker="x")
plt.xlabel('y-Position der Kontrollpunkte [m]', fontsize=15)
plt.ylabel('Gleitzahl $\\epsilon$', fontsize=15)
plt.grid(True)
plt.savefig(check_path('results/vtk_opt/glide.png'),  bbox_inches='tight')
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
