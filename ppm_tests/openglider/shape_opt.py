from copy import deepcopy
import numpy as np
from matplotlib import pyplot as plt

from openglider.glider import glider_2d
from openglider.jsonify import dump, load
from openglider.glider.in_out.export_3d import ppm_Panels
from openglider.airfoil import Profile2D

import ppm
from ppm.pan3d import DirichletDoublet0Source0Case3 as Case
from ppm.vtk_export import CaseToVTK
from ppm.utils import check_path

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
    line = list(shape[1]) + list(shape[2][::-1]) + [shape[1][0]]
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
    print(_glider2d.aspect_ratio)
    print(_glider2d.flat_area)
    panels = ppm_Panels(glider3d,
                        midribs=0,
                        profile_numpoints=20,
                        symmetric=False,
                        distribution=Profile2D.nose_cos_distribution(0.5),
                        num_average=0)
    case = Case(panels[1], panels[2])
    case.farfield = 5
    case.drag_calc = "trefftz"
    case.A_ref = glider3d.area

    case.create_wake(length=10000, count=5)
    case.v_inf = ppm.Vector3(8, 0, 1)
    case.trefftz_cut_pos = case.v_inf * 50
    case.run()
    vtk_writer = CaseToVTK(case, "results/vtk_opt", suffix=str(count))
    vtk_writer.write_panels(data_type="point")
    return x, case.cL, case.cD * 20, glider3d.aspect_ratio, glider3d.projected_area


cL_cD_table = [min_func(i) for i in np.linspace(-3, 4.5, 10)]
x, cL, cD, ar, area = np.array(list(zip(*cL_cD_table)))


glide = cL / cD
aero_eff = (cL**2 / cD / ar / np.pi)
best_ind = list(aero_eff).index(max(aero_eff))
best_x = x[best_ind]
plt.figure(figsize=(10,4))
plt.plot(x, cL)
plt.plot(x, cD)
plt.xlabel('y-Position der Kontrollpunkte [m]', fontsize=15)
plt.ylabel('aerodynamische Effizienz', fontsize=15)
plt.grid(True)
plt.savefig(check_path('results/vtk_opt/aero_eff.png'))
plt.show()
plt.close()

# export best, min, max

glider_set_controlpoint(glider2d, best_x)
with open(check_path("results/glider/optimized.json"), "w") as _file:
    dump(glider2d, _file)

plt.figure(figsize=(10, 4))
plt.axes().set_aspect("equal", "datalim")
plt.grid(True)
plt.plot(*shape_plot(glider2d), color="0", label="optimal")

glider_set_controlpoint(glider2d, -3)
with open(check_path("results/glider/min.json"), "w") as _file:
    dump(glider2d, _file)
plt.plot(*shape_plot(glider2d), color="r", label="untere Grenze")


glider_set_controlpoint(glider2d, 4.5)
with open(check_path("results/glider/max.json"), "w") as _file:
    dump(glider2d, _file)

plt.plot(*shape_plot(glider2d), color="b", label="obere Grenze")
plt.legend()
plt.savefig(check_path('results/vtk_opt/shapes.png'))
