# -*- coding: utf-8 -*-
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from openglider.jsonify import load
from openglider.utils.distribution import Distribution
from openglider.glider.in_out.export_3d import ppm_Panels

import ppm
from ppm.pan3d import DirichletDoublet0Source0Case3 as Case
from ppm.vtk_export import CaseToVTK
from ppm.utils import v_inf_deg_range3, check_path

with open("glider/referenz_schirm_berg.json") as _file:
    parametricGlider = load(_file)["data"]
    parametricGlider.shape.set_const_cell_dist()
    glider = parametricGlider.get_glider_3d()

_, panels, trailing_edge = ppm_Panels(
    glider,
    midribs=0,
    profile_numpoints=20,
    distribution=Distribution.cos_2_distribution,
    num_average=0,
    symmetric=False)

v_inf = [8, 0, 1]

case = Case(panels, trailing_edge)
case.mom_ref_point = ppm.Vector3(1.25, 0, -6)
case.A_ref = parametricGlider.shape.area
case.farfield = 5
case.drag_calc = "on_body"
case.trefftz_cut_pos = case.v_inf * 800

case.v_inf = ppm.Vector3(v_inf)
case.create_wake(length=10000, count=4)    # length, count
polars = case.polars(v_inf_deg_range3(case.v_inf, -5, 12, 5))


vtk_writer = CaseToVTK(case, "results/vtk_glider_case")
vtk_writer.write_panels(data_type="cell")
vtk_writer.write_wake_panels()
vtk_writer.write_body_stream(panels, 100)

p = []
alpha = []
for i in polars.values:
    alpha.append(np.rad2deg(i.alpha))
    p.append((i.cL, i.cD, i.cP, i.cop.x, i.cop.z))
plt.figure(figsize=(10, 4))
# plt.title(u"Beiwerte bei variiertem Anströmwinkel")
plt.ylabel(u"$\\alpha$", rotation=0, fontsize=20)
plt.xlabel(u"$c_{Wi}$, $c_N$, $c_A$", rotation=0, fontsize=20)

plt.plot(np.array(p).T[0], alpha, label=u"$c_A$")
plt.plot(np.array(p).T[1], alpha, label=u"$c_{Wi}$")
plt.plot(-np.array(p).T[2], alpha, label=u"- $c_N$")
plt.legend(fontsize=15)
plt.grid()
plt.savefig(check_path("results/cLcDcM.png"), bbox_inches="tight")
print(case.center_of_pressure)


p = []
alpha = []
for i in polars.values:
    alpha.append(np.rad2deg(i.alpha))
    p.append((i.cL, i.cD, i.cP, i.cop.x, i.cop.z))
plt.figure(figsize=(10, 4))
# plt.title(u"Beiwerte bei variiertem Anströmwinkel")
plt.ylabel(u"$\\alpha$", rotation=0, fontsize=20)
plt.xlabel(u"$c_{Wi}$, $c_N$, $c_A$", rotation=0, fontsize=20)

plt.plot(np.array(p).T[0], alpha, label=u"$c_A$")
plt.plot(np.array(p).T[1], alpha, label=u"$c_{Wi}$")
plt.plot(-np.array(p).T[2], alpha, label=u"- $c_N$")
plt.legend(fontsize=15)
plt.grid()
plt.savefig(check_path("results/cLcDcM.png"), bbox_inches="tight")
