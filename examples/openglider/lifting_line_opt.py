# -*- coding: utf-8 -*-
from __future__ import division
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from openglider.jsonify import load
from openglider.glider.in_out.export_3d import paraBEM_Panels
from openglider.utils.distribution import Distribution

import paraBEM
from paraBEM._paraBEM import Case3
from paraBEM.liftingline import LiftingLine
from paraBEM.pan3d import DirichletDoublet0Source0Case3
from paraBEM.utils import check_path


count = 0
#   load the glider
with open("results/glider/optimized.json") as _file:
    glider2d = load(_file)["data"]
    area = glider2d.shape.area
    aspect_ratio = glider2d.shape.aspect_ratio

glider3d = glider2d.get_glider_3d()
verts, panels, trailing_edge = paraBEM_Panels(glider3d,
                        midribs=0,
                        profile_numpoints=20,
                        symmetric=False,
                        distribution=Distribution.from_cos_2_distribution,
                        num_average=0)

case = DirichletDoublet0Source0Case3(panels, trailing_edge)
case.v_inf = paraBEM.Vector3(*glider2d.v_inf)
case.create_wake(1000, 100)
case.trefftz_cut_pos = case.v_inf * 20
case.run()
gamma_pan = -np.array([edge.vorticity for edge in case.trailing_edge])



######################LiftingLine############################
te_line = [paraBEM.Vector3(0, i.y, i.z) for i in trailing_edge]
lifting_line = LiftingLine(te_line)
lifting_line.v_inf = case.v_inf
lifting_line.solve_for_best_gamma(1)
gamma = [line.best_gamma for line in lifting_line.segments]
gamma_max = max(gamma)
gamma_pan *= gamma_max / max(gamma_pan)

line_segment_length = np.array([0] + [line.b for line in lifting_line.segments])
line_segment_length = line_segment_length.cumsum()
spw = max(line_segment_length)
spw_proj = abs(lifting_line.segments[0].v1.y) * 2
line_segment_length -= spw / 2


#####################Plottin#################################


x, y, z, nx, ny, nz = np.transpose(
    np.array([np.array([i.mids.x, i.mids.y, i.mids.z, i.n.x, i.n.y, i.n.z])
    for i in lifting_line.segments]))
l = np.sqrt(np.diff(y) ** 2 + np.diff(z) ** 2)
l = l.cumsum()
l = np.insert(l, 0, 0)
l -= max(l) / 2

gamma_el_1 = lambda y: gamma_max * (1 - (y / spw *2)**2)**(1/2)
gamma_el_2 = lambda y: gamma_max * (1 - (y / spw_proj *2)**2)**(1/2)

np_gamma_1 = np.vectorize(gamma_el_1)
np_gamma_2 = np.vectorize(gamma_el_2)

gamma_ell_1 = np_gamma_1(l)
gamma_ell_2 = np_gamma_2(y)

plt.figure(figsize=(10,6))

label = (
    u"Hinterkante des Schirms",
    u"Optimale Zirkulationsverteilung nach dem Traglinien-Verfahren",
    u"Zirkulationsverteilung für den optimierten Schirm",
    u"Elliptische Zirkulation für den flach ausgebreiteten Schirm")

plt.axes().set_aspect("equal", "datalim")
plt.ylim([-2, 4])
plt.plot(y,z)
plt.plot(y + ny * gamma * 10, z + nz * gamma * 10, color="magenta", marker="x")
plt.plot(y + ny * gamma_pan * 10, z + nz * gamma_pan * 10, color="green", marker="+")
plt.plot(y + ny * gamma_ell_1 * 10, z + nz * gamma_ell_1 * 10, color="black", marker="h")
plt.grid(True)
plt.legend(label)
plt.savefig(check_path("results/vtk_opt/circulation_comparison.png"))
