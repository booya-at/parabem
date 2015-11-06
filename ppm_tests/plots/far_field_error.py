# -*- coding: utf-8 -*-
from __future__ import division
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import ppm
from ppm.pan3d import doublet_3_0_n0, doublet_3_0_sphere, doublet_3_0_vsaero
from ppm.utils import check_path

pnt1 = ppm.PanelVector3(-0.5, -0.5, 0)
pnt2 = ppm.PanelVector3(0.5, -0.5, 0)
pnt3 = ppm.PanelVector3(0.5, 0.5, 0)
pnt4 = ppm.PanelVector3(-0.5, 0.5, 0)

source = ppm.Panel3([pnt1, pnt2, pnt3, pnt4])

x = np.linspace(-0, 5, 500)
y = []
for xi in x:
    target = ppm.PanelVector3(xi, 0., 0.1)
    panel_infl = doublet_3_0_vsaero(target, source)
    point_infl = doublet_3_0_n0(target, source)
    y.append([panel_infl, point_infl, abs(panel_infl - point_infl)])
y = list(zip(*y))


plt.figure(figsize=(8,3))
plt.gcf().subplots_adjust(bottom=0.15)
plt.ylabel("Einfluss")
plt.xlabel("x")
plt.plot(x, y[0], label=u"Exakte Lösung")
plt.plot(x, y[1], label=u"Näherungslösung")
plt.ylim(-0.1, 0.6)
plt.xlim(0, None)
plt.legend()
plt.grid(True)
plt.savefig(check_path("results/3d/far_vs_near_doublet.png"))
plt.close()

plt.figure(figsize=(8,3))
plt.gcf().subplots_adjust(bottom=0.15)
plt.ylabel("Fehler")
plt.xlabel("x")
plt.yscale('log')
plt.plot(x, y[2], label="Fehler durch Fernfeldmethode")
plt.ylim(None, 1)
plt.xlim(0, None)
plt.legend()
plt.grid(True)
plt.savefig(check_path("results/3d/far_error_doublet.png"))
