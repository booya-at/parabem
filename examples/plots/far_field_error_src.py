# -*- coding: utf-8 -*-
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import parabem
from parabem.pan3d import src_3_0_vsaero, src_3_0_n0
from parabem.utils import check_path

pnt1 = parabem.PanelVector3(-0.5, -0.5, 0)
pnt2 = parabem.PanelVector3(0.5, -0.5, 0)
pnt3 = parabem.PanelVector3(0.5, 0.5, 0)
pnt4 = parabem.PanelVector3(-0.5, 0.5, 0)

source = parabem.Panel3([pnt1, pnt2, pnt3, pnt4])

x = np.linspace(0, 5, 500)
y = []
for xi in x:
    target = parabem.PanelVector3(xi, 0.0, 0.1)
    panel_infl = src_3_0_vsaero(target, source)
    point_infl = src_3_0_n0(target, source)
    y.append([panel_infl, point_infl, abs(panel_infl - point_infl)])

y = list(zip(*y))

plt.figure(figsize=(8,3))
plt.gcf().subplots_adjust(bottom=0.15)
plt.plot(x, y[0], label=u"Exakte Lösung")
plt.plot(x, y[1], label=u"Näherungslösung")
plt.grid(True)
plt.legend()
plt.ylabel("Einfluss")
plt.xlabel("x")
plt.ylim(None, 0.6)

plt.savefig(check_path("results/3d/far_vs_near_source.png"))
plt.close()


plt.figure(figsize=(8,3))
plt.gcf().subplots_adjust(bottom=0.15)
plt.plot(x, y[2], label="Fehler durch Fernfeldmethode")
plt.yscale('log')
plt.grid(True)
plt.legend()
plt.ylim(None, 0.6)

plt.xlabel("x")
plt.ylabel("Fehler")
plt.savefig(check_path("results/3d/far_error_source.png"))
