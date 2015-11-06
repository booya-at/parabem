import numpy as np
import matplotlib.pyplot as plt

import ppm
from ppm.pan3d import src_3_0_vsaero, src_3_0_n0
from ppm.utils import check_path

pnt1 = ppm.PanelVector3(-0.5, -0.5, 0)
pnt2 = ppm.PanelVector3(0.5, -0.5, 0)
pnt3 = ppm.PanelVector3(0.5, 0.5, 0)
pnt4 = ppm.PanelVector3(-0.5, 0.5, 0)

source = ppm.Panel3([pnt1, pnt2, pnt3, pnt4])

x = np.linspace(0, 5, 500)
y = []
for xi in x:
    target = ppm.PanelVector3(xi, 0.0, 0.0)
    panel_infl = src_3_0_vsaero(target, source)
    point_infl = src_3_0_n0(target, source)
    y.append([panel_infl, point_infl, abs(panel_infl - point_infl) / abs(panel_infl)])

y = list(zip(*y))

plt.figure(figsize=(8,3))
plt.gcf().subplots_adjust(bottom=0.15)
plt.plot(x, y[0], label="Nahbereich")
plt.plot(x, y[1], label="Fernfeld")
plt.grid(True)
plt.legend()
plt.ylabel("Einfluss")
plt.xlabel("Abstand")
plt.ylim(-0.1, 0.6)

plt.savefig(check_path("results/3d/far_vs_near_source.png"))
plt.close()


plt.figure(figsize=(8,3))
plt.gcf().subplots_adjust(bottom=0.15)
plt.plot(x, y[2], label="Fehler durch Fernfeldmethode")
plt.yscale('log')
plt.grid(True)
plt.legend()
plt.ylim(-0.1, 0.6)

plt.xlabel("Abstand")
plt.ylabel("relativer Fehler")
plt.savefig(check_path("results/3d/far_error_source.png"))
