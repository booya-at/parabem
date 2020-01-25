import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import parabem
from parabem.pan3d import src_3_0_n0
from parabem.utils import check_path

pnt1 = parabem.PanelVector3(-0.5, -0.5, 0)
pnt2 = parabem.PanelVector3(0.5, -0.5, 0)
pnt3 = parabem.PanelVector3(0.5, 0.5, 0)
pnt4 = parabem.PanelVector3(-0.5, 0.5, 0)

source = parabem.Panel3([pnt1, pnt2, pnt3, pnt4])

x = np.arange(-4, 4, 0.01)
y = []
for xi in x:
    trg1 = parabem.Vector3(xi, 0, 0.0)
    trg2 = parabem.Vector3(xi, 0, 0.5)
    trg3 = parabem.Vector3(xi, 0, 1)
    val1 = src_3_0_n0(trg1, source)
    val2 = src_3_0_n0(trg2, source)
    val3 = src_3_0_n0(trg3, source)
    y.append([val1, val2, val3])

fig = plt.figure()
ax1 = fig.add_subplot(131)
fig.gca().set_ylim([-5, 1])
ax1.plot(x, y)

y = []
for xi in x:
    trg1 = parabem.Vector3(0, xi, 0.0)
    trg2 = parabem.Vector3(0, xi, 0.5)
    trg3 = parabem.Vector3(0, xi, 1)
    val1 = src_3_0_n0(trg1, source)
    val2 = src_3_0_n0(trg2, source)
    val3 = src_3_0_n0(trg3, source)
    y.append([val1, val2, val3])
ax2 = fig.add_subplot(132)
fig.gca().set_ylim([-5, 1])
ax2.plot(x, y)

y = []
for xi in x:
    trg1 = parabem.Vector3(0, 0, xi)
    trg2 = parabem.Vector3(0.5, 0, xi)
    trg3 = parabem.Vector3(1, 0, xi)
    val1 = src_3_0_n0(trg1, source)
    val2 = src_3_0_n0(trg2, source)
    val3 = src_3_0_n0(trg3, source)
    y.append([val1, val2, val3])
ax3 = fig.add_subplot(133)
fig.gca().set_ylim([-5, 1])
ax3.plot(x, y)

plt.savefig(check_path("results/3d/source_far.png"))
