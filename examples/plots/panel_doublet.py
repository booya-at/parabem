import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import paraBEM
from paraBEM.pan3d import doublet_3_0_sphere as doublet
from paraBEM.utils import check_path

pnt1 = paraBEM.PanelVector3(-0.5, -0.5, 0)
pnt2 = paraBEM.PanelVector3(0.5, -0.5, 0)
pnt3 = paraBEM.PanelVector3(0.5, 0.5, 0)
pnt4 = paraBEM.PanelVector3(-0.5, 0.5, 0)

source = paraBEM.Panel3([pnt1, pnt2, pnt3, pnt4])

x = np.linspace(-4, 4, 100)
y = []
for xi in x:
    target1 = paraBEM.PanelVector3(xi, 0., 0.0)
    target2 = paraBEM.PanelVector3(xi, 0, 0.5)
    target3 = paraBEM.PanelVector3(xi, 0, 1.0)
    val1 = doublet(target1, source)
    val2 = doublet(target2, source)
    val3 = doublet(target3, source)
    y.append([val1, val2, val3])

fig = plt.figure()
ax1 = fig.add_subplot(131)
ax1.plot(x, y)


y = []
for xi in x:
    target1 = paraBEM.Vector3(0, xi, 0.0)
    target2 = paraBEM.Vector3(0, xi, 0.5)
    target3 = paraBEM.Vector3(0, xi, 1.0)
    val1 = doublet(target1, source)
    val2 = doublet(target2, source)
    val3 = doublet(target3, source)
    y.append([val1, val2, val3])

ax2 = fig.add_subplot(132)
ax2.plot(x, y)

y = []
for xi in x:
    target1 = paraBEM.Vector3(0., 0, xi)
    target2 = paraBEM.Vector3(0.5, 0, xi)
    target3 = paraBEM.Vector3(1.0, 0, xi)
    val1 = doublet(target1, source)
    val2 = doublet(target2, source)
    val3 = doublet(target3, source)
    y.append([val1, val2, val3])

ax3 = fig.add_subplot(133)
ax3.plot(x, y)
plt.savefig(check_path("results/3d/doublet.png"))
