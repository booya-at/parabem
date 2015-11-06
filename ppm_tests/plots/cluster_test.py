from __future__ import division

import numpy as np
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt

import ppm
from ppm.pan3d import doublet_3_0_sphere, doublet_3_0_n0

x_range = [-1, 0, 1]
y_range = [-1, 0, 1]
points = [ppm.PanelVector3(x, y, 0) for y in y_range for x in x_range]
panel_indices = [
    [0, 1, 4, 3],
    [1, 2, 5, 4],
    [3, 4, 7, 6],
    [4, 5, 8, 7]]
panels = [ppm.Panel3([points[index] for index in indices]) for indices in panel_indices]
mue = [1, 1, 1, 100]
mue_mid = sum(mue) / 4
mid_pan = ppm.Panel3([
    points[0],
    points[2],
    points[8],
    points[6]])

def infl_near(point):
    out = 0
    target = ppm.Vector3(*point)
    for i, pan in enumerate(panels):
        out += mue[i] * doublet_3_0_sphere(target, pan)
    return out

def infl_far(point):
    return mue_mid * doublet_3_0_n0(ppm.Vector3(*point), mid_pan)

val_near = []
val_far = []
x_vals = np.linspace(-20, 20, 1000)
for x in x_vals:
    target = [x, x, 0.5]
    val_near.append(infl_near(target))
    val_far.append(infl_far(target))
val_near = np.array(val_near)
val_far = np.array(val_far)
error = (val_near - val_far) / val_near
plt.plot(x_vals, val_near, label="near")
plt.plot(x_vals, val_far, label="far")
plt.legend()

# plt.show()

plt.plot(x_vals, error)
# plt.show()

