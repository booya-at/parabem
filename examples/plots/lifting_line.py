from __future__ import division
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import paraBEM

from paraBEM.liftingline import LiftingLine
from paraBEM.utils import check_path


# WingGeometry
spw = 2
numpos = 50
z_fac_1 = -0.3
z_fac_2 = -0.7
y = np.sin(np.linspace(0, np.pi/2, numpos)) * spw/2
x = [0. for _ in y]
z = [i**2 * z_fac_1 + i**6 * z_fac_2 for i in y]

mirror = lambda xyz: [xyz[0], -xyz[1], xyz[2]]
wing = list(zip(x, y, z))
wing = list(map(mirror, wing))[::-1] + list(wing)[1:]
wing = [paraBEM.Vector3(*i) for i in wing]

# LiftingLine
lifting_line = LiftingLine(wing)
lifting_line.v_inf = paraBEM.Vector3(1, 0, 0)
lifting_line.solve_for_best_gamma(1)
gamma = [i.best_gamma for i in lifting_line.segments]
gamma_max = max(gamma)

# Plot
gamma_el = lambda y: gamma_max * (1 - (y / spw * 2)**2)**(1 / 2)
mids = [[i.mids.x, i.mids.y, i.mids.z] for i in lifting_line.segments]
x, y, z = zip(*mids)

fig = plt.figure()
ax1 = fig.add_subplot(3, 1, 1)
ax1.plot(y, z)

ax2 = fig.add_subplot(3, 1, 2)
ax2.plot(y, x, marker="x")

ax3 = fig.add_subplot(3, 1, 3)
y_el = np.linspace(-1, 1, 400)
ax3.plot([-spw/2] + list(y) + [spw/2], [0] + gamma + [0], marker="x")
ax3.plot(y_el, list(map(gamma_el, y_el)))
plt.savefig(check_path("results/2d/liftingline.png"))

total = 0
for i in lifting_line.segments:
    total += i.lift_factor * i.best_gamma
print(total)
