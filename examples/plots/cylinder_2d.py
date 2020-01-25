import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import numpy as np
import parabem
from parabem.pan2d import DirichletDoublet0Case2 as Case
from parabem.utils import check_path

# geometry
numpoints = 30
phi = np.linspace(0, 2 * np.pi, numpoints + 1)
x = np.cos(phi)[:-1]
y = np.sin(phi)[:-1]
xy = np.transpose(np.array([x, y]))

# mapping the geometry
coordinates = [parabem.PanelVector2(*i) for i in xy]
coordinates += [coordinates[0]]
panels = [parabem.Panel2([vec, coordinates[i+1]]) for i, vec in enumerate(coordinates[:-1])]

# setting up the case
case = Case(panels)
case.v_inf = parabem.Vector2(1, 0)
case.run()

# visualisation
x1 = [list(i.center) for i in case.panels]
x2 = [[i.center.x, i.velocity.norm()] for i in case.panels]

plt.axes().set_aspect("equal", "datalim")
plt.grid=True
plt.plot(*zip(*(x1 + [x1[0]])))
plt.plot(*zip(*x2))
plt.savefig(check_path("results/2d/cylinder_cp.png"))

# plt.show()
