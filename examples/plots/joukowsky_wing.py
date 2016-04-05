import matplotlib
matplotlib.use('Agg')import matplotlib.pyplot as plt
import numpy as np
import paraBEM
from paraBEM.vtk_export import CaseToVTK
from paraBEM.pan3d import DirichletDoublet0Source0Case3 as Case
from paraBEM.airfoil import Airfoil

def rib3d(airfoil, y_pos):
    out = [paraBEM.PanelVector3(coo[0], y_pos, coo[1]) for coo in airfoil.coordinates[:-1]]
    out.append(out[0])
    out[0].wake_vertex = True
    return out


n_x = 50
n_y = 10


a = Airfoil.joukowsky(-0.01+1j)
a.numpoints = n_x
print(a.coordinates)
ribs = [rib3d(a, y) for y in np.linspace(-5, 5, n_y)]
panels = []
for i in range(n_y)[:-1]:
    for j in range(n_x):
        panels.append(paraBEM.Panel3([ribs[i][j], ribs[i + 1][j], ribs[i + 1][j + 1], ribs[i][j + 1]]))
te = [rib[0] for rib in ribs]
print(te)
case = Case(panels, te)
case.farfield = 5
case.v_inf = paraBEM.Vector3(1, 0, 0.0)
case.create_wake(length=10000, count=3)    # length, count
case.run()

for i in range(n_y):
    plt.plot(*zip(*[[pan.center.x, pan.cp] for pan in case.panels[i * n_x : (i+1) * n_x]]),
        marker="x")
# plt.show()

plt.plot(*zip(*a.coordinates))
# plt.show()

vtk_writer = CaseToVTK(case, "results/joukowsky3_d")
vtk_writer.write_panels(data_type="point")