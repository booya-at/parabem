import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import ppm
from ppm.vtk_export import CaseToVTK
from ppm.pan2d import DirichletDoublet0Source0Case2 as Case
from ppm.airfoil import Airfoil


n_x = 2000
n_y = 10


a = Airfoil.import_from_dat("/home/lo/Copy/modellbau/profile/jwl/jwl044.dat")
a.numpoints = n_x

case = Case(a.panels)
case.v_inf = ppm.Vector2(1, 0.2)
case.run()

plt.plot(*zip(*[[pan.center.x, pan.cp] for pan in case.panels]), marker="x")
# plt.show()
