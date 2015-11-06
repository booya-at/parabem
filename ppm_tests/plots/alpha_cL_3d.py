import numpy as np
from matplotlib import pyplot as plt

import ppm
from ppm import pan3d
from ppm.mesh import mesh_object
from ppm.vtk_export import CaseToVTK
from ppm.utils import check_path, vinf_deg_range3

mesh = mesh_object.from_OBJ("../mesh/wing_lift.obj")
alpha_0 = np.deg2rad(5)
v_inf = ppm.Vector3(np.cos(alpha_0), 0, np.sin(alpha_0))
v_inf_range = vinf_deg_range3(v_inf, -5, 10, 20)

fz = []
cmy = []
xcp = []
cL = []
cD = []

case = pan3d.DirichletDoublet0Source0Case3(mesh.panels, mesh.trailing_edges)
case.farfield = 5
case.create_wake(10, 30)

case.v_inf = v_inf
case.create_wake(length=10000, count=4)    # length, count
polars = case.polars(v_inf_range)

p = []
alpha = []
for i in polars.values:
    alpha.append(np.rad2deg(i.alpha))
    p.append((i.cL, i.cD, i.cP, i.cop.x, i.cop.z))
plt.figure(figsize=(10,4))
plt.ylabel(u"alpha")
plt.xlabel(u"cW, cN, cA")

plt.plot(np.array(p).T[0], alpha, label=u"cA")
plt.plot(np.array(p).T[1], alpha, label=u"cW_i")
plt.plot(-np.array(p).T[2], alpha, label=u"- cN")
plt.legend()
plt.grid()
plt.savefig(check_path("results/3d/cLcD.png"))
