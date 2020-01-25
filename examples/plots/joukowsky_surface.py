# -*- coding: utf-8 -*-
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import numpy as np

from parabem.pan2d import *
from parabem import PanelVector2, Vector2, Panel2
from parabem.airfoil.conformal_mapping import *
from parabem.utils import check_path

# GEOMETRY:
airfoil = JoukowskyAirfoil(midpoint=-0.1 + 0.05j)
airfoil_label = "Joukowsky Profil (m=-0.1 + 0.05j)"

# airfoil = VanDeVoorenAirfoil(tau=np.deg2rad(15), epsilon=0.03)
# airfoil_label = "Trefftz-Kutta Profil (m=-0.1 + 0.05j, tau=30)"

# airfoil = TrefftzKuttaAirfoil(midpoint=-0.05 + 0.05j, tau=np.deg2rad(10))
# airfoil_label = "VanDeVooren Profil(tau=20, epsilon=0.05)"

# INPUTS:
num_konform = 60000
num_pan = 40
_alpha = 10
alpha = np.deg2rad(_alpha)

# KONFORMAL MAPPING:
vel = airfoil.surface_velocity(alpha, num=num_konform)
vel = np.sqrt(vel.imag ** 2 + vel.real ** 2)
cp = airfoil.surface_cp(alpha, num=num_konform)

# PANELMETHODE:
coordiantes = list(zip(airfoil.coordinates(num=num_pan).real, airfoil.coordinates(num=num_pan).imag))
vertices = [PanelVector2(*v) for v in coordiantes[:-1]]
vertices[0].wake_vertex = True
panels = [Panel2([vertices[i], vertices[i + 1]]) for i in range(len(vertices[:-1]))]
panels.append(Panel2([vertices[-1], vertices[0]]))

pan_vel = []
for Case in [NeumannDoublet0Case2,
             DirichletDoublet0Case2,
             DirichletDoublet0Source0Case2,
             DirichletDoublet1Case2]:
    case = Case(panels)
    case.v_inf = Vector2(np.cos(alpha), np.sin(alpha))
    case.run()
    pan_vel.append([i.velocity.norm() for i in panels])
    print(case.cl)
    del(case)


# VISUALISATION
pan_center_x = [i.center.x for i in panels]


fig=plt.figure(figsize=(10,6))
fig.suptitle('alpha = '+ str(_alpha) +
    ', Panelanzahl = '+ str(num_pan), fontsize=17)

plt.plot(airfoil.coordinates(num_konform).real, airfoil.coordinates(num_konform).imag, label=airfoil_label)


plt.plot(airfoil.x(num_konform), np.array(vel), label="Konforme Abbildung")

label = "Dirichlet-Randbedingung konstante Dipol-Panel"
plt.plot(pan_center_x, pan_vel[1], marker='+', ls=':', label=label)

label = "Dirichlet-Randbedingung konstante Dipol-Quellen-panel "
plt.plot(pan_center_x, pan_vel[2], marker='1', ls=':', label=label)

# label = "Neumann Randbedingung konstante Dipol-Panel"
# plt.plot(pan_center_x, pan_vel[0], marker='x', ls=':', label=label)

# label = "Dirichlet Randbedingung linearen Dipol-panel "
# plt.plot(pan_center_x, pan_vel[3], marker='1', ls=':', label=label)


# plt.plot(airfoil.x(num_konform), cp)
# plt.plot(pan_center_x, pan_cp, marker="x", ls='None')

plt.legend()
plt.axes().set_aspect("equal", "datalim")
plt.xlim([-1,1])
plt.ylim([-0.5,3])
plt.xlabel('$x$', fontsize=25)
plt.ylabel('$\\frac{u}{u_{\\infty}}$', fontsize=25 , rotation=0)
plt.grid(True)
# plt.savefig("results/2d/trefftz-kutta.png")
plt.savefig(check_path("results/2d/joukowsky.png"))
# plt.savefig("results/2d/vandevooren_100.png")
# plt.show()
