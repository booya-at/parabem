import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import parabem
from parabem import pan3d
from parabem.mesh import mesh_object
from parabem.utils import check_path

directory = os.path.dirname(__file__)
# plots the values of the numerical solution as a scatter plot.

# running the case
mesh = mesh_object.from_OBJ(os.path.join(directory, "..", "mesh", "sphere_low_tri.obj"))
case = pan3d.DirichletDoublet0Case3(mesh.panels)
case.v_inf = parabem.Vector3(1, 0, 0)
case.run()

#  plotting the analytic solution
fig = plt.figure()
pot_f = lambda z: z * 3/2
pot_n_x, pot_n_y = np.array([[i.center.x, i.potential] for i in case.panels]).T

pot_x = np.arange(-1., 1.001, 0.1)
pot_y = list(map(pot_f, pot_x))

ax1 = fig.add_subplot(311)
ax1.plot(pot_x, pot_y, c=[0, 0, 0], label="Analytisch")
ax1.scatter(pot_n_x, pot_n_y, marker="+", label="Panelmethode")
ax1.set_xlabel('x')
ax1.set_ylabel('Potential')
ax1.legend(loc="upper left", fontsize=8)


vel_f = lambda z: np.sin(np.arccos(z)) * 3/2
vel_n_x, vel_n_y = np.array([[i.center.x, i.velocity.norm()] for i in case.panels]).T
vel_x = np.arange(-1, 1, 0.01)
vel_y = list(map(vel_f, vel_x))

ax2 = fig.add_subplot(312)
ax2.plot(vel_x, vel_y, c=[0, 0, 0], label="Analytisch")
vel_scat = ax2.scatter(vel_n_x, vel_n_y, marker="+", label="Panelmethode")
ax2.set_xlabel('x')
ax2.set_ylabel('Geschwindigkeit')
ax2.legend(loc="upper left", fontsize=8)


cp_f = lambda z: 1 - 9./4. * np.sin(np.arccos(z)) ** 2
cp_n_x, cp_n_y = np.array([[i.center.x, i.cp] for i in case.panels]).T
cp_x = np.arange(-1, 1, 0.01)
cp_y = list(map(cp_f, cp_x))

ax3 = fig.add_subplot(313)
ax3.plot(cp_x, cp_y, c=[0, 0, 0], label="Analytisch")
cp_scat = ax3.scatter(cp_n_x, cp_n_y, marker="+", label="Panelmethode")
ax3.set_xlabel('x')
ax3.set_ylabel('Druckbeiwert')
ax3.legend(loc="upper left", fontsize=8)

plt.savefig(check_path("results/3d/sphere_error.png"), bbox_inches='tight', dpi=100)
