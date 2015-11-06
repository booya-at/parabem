import numpy as np
import matplotlib.pyplot as plt
import ppm
from ppm import pan3d
from ppm.mesh import mesh_object
from ppm.utils import check_path

# plots the values of the numerical solution as a scatter plot.

# running the case
mesh = mesh_object.from_OBJ("../mesh/sphere_low_tri.obj")
case = pan3d.DirichletDoublet0Case3(mesh.panels)
case.v_inf = ppm.Vector3(1, 0, 0)
case.run()

#  plotting the analytic solution
pot_f = lambda z: z * 3/2
pot_n_x, pot_n_y = np.array([[i.center.x, i.potential] for i in case.panels]).T

pot_x = np.cos(np.linspace(0, np.pi, 100))
pot_y = list(map(pot_f, pot_x))


vel_f = lambda z: np.sin(np.arccos(z)) * 3/2
vel_n_x, vel_n_y = np.array([[i.center.x, i.velocity.norm()] for i in case.panels]).T
vel_x = np.cos(np.linspace(0, np.pi, 100))
vel_y = list(map(vel_f, vel_x))


cp_f = lambda z: 1 - 9./4. * np.sin(np.arccos(z)) ** 2
cp_n_x, cp_n_y = np.array([[i.center.x, i.cp] for i in case.panels]).T
cp_x = np.cos(np.linspace(0, np.pi, 100))
cp_y = list(map(cp_f, cp_x))

plt.rcParams['figure.figsize'] = 10, 5

ax = plt.gca()
ax.set_xlabel("$\\frac{x}{r_K}$", fontsize=30, labelpad=20)
ax.set_ylabel("$\\frac{\phi}{u_{\infty} r_K}$", fontsize=30, rotation=0, labelpad=20)
plt.plot(pot_x, pot_y, label="Analytisch")
plt.title("Potential", fontsize=20)
plt.scatter(pot_n_x, pot_n_y, marker="+", label="Panelmethode")
plt.legend(fontsize=12)
plt.grid()
plt.savefig(check_path("results/3d/sphere_error_pot.png"), bbox_inches='tight', dpi=100)
plt.close()

ax = plt.gca()
ax.set_xlabel("$\\frac{x}{r_K}$", fontsize=30, labelpad=20)
ax.set_ylabel("$\\frac{u}{u_{\infty}}$", fontsize=30, rotation=0, labelpad=20)
plt.title("Geschwindigkeit", fontsize=20)
plt.plot(vel_x, vel_y, label="Analytisch")
plt.scatter(vel_n_x, vel_n_y, marker="+", label="Panelmethode")
plt.legend(fontsize=12)
plt.grid()
plt.savefig(check_path("results/3d/sphere_error_vel.png"), bbox_inches='tight', dpi=100)
plt.close()

ax = plt.gca()
ax.set_xlabel("$\\frac{x}{r_K}$", fontsize=30, labelpad=20)
ax.set_ylabel("$c_p$", fontsize=30, rotation=0, labelpad=20)
plt.title("Druckbeiwert", fontsize=20)
plt.plot(cp_x, cp_y, label="Analytisch")
plt.scatter(cp_n_x, cp_n_y, marker="+", label="Panelmethode")
plt.legend(fontsize=12)
plt.grid()
plt.savefig(check_path("results/3d/sphere_error_cp.png"), bbox_inches='tight', dpi=100)
plt.close()

# plt.savefig(check_path("results/3d/sphere_error_new.png"), bbox_inches='tight', dpi=100)
