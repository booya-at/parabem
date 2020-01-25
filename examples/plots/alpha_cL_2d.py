import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import numpy as np

import parabem
from parabem.pan2d import DirichletDoublet0Source0Case2 as Case
from parabem.airfoil import Airfoil
from parabem.utils import check_path

a = Airfoil.trefftz_kutta(-0.1+0.01j, np.deg2rad(2), 51)

alpha_list = np.deg2rad(np.linspace(-15, 30, 30))
cl = []
cd = []
cm = []
xcp = []

for alpha in alpha_list:
    case = Case(a.panels)
    case.v_inf = parabem.Vector2(np.cos(alpha), np.sin(alpha))
    case.mom_ref_point = parabem.Vector2(-0, -3)
    case.run()
    cl.append(case.cl)
    cd.append(case.force.dot(case.v_inf) * 10)
    cm.append(case.cm)
    # xcp.append(case.center_of_pressure.x)
    del(case)

plt.plot(cl, alpha_list,  color="black", linestyle="-", label="cl")
plt.plot(cm, alpha_list, color="black", dashes=[8, 4, 2, 4, 2, 4], label="cm")
# plt.plot(xcp, alpha_list, color="black", linestyle="dashed", label="x(cm=0)")
plt.legend()
plt.plot(*zip(*a.coordinates), marker="*")
plt.grid(True)
# plt.axes().set_aspect("equal", "datalim")

axes = plt.gca()
axes.set_xlim([-0.25, 1.25])

plt.savefig(check_path("results/2d/alpha_ca_cm_xcp_2d.png"))
# plt.show()
