from __future__ import division

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import parabem
from parabem.pan3d import vortex_3_0_half_infinity_v
from parabem.utils import check_path

# infl = -1 / 4 / pi / r on a line passing through
# v1 and normal to the direction t of the vortex line
v1 = parabem.Vector3(0, 0, 0)
t = parabem.Vector3(1, 0, 0)
infl_f_test = lambda y: vortex_3_0_half_infinity_v(v1, t, parabem.Vector3(y, 1, 0)).z
infl_f = lambda y: - 1 / 4 / np.pi / y
y = [i for i in np.linspace(-2, 10, 100)]
infl_test = list(map(infl_f_test, y))
infl = list(map(infl_f, y))

# plt.plot(y, infl)
plt.plot(y, infl_test)
plt.savefig(check_path("results/3d/vortex_3_plot.png"))