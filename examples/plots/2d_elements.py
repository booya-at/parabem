import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import numpy as np

import paraBEM
from paraBEM.pan2d import doublet_2_0, source_2_0, doublet_2_0_v
from paraBEM.utils import check_path

pnt1 = paraBEM.PanelVector2(-1, 0)
pnt2 = paraBEM.PanelVector2(1, 0)

source = paraBEM.Panel2([pnt1, pnt2])

y = np.linspace(-3, 3, 100)
val = [source_2_0(paraBEM.Vector2(yi, 8), source) for yi in y]
plt.plot(y, val)
val = [source_2_0(paraBEM.Vector2(yi, 0.01), source) for yi in y]
plt.plot(y, val)
val = [source_2_0(paraBEM.Vector2(yi, 0.0), source) for yi in y]
plt.plot(y, val)
val = [source_2_0(paraBEM.Vector2(yi, 3), source) for yi in y]
plt.plot(y, val)
plt.savefig(check_path("results/2d/source.png"))
plt.close()

y = np.linspace(-3, 3, 100)
val = [doublet_2_0(paraBEM.Vector2(yi, 7), source) for yi in y]
plt.plot(y, val)
val = [doublet_2_0(paraBEM.Vector2(yi, 0.01), source) for yi in y]
plt.plot(y, val)
val = [doublet_2_0(paraBEM.Vector2(yi, 0.0), source) for yi in y]
plt.plot(y, val)
val = [doublet_2_0(paraBEM.Vector2(yi, 3), source) for yi in y]
plt.plot(y, val)
plt.savefig(check_path("results/2d/doublet.png"))
plt.close()

y = np.linspace(-3, 3, 100)
val = [doublet_2_0_v(paraBEM.Vector2(yi, 7), source).x for yi in y]
plt.plot(y, val)
val = [doublet_2_0_v(paraBEM.Vector2(yi, 0.2), source).x for yi in y]
plt.plot(y, val)
val = [doublet_2_0_v(paraBEM.Vector2(yi, 3), source).x for yi in y]
plt.savefig(check_path("results/2d/doublet_v.png"))
