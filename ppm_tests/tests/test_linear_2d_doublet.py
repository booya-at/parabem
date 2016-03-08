from ppm.pan2d import doublet_2_1
import ppm
import numpy as np

v1 = ppm.PanelVector2(-1, 0)
v2 = ppm.PanelVector2(1,  0)
panel = ppm.Panel2([v2, v1])

vals = ([doublet_2_1(ppm.Vector2(x, 0), panel, True) for x in np.linspace(-2, 2, 20)])
print(vals)
