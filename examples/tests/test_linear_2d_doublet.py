from parabem.pan2d import doublet_2_1
import parabem
import numpy as np

v1 = parabem.PanelVector2(-1, 0)
v2 = parabem.PanelVector2(1,  0)
panel = parabem.Panel2([v2, v1])

vals = ([doublet_2_1(parabem.Vector2(x, 0), panel, True) for x in np.linspace(-2, 2, 20)])
print(vals)
