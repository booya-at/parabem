from paraBEM.pan2d import doublet_2_1
import paraBEM
import numpy as np

v1 = paraBEM.PanelVector2(-1, 0)
v2 = paraBEM.PanelVector2(1,  0)
panel = paraBEM.Panel2([v2, v1])

vals = ([doublet_2_1(paraBEM.Vector2(x, 0), panel, True) for x in np.linspace(-2, 2, 20)])
print(vals)
