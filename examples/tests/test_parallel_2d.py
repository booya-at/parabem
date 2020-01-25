import numpy as np

import parabem
from parabem.pan2d import doublet_2_0, doublet_2_0_v

v1 = parabem.PanelVector2(-1, 0)
v2 = parabem.PanelVector2(1, 0)
p = parabem.Panel2([v1, v2])

v = parabem.Vector2(0, 0)
print("doublet influence on panel center:")
print(doublet_2_0_v(v, p).y)
print("")

y = 3
pos = [[x, y] for x in np.linspace(-2, 2, 100)]
_pos = map(parabem.Vector2, pos)
vals = [doublet_2_0(i, p) for i in _pos]

print("x, y, doublet_value")
print("==============================\n")
print(np.array([[pos[0], pos[1], vals[i]] for i, pos in enumerate(pos)]))
