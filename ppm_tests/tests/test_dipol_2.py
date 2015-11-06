import ppm
from ppm.pan2d import doublet_2, doublet_2_v

t = ppm.Vector2(2, 2)
s = ppm.Vector2(0, 0)
d = ppm.Vector2(1, 0)
print(doublet_2_v(t, s, d))