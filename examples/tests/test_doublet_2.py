import parabem
from parabem.pan2d import doublet_2, doublet_2_v

t = parabem.Vector2(2, 2)
s = parabem.Vector2(0, 0)
d = parabem.Vector2(1, 0)
print(doublet_2_v(t, s, d))