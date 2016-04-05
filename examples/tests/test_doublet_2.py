import paraBEM
from paraBEM.pan2d import doublet_2, doublet_2_v

t = paraBEM.Vector2(2, 2)
s = paraBEM.Vector2(0, 0)
d = paraBEM.Vector2(1, 0)
print(doublet_2_v(t, s, d))