import paraBEM
from paraBEM.pan3d import doublet_3_0_vsaero_v, vortex_3_0_v


v1 = paraBEM.PanelVector3(-0.5, -0.5, 0)
v2 = paraBEM.PanelVector3(0.5, -0.5, 0)
v3 = paraBEM.PanelVector3(0.5, 0.5, 0)
v4 = paraBEM.PanelVector3(-0.5,  0.5, 0)

source = paraBEM.Panel3([v1, v2, v3, v4])
targets = [paraBEM.Vector3([0, 0, 0])]

for target in targets:
    print(doublet_3_0_vsaero_v(target, source))
    vortex_ring = (vortex_3_0_v(target, v2, v3),
                   vortex_3_0_v(target, v3, v4),
                   vortex_3_0_v(target, v4, v1),
                   vortex_3_0_v(target, v1, v2))
    print(vortex_ring)
