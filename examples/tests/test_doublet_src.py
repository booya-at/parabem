import numpy as np

import paraBEM
from paraBEM import pan3d

v1 = paraBEM.PanelVector3(-0.5,-0.5, 0)
v2 = paraBEM.PanelVector3( 0.5,-0.5, 0)
v3 = paraBEM.PanelVector3( 0.5, 0.5, 0)
v4 = paraBEM.PanelVector3(-0.5, 0.5, 0)

p = paraBEM.Panel3([v1, v2, v3, v4])

v5 = paraBEM.Vector3(0.5, 0, 0.5)
v6 = paraBEM.Vector3(0.4999, 0 , 0)
v7 = paraBEM.Vector3(0.5, 0.0 , 0)
v8 = paraBEM.Vector3(0.5001, 0 , 0)

checklist = [v1, v2, v3, v4, v5, v6, v7, v8, p.center]
for v in checklist:
    dip, src = pan3d.doublet_src_3_0_vsaero(v, p)
    print(v, ": doublet:", dip, "source:", src)
