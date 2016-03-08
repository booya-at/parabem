import numpy as np

import ppm
from ppm import pan3d

v1 = ppm.PanelVector3(-0.5,-0.5, 0)
v2 = ppm.PanelVector3( 0.5,-0.5, 0)
v3 = ppm.PanelVector3( 0.5, 0.5, 0)
v4 = ppm.PanelVector3(-0.5, 0.5, 0)

p = ppm.Panel3([v1, v2, v3, v4])

v5 = ppm.Vector3(0.5,0,0.5)
v6 = ppm.Vector3(0.4999, 0 , 0)
v7 = ppm.Vector3(0.5, 0.0 , 0)
v8 = ppm.Vector3(0.5001, 0 , 0)

checklist = [v1, v2, v3, v4, v5, v6, v7, v8, p.center]
for v in checklist:
    dip, src = pan3d.doublet_src_3_0_vsaero(v, p)
    print(v, ": doublet:", dip, "source:", src)
