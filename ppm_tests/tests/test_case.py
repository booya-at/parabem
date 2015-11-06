import ppm
from ppm import pan3d

v1 = ppm.PanelVector3(-0.5, -0.5, 0)
v2 = ppm.PanelVector3(0.5, -0.5, 0)
v3 = ppm.PanelVector3(0.5, 0.5, 0)
v4 = ppm.PanelVector3(-0.5, 0.5, 0)
v5 = ppm.PanelVector3(0.5,  1.5, 1)
v6 = ppm.PanelVector3(-0.5,  1.5, 1)

p1 = ppm.Panel3([v1, v2, v3, v4])
p2 = ppm.Panel3([v4, v3, v5, v6])

v_inf = ppm.Vector3(0, 1, 0)

c = pan3d.DirichletDoublet0Case3([p1, p2], [v1, v2, v3])
c.run()
