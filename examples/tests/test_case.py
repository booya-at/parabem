import parabem
from parabem import pan3d

v1 = parabem.PanelVector3(-0.5, -0.5, 0)
v2 = parabem.PanelVector3(0.5, -0.5, 0)
v3 = parabem.PanelVector3(0.5, 0.5, 0)
v4 = parabem.PanelVector3(-0.5, 0.5, 0)
v5 = parabem.PanelVector3(0.5,  1.5, 1)
v6 = parabem.PanelVector3(-0.5,  1.5, 1)

p1 = parabem.Panel3([v1, v2, v3, v4])
p2 = parabem.Panel3([v4, v3, v5, v6])

v_inf = parabem.Vector3(0, 1, 0)

c = pan3d.DirichletDoublet0Case3([p1, p2], [v1, v2, v3])
c.run()
