import os
import parabem
from parabem import pan3d
from parabem.mesh import mesh_object

directory = os.path.dirname(__file__)
mesh = mesh_object.from_OBJ(os.path.join(directory, "..", "mesh", "box_minimal.obj"))

case = pan3d.DirichletDoublet0Case3(mesh.panels)
case.v_inf = parabem.Vector3(1, 0, 0)

a = case.panels[0]
b = case.panels[1]

print(a.center, " ", a.n)
print(b.center, " ", b.n)
print(pan3d.doublet_3_0_vsaero(a.center, b))
