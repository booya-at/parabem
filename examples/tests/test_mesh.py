import os
import parabem
from parabem import pan3d
from parabem.mesh import mesh_object

directory = os.path.dirname(__file__)
mesh = mesh_object.from_OBJ(os.path.join(directory, "..", "mesh", "wing_lift.obj"))

case = pan3d.DirichletDoublet0Case3(mesh.panels, mesh.trailing_edges)
case.v_inf = parabem.Vector3(10, 0, 0)
case.create_wake(length=100, count=2)
case.run()
