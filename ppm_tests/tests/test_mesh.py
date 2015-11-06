import ppm
from ppm import pan3d
from ppm.mesh import mesh_object

mesh = mesh_object.from_OBJ("../mesh/wing_lift.obj")

case = pan3d.DirichletDoublet0Case3(mesh.panels, mesh.trailing_edges)
case.v_inf = ppm.Vector3(10, 0, 0)
case.create_wake(length=100, count=2)
case.run()
