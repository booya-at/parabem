import os
from parabem.mesh import mesh_object
from parabem.pan3d import DirichletDoublet0Source0Case3 as Case
from parabem.vtk_export import CaseToVTK

directory = os.path.dirname(__file__)

mesh = mesh_object.from_OBJ(os.path.join(directory, "..", "mesh", "sphere_half.obj"))
for panel in mesh.panels:
    panel.set_symmetric()

case = Case(mesh.panels)
case.run()

writer = CaseToVTK(case, "results/symmetric_test")
writer.write_panels(data_type="point")
writer.write_field([-2, 2, 20], [-2, 2, 20], [-2, 2, 20])
