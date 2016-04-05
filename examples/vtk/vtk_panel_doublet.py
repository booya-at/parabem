import paraBEM
from paraBEM.pan3d import doublet_3_0_vsaero, doublet_3_0_vsaero_v, doublet_3_0_sphere
from paraBEM.vtk_export import VtkWriter
import numpy
from paraBEM.utils import check_path

v1 = paraBEM.PanelVector3(-0.5, -0.5, 0)
v2 = paraBEM.PanelVector3(0.5, -0.5, 0)
v3 = paraBEM.PanelVector3(0.5, 0.5, 0)
v4 = paraBEM.PanelVector3(-0.5,  0.5, 0)

p = paraBEM.Panel3([v1, v2, v3, v4])
# p = paraBEM.Panel3([v1, v2, v3])   # triangle
p.potential = 1.

n = 50

a = numpy.linspace(-1, 1, n).tolist()
b = [paraBEM.PanelVector3(i, j, k) for i in a for j in a for k in a]

pot = [doublet_3_0_sphere(i, p) for i in b]
vel = [doublet_3_0_vsaero_v(i, p) for i in b]

writer = VtkWriter()
with open(check_path("results/panel_influence.vtk"), "w") as _file:
    writer.structed_grid(_file, "duplet", [n, n, n])
    writer.points(_file, b)
    writer.data(_file, pot, name="potential", _type="SCALARS", data_type="POINT_DATA")
    writer.data(_file, vel, name="velocity", _type="VECTORS", data_type="POINT_DATA")
