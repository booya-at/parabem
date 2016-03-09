import ppm
from ppm.pan3d import doublet_3_0_vsaero, doublet_3_0_vsaero_v, doublet_3_0_sphere
from ppm.vtk_export import VtkWriter
import numpy
from ppm.utils import check_path

v1 = ppm.PanelVector3(-0.5, -0.5, 0)
v2 = ppm.PanelVector3(0.5, -0.5, 0)
v3 = ppm.PanelVector3(0.5, 0.5, 0)
v4 = ppm.PanelVector3(-0.5,  0.5, 0)

p = ppm.Panel3([v1, v2, v3, v4])
# p = ppm.Panel3([v1, v2, v3])   # triangle
p.potential = 1.

n = 50

a = numpy.linspace(-1, 1, n).tolist()
b = [ppm.PanelVector3(i, j, k) for i in a for j in a for k in a]

pot = [doublet_3_0_sphere(i, p) for i in b]
vel = [doublet_3_0_vsaero_v(i, p) for i in b]

writer = VtkWriter()
with open(check_path("results/panel_influence.vtk"), "w") as _file:
    writer.structed_grid(_file, "duplet", [n, n, n])
    writer.points(_file, b)
    writer.data(_file, pot, name="potential", _type="SCALARS", data_type="POINT_DATA")
    writer.data(_file, vel, name="velocity", _type="VECTORS", data_type="POINT_DATA")
