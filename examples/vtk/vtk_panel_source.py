import parabem
from parabem.pan3d import doublet_src_3_0_vsaero, src_3_0_vsaero_v
from parabem.vtk_export import VtkWriter
import numpy
from parabem.utils import check_path

v1 = parabem.PanelVector3(-0.5, -0.5, 0)
v2 = parabem.PanelVector3(0.5, -0.5, 0)
v3 = parabem.PanelVector3(0.5, 0.5, 0)
v4 = parabem.PanelVector3(-0.5, 0.5, 0)

p = parabem.Panel3([v1, v2, v3, v4])
p.potential = 1.

n = 50

a = numpy.linspace(-1, 1, n).tolist()
b = [parabem.PanelVector3(i, j, k) for i in a for j in a for k in a]

pot = [doublet_src_3_0_vsaero(i, p)[1] for i in b]
vel = [src_3_0_vsaero_v(i, p) for i in b]

writer = VtkWriter()
with open(check_path("results/src_panel.vtk"), "w") as _file:
    writer.structed_grid(_file, "duplet", [n, n, n])
    writer.points(_file, b)
    writer.data(_file, pot, name="potential", _type="SCALARS", data_type="POINT_DATA")
    writer.data(_file, vel, name="velocity", _type="VECTORS", data_type="POINT_DATA")
