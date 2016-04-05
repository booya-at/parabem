import numpy
import paraBEM
from paraBEM.pan3d import vortex_3_0_v
from paraBEM.vtk_export import VtkWriter
from paraBEM.utils import check_path


v1 = paraBEM.Vector3(-1, 0, 0)
v2 = paraBEM.Vector3(1, 0, 0)

n = 80

a = numpy.linspace(-10, 10, n).tolist()
b = [paraBEM.Vector3(i, j, k) for i in a for j in a for k in a]

vel = [vortex_3_0_v(v1, v2, i) for i in b]

writer = VtkWriter()
with open(check_path("results/vortex_field.vtk"), "w") as _file:
    writer.structed_grid(_file, "vortex", [n, n, n])
    writer.points(_file, b)
    writer.data(_file, vel, name="velocity", _type="VECTORS", data_type="POINT_DATA")
