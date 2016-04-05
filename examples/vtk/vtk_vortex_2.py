import numpy as np
import paraBEM
from paraBEM.pan2d import vortex_2, vortex_2_v
from paraBEM.vtk_export import VtkWriter
from paraBEM.utils import check_path

source = paraBEM.Vector2(0, 0)
direction = paraBEM.Vector2(1, 0)
n = 20

a = np.linspace(-10, 10, n).tolist()
b = [paraBEM.Vector2(i, j) for i in a for j in a]

vel = [vortex_2_v(target, source) for target in b]
pot = [vortex_2(target, source, direction) for target in b]

writer = VtkWriter()
with open(check_path("results/vortex_2_field.vtk"), "w") as _file:
    writer.structed_grid(_file, "vortex", [n, n, 1])
    writer.points(_file, b)
    writer.data(_file, vel, name="velocity", _type="VECTORS", data_type="POINT_DATA")
    writer.data(_file, pot, name="potential", _type="SCALARS", data_type="POINT_DATA")
