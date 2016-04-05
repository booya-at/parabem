import numpy

import ppm
from ppm.pan2d import doublet_2_0, source_2_0, source_2_0_v
from ppm.vtk_export import VtkWriter
from ppm.utils import check_path

a = list(map(ppm.PanelVector2, [[-1, -1], [-1, 1]]))
b = list(map(ppm.PanelVector2, [[1, -1], [1, 1]]))
pan_a = ppm.Panel2(a)
pan_b = ppm.Panel2(b)

n = 100
space = numpy.linspace(-2, 2, n)
grid = [ppm.Vector2(x, y) for y in space for x in space]
pot = [
    source_2_0(v, pan_a) -
    source_2_0(v, pan_b)
    for v in grid]

vel = [
    source_2_0_v(v, pan_a) -
    source_2_0_v(v, pan_b)
    for v in grid]


writer = VtkWriter()
with open(check_path("results/parallel_flow.vtk"), "w") as _file:
    writer.structed_grid(_file, "doublet_2", [n, n, 1])
    writer.points(_file, grid)
    writer.data(_file, pot, name="potential", _type="SCALARS", data_type="POINT_DATA")
    writer.data(_file, vel, name="velocity", _type="VECTORS", data_type="POINT_DATA")
