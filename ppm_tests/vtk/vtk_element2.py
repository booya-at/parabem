import ppm
from ppm.pan2d import doublet_2_0, source_2_0, doublet_2_0_v, source_2_0_v, doublet_2_1, doublet_2_1_v
from ppm.vtk_export import VtkWriter
import numpy as np
from ppm.utils import check_path

v1 = ppm.PanelVector2(-2, 0)
v2 = ppm.PanelVector2(2, 0)
panel = ppm.Panel2([v1, v2])

n = 500
space = np.linspace(-5, 5, n)

grid = [ppm.Vector2([x, y]) for y in space for x in space]
dub_vals = [doublet_2_0(target, panel) for target in grid]
dub_vals_lin_1 = [doublet_2_1(target, panel, True) for target in grid]
dub_vals_lin_2 = [doublet_2_1(target, panel, False) for target in grid]
src_vals = [source_2_0(target, panel) for target in grid]
dublinv1_vals = [doublet_2_1_v(target, panel, True) for target in grid]
dublinv2_vals = [doublet_2_1_v(target, panel, False) for target in grid]
dubv_vals = [doublet_2_0_v(target, panel) for target in grid]
srcv_vals = [source_2_0_v(target, panel) for target in grid]


writer = VtkWriter()
with open(check_path("results/element_2.vtk"), "w") as _file:
    writer.structed_grid(_file, "element_2", [n, n, 1])
    writer.points(_file, grid)
    writer.data(_file, dub_vals, name="doublet", _type="SCALARS", data_type="POINT_DATA")
    writer.data(_file, dub_vals_lin_1, name="doublet_lin_1", _type="SCALARS", data_type="POINT_DATA")
    writer.data(_file, dub_vals_lin_2, name="doublet_lin_2", _type="SCALARS", data_type="POINT_DATA")
    writer.data(_file, src_vals, name="source", _type="SCALARS", data_type="POINT_DATA")
    writer.data(_file, dublinv1_vals, name="doublet_lin_1_v", _type="VECTORS", data_type="POINT_DATA")
    writer.data(_file, dublinv2_vals, name="doublet_lin_2_v", _type="VECTORS", data_type="POINT_DATA")
    writer.data(_file, dubv_vals, name="doublet_v", _type="VECTORS", data_type="POINT_DATA")
    writer.data(_file, srcv_vals, name="source_v", _type="VECTORS", data_type="POINT_DATA")
