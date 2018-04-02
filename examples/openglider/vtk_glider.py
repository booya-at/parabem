import numpy as np

from openglider.jsonify import load
from openglider.utils.distribution import Distribution
from openglider.glider.in_out.export_3d import paraBEM_Panels

import paraBEM
from paraBEM.pan3d import DirichletDoublet0Source0Case3 as Case
from paraBEM.vtk_export import CaseToVTK
from paraBEM.utils import v_inf_deg_range3

n_x = 30

with open("glider/referenz_schirm_berg.json") as _file:
    glider_2d = load(_file)["data"]
    glider_2d.shape.set_const_cell_dist()
    glider = glider_2d.get_glider_3d()

_, panels, trailing_edge = paraBEM_Panels(
    glider,
    midribs=3,
    profile_numpoints=n_x,
    distribution=Distribution.from_cos_2_distribution(n_x),
    num_average=0,
    symmetric=False)

v_inf = [8, 0, 1]

case = Case(panels, trailing_edge)
case.mom_ref_point = paraBEM.Vector3(1.25, 0, 0)
case.A_ref = glider_2d.shape.area
case.farfield = 5
case.drag_calc = "on_body"

case.v_inf = paraBEM.Vector3(v_inf)
case.create_wake(length=10000, count=4)    # length, count
polars = case.polars(v_inf_deg_range3(case.v_inf, -5, 15, 50))


vtk_writer = CaseToVTK(case, "results/vtk_glider_case")
vtk_writer.write_panels(data_type="cell")
# vtk_writer.write_wake_panels()
vtk_writer.write_body_stream(panels, 100)