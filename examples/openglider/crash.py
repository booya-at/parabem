import os

from openglider.jsonify import dump, load
from openglider.glider.in_out.export_3d import parabem_Panels

import parabem
from parabem.pan3d import DirichletDoublet0Source0Case3 as Case
from parabem.utils import v_inf_deg_range3

directory = os.path.dirname(__file__)

count = 0
with open(os.path.join(directory, "glider", "referenz_schirm_berg.json"), "r") as _file:
    glider2d = load(_file)["data"]
glider3d = glider2d.get_glider_3d()
panels = parabem_Panels(glider3d,
                    midribs=0,
                    profile_numpoints=40,
                    symmetric=True,
                    num_average=0)
case = Case(panels[1], panels[2])
case.farfield = 5
case.drag_calc = "trefftz"
case.A_ref = 23

case.create_wake(length=10000, count=5)
case.v_inf = parabem.Vector3(8, 0, 1)
case.trefftz_cut_pos = case.v_inf * 50
case.run()
case.create_wake(length=10000, count=5)
case.run()
case.relax_wake(3, 1)
case.run()