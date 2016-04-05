from openglider.glider.in_out.export_3d import paraBEM_Panels
from openglider.jsonify import load
import paraBEM
from paraBEM.pan3d import DirichletDoublet0Source0Case3 as Case
from paraBEM.vtk_export import CaseToVTK

with open("../openglider/glider/referenz_schirm_berg.json") as _file:
    glider = load(_file)["data"]
    panels = paraBEM_Panels(glider.get_glider_3d() , 0, 30, 0, False)

case = Case(panels[1], panels[2])
case.v_inf = paraBEM.Vector(glider.v_inf)
case.drag_calc = "on_body"
case.A_ref = glider.shape.area

#wake parallel to the x axis
print("*    wake parallel to the x-axis")
case.create_wake(50, 100, paraBEM.Vector3(1., 0., 0.))
case.run()
writer = CaseToVTK(case, "results/wake/x-parallel")
writer.write_panels()
writer.write_wake_panels()
print("F_l = ", case.cL * case.A_ref * 1.2 * 11**2 / 2.,
      "F_w = ", case.cD * case.A_ref * 1.2 * 11**2 / 2.)


#wake parallel to v_inf
print("*    wake parallel to v_inf")
case.create_wake(50, 100, case.v_inf)
case.run()
writer = CaseToVTK(case, "results/wake/v_inf-parallel")
writer.write_panels()
writer.write_wake_panels()
print("F_l = ", case.cL * case.A_ref * 1.2 * 11**2 / 2.,
      "F_w = ", case.cD * case.A_ref * 1.2 * 11**2 / 2.)


#wake rollup to the x axis
print("*    rolled up wake")
case.relax_wake(3, 1)
case.run()
writer = CaseToVTK(case, "results/wake/roll_up")
writer.write_panels()
writer.write_wake_panels()
print("F_l = ", case.cL * case.A_ref * 1.2 * 11**2 / 2.,
      "F_w = ", case.cD * case.A_ref * 1.2 * 11**2 / 2.)
