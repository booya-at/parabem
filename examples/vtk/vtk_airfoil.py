import numpy as np

import paraBEM
from paraBEM.vtk_export import VtkWriter
from paraBEM.pan2d import DirichletDoublet0Case2 as Case
from paraBEM.airfoil import Airfoil
from paraBEM.utils import check_path

airfoil = Airfoil.trefftz_kutta(-0.1 + 0.00j, np.deg2rad(10))
airfoil.numpoints = 100

coordinates_1 = np.array(airfoil.coordinates[:-1])
coordinates_1 = list(map(list, coordinates_1))
pan_vectors_1 = list(map(paraBEM.PanelVector2, coordinates_1))
pan_vectors_1[0].wake_vertex = True

coordinates_2 = np.array(airfoil.coordinates[:-1]) + np.array([2.2, -0.3])
coordinates_2 = list(map(list, coordinates_2 * 0.5))
pan_vectors_2 = list(map(paraBEM.PanelVector2, coordinates_2))
pan_vectors_2[0].wake_vertex = True
panels = []
panels += [paraBEM.Panel2([pan_vec, pan_vectors_1[i+1]])
    for i, pan_vec 
    in enumerate(pan_vectors_1[:-1])]
panels.append(paraBEM.Panel2([pan_vectors_1[-1], pan_vectors_1[0]]))

panels += [paraBEM.Panel2([pan_vec, pan_vectors_2[i+1]])
    for i, pan_vec 
    in enumerate(pan_vectors_2[:-1])]
panels.append(paraBEM.Panel2([pan_vectors_2[-1], pan_vectors_2[0]]))
case = Case(panels)
case.v_inf = paraBEM.Vector2(1, 0.05)
case.run()

# nx = 200
# ny = 200
# space_x = np.linspace(-0.1, 2, nx)
# space_y = np.linspace(-0.4, 0.4, ny)
# grid = [paraBEM.Vector2(x, y) for y in space_y for x in space_x]

# velocity = list(map(case.off_body_velocity, grid))
# pot = list(map(case.off_body_potential, grid))

# with open(check_path("results/airfoil_2d.vtk"), "w") as _file:
#     writer = VtkWriter()
#     writer.structed_grid(_file, "airfoil", [nx, ny, 1])
#     writer.points(_file, grid)
#     writer.data(_file, velocity, name="velocity", _type="VECTORS", data_type="POINT_DATA")
#     writer.data(_file, pot, name="pot", _type="SCALARS", data_type="POINT_DATA")

# with open(check_path("results/airfoil.vtk"), "w") as _file:
#     writer = VtkWriter()
#     writer.unstructed_grid(_file, "airfoil")
#     writer.points(_file, [[i[0], i[1], 0]for i in airfoil.coordinates])
#     writer.lines(_file, [range(len(airfoil.coordinates))])
