# -*- coding: utf-8 -*-
import numpy as np

import paraBEM
from paraBEM.vtk_export import VtkWriter
from paraBEM.pan2d import *
from paraBEM import PanelVector2, Vector2, Panel2
from paraBEM.airfoil.conformal_mapping import *
from paraBEM.pan2d import DirichletDoublet0Case2 as Case
from paraBEM.utils import check_path

airfoil = JoukowskyAirfoil(midpoint=-0.1 + 0.05j)

# airfoil = VanDeVoorenAirfoil(tau=np.deg2rad(20), epsilon=0.05)

# airfoil = TrefftzKuttaAirfoil(midpoint=-0.1 + 0.05j, tau=np.deg2rad(10))

num_pan = 30
alpha = np.deg2rad(10)

# panelmethode
coordiantes = list(zip(airfoil.coordinates(num=num_pan).real, airfoil.coordinates(num=num_pan).imag))
vertices = [PanelVector2(*v) for v in coordiantes[:-1]]
vertices[0].wake_vertex = True
panels = [Panel2([vertices[i], vertices[i + 1]]) for i in range(len(vertices[:-1]))]
panels.append(Panel2([vertices[-1], vertices[0]]))


case = Case(panels)
case.v_inf = Vector2(np.cos(alpha), np.sin(alpha))
case.run()

nx = 300
ny = 200

space_x = np.linspace(-3, 7, nx)
space_y = np.linspace(-2, 2, ny)
vec = lambda x: paraBEM.Vector2(x[0], x[1])
vec3 = lambda x: [x[0], x[1], 0]
grid = [[x, y, 0] for y in space_y for x in space_x]
_grid = list(map(vec, grid))

velocity = list(map(vec3, list(map(case.off_body_velocity, _grid))))
vel1 = [(i[0]**2 + i[1]**2)**(0.5) for i in velocity]
pot = list(map(case.off_body_potential, _grid))

with open(check_path("results/coordinates_not_aligned/field.vtk"), "w") as _file:
    writer = VtkWriter()
    writer.structed_grid(_file, "airfoil", [nx, ny, 1])
    writer.points(_file, grid)
    writer.data(_file, velocity, name="velocity", _type="VECTORS", data_type="POINT_DATA")
    writer.data(_file, pot, name="pot", _type="SCALARS", data_type="POINT_DATA")
    writer.data(_file, vel1, name="vel", _type="SCALARS", data_type="POINT_DATA")

with open(check_path("results/coordinates_not_aligned/airfoil.vtk"), "w") as _file:
    writer = VtkWriter()
    writer.unstructed_grid(_file, "airfoil")
    writer.points(_file, [[i[0], i[1], 0] for i in coordiantes])
    writer.lines(_file, [range(len(coordiantes))])