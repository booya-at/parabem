import paraBEM
from paraBEM.pan2d import DirichletDoublet0Case2 as Case
from paraBEM.airfoil import Airfoil

airfoil = Airfoil.vandevooren(0.5)


airfoil.numpoints = 40
case = Case(airfoil.panels)
case.v_inf = paraBEM.Vector2(1, 0.1)
case.run()

print(airfoil.coordinates)
