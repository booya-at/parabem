"""
ppm (python panel method) = python bindings to pm (panel method)
=================================================================

This package provides 2d and 3d panelmethodes and other utils to
calculate potential flow problems. To see how it works try the
following example:


Code snippets are indicated by three greater-than signs:

>>> import ppm


Use the built-in ``help`` function to view a function's docstring:

>>> help(ppm)

The code below shows the way how the package works. The example calculate the
potential flow arround a cylinder computed with a Panelmethode which uses
Dirichlet boundary condition with constant doublet (dipol) panels.


1. create the geometry of the flow problem with numpy:

>>> import numpy as np
>>> phi = np.linspace(0, 2 * np.pi, 30 + 1)
>>> x = np.cos(phi)[:-1]
>>> y = np.sin(phi)[:-1]
>>> xy = np.transpose(np.array([x, y]))


2. now map the geometry to PanelVecvtors and create Panels from the PanelVectors:

>>> import ppm
>>> coordinates = [ppm.PanelVector2(*i) for i in xy]
>>> coordinates += [coordinates[0]]
>>> panels = [ppm.Panel2([vec, coordinates[i+1]]) for i, vec in enumerate(coordinates[:-1])]


3. setup the Case and run it:

>>> from ppm.pan2d import DirichletDoublet0Case2 as Case
>>> case = Case(panels)
>>> case.v_inf = ppm.Vector2(1, 0)
>>> case.run()


4. displaying the results using matplotlib:

>>> import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt

>>> x1 = [list(i.center) for i in case.panels]
>>> x2 = [[i.center.x, i.velocity.norm()] for i in case.panels]

>>> plt.axes().set_aspect("equal", "datalim")
>>> plt.grid=True
>>> plt.plot(*zip(*x1))
>>> plt.plot(*zip(*x2))
>>> plt.show()

for more information look into the ipython notebook (tutorial.ipynb) which is
located in the docs directory.
"""

from ._ppm import PanelVector2, PanelVector3, Panel3, Panel2
from ._ppm import vector2 as Vector2
from ._ppm import vector3 as Vector3
from .utils import Vector
