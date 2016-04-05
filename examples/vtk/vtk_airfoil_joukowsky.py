import numpy as np

from paraBEM.airfoil.conformal_mapping import JoukowskyAirfoil
from paraBEM.vtk_export import VtkWriter
from paraBEM.utils import check_path

####################################################
#      analytic solution of joukowsky-airfoils     #
####################################################

# -inputparameter for the joukowsky airfoil
midpoint = -0.2 + 0.0j  # m is the complex mid point of a circle passing (1 + 0j)
alpha = np.deg2rad(10)  # alpha is the angle of attack in rad
num_x = 300             # number of plot-points in x-direction
num_y = 300             # number of plot-points in y-direction

# -create joukowsky object
airfoil = JoukowskyAirfoil(midpoint)

# -helper functions


def complex_to_3vec(z):
    return [z.real, z.imag, 0]


def zeta_velocity(z):
    vel = airfoil.velocity(z, alpha)
    return [vel.real, -vel.imag, 0.]


def z_velocity(z):
    vel = airfoil.z_velocity(z, alpha)
    return [vel.real, -vel.imag, 0.]


def potential(z):
    pot = airfoil.potential(z, alpha)
    return pot.real


def stream(z):
    stream = airfoil.potential(z, alpha)
    return stream.imag


# complex z-plane with circle
# ----------------------------------------------------------
z_range_x = np.linspace(-3, 3, num_x)
z_range_y = np.linspace(-3, 3, num_y)
z_grid = [x + 1j * y for y in z_range_y for x in z_range_x]
z_vel = list(map(z_velocity, z_grid))
z_pot = list(map(potential, z_grid))
z_stream = list(map(stream, z_grid))

circle = list(map(complex_to_3vec, airfoil.circle()))
# ----------------------------------------------------------

# complex zeta-plane with mapped circle (=joukowsky airfoil)
# ----------------------------------------------------------
zeta_range_x = np.linspace(-3, 3, num_x)
zeta_range_y = np.linspace(-3, 3, num_y)
zeta_grid = [x + 1j * y for y in zeta_range_y for x in zeta_range_x]
zeta_to_z = list(map(airfoil.z, zeta_grid))
zeta_vel = list(map(zeta_velocity, zeta_to_z))
zeta_pot = list(map(potential, zeta_to_z))
zeta_stream = list(map(stream, zeta_to_z))

airfoil = list(map(complex_to_3vec, airfoil.coordinates()))
# ----------------------------------------------------------


with open(check_path("results/conformal_mapping/joukowsky_z.vtk"), "w") as _file:
    writer = VtkWriter()
    writer.structed_grid(_file, "z_plane", [num_x, num_y, 1])
    writer.points(_file, list(map(complex_to_3vec, z_grid)))
    writer.data(_file, z_stream, name="stream", _type="SCALARS", data_type="POINT_DATA")
    writer.data(_file, z_pot, name="pot", _type="SCALARS", data_type="POINT_DATA")
    writer.data(_file, z_vel, name="velocity", _type="VECTORS", data_type="POINT_DATA")

with open(check_path("results/conformal_mapping/joukowsky_zeta.vtk"), "w") as _file:
    writer = VtkWriter()
    writer.structed_grid(_file, "zeta_plane", [num_x, num_y, 1])
    writer.points(_file, list(map(complex_to_3vec, zeta_grid)))
    writer.data(_file, zeta_stream, name="stream", _type="SCALARS", data_type="POINT_DATA")
    writer.data(_file, zeta_pot, name="pot", _type="SCALARS", data_type="POINT_DATA")
    writer.data(_file, zeta_vel, name="velocity", _type="VECTORS", data_type="POINT_DATA")

with open(check_path("results/conformal_mapping/joukowsky_circle.vtk"), "w") as _file:
    writer = VtkWriter()
    writer.unstructed_grid(_file, "circle")
    writer.points(_file, circle)
    writer.lines(_file, [range(len(circle))])

with open(check_path("results/conformal_mapping/joukowsky_airfoil.vtk"), "w") as _file:
    writer = VtkWriter()
    writer.unstructed_grid(_file, "airfoil")
    writer.points(_file, airfoil)
    writer.lines(_file, [range(len(airfoil))])
