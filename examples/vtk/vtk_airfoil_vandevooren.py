import numpy as np

from paraBEM.airfoil.conformal_mapping import VanDeVoorenAirfoil
from paraBEM.vtk_export import VtkWriter
from paraBEM.utils import check_path

####################################################
#      analytic solution of vandevooren-airfoils     #
####################################################

# -inputparameter for the vandevooren airfoil
alpha = np.deg2rad(10)  # alpha is the angle of attack in rad
tau = np.deg2rad(1)     # tau is the angle of the trailing edge
epsilon = 0.10          # epsilon is the thickness-parameter
num_x = 300             # number of plot-points in x-direction
num_y = 300             # number of plot-points in y-direction

# -create joukowsky object
airfoil = VanDeVoorenAirfoil(tau=tau, epsilon=epsilon)

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

# remove all points lying inside the circle with radius r
r = airfoil.radius
z_grid = [z for z in z_grid if abs(z) > r]
# ----------------------------------------------------------

# complex zeta-plane with mapped circle (=joukowsky airfoil)
# ----------------------------------------------------------
z_to_zeta = list(map(airfoil.zeta, z_grid))
zeta_vel = list(map(zeta_velocity, z_grid))
zeta_pot = list(map(potential, z_grid))
zeta_stream = list(map(stream, z_grid))

airfoil = list(map(complex_to_3vec, airfoil.coordinates()))
# ----------------------------------------------------------



with open(check_path("results/conformal_mapping/vandevooren_zeta.vtk"), "w") as _file:
    writer = VtkWriter()
    writer.structed_grid(_file, "zeta_plane", [num_x, num_y, 1])
    writer.points(_file, list(map(complex_to_3vec, z_to_zeta)))
    writer.data(_file, zeta_stream, name="stream", _type="SCALARS", data_type="POINT_DATA")
    writer.data(_file, zeta_pot, name="pot", _type="SCALARS", data_type="POINT_DATA")
    writer.data(_file, zeta_vel, name="velocity", _type="VECTORS", data_type="POINT_DATA")


with open(check_path("results/conformal_mapping/vandevooren_airfoil.vtk"), "w") as _file:
    writer = VtkWriter()
    writer.unstructed_grid(_file, "airfoil")
    writer.points(_file, airfoil)
    writer.lines(_file, [range(len(airfoil))])
