from parabem_cpp import (DirichletDoublet0Case3, DirichletDoublet0Source0Case3,
                      doublet_3_0_vsaero, doublet_3_0_n0, src_3_0_vsaero_v,
                      doublet_3_0_sphere, doublet_src_3_0_vsaero_v,
                      doublet_src_3_0_vsaero, doublet_src_3_0_n0, vortex_3_0_v,
                      doublet_3, doublet_3_v,
                      vortex_3_0_half_infinity_v, doublet_3_0_vsaero_v)


def src_3_0_vsaero(pan, position):
    return doublet_src_3_0_vsaero(pan, position)[1]


def src_3_0_n0(pan, position):
    return doublet_src_3_0_n0(pan, position)[1]