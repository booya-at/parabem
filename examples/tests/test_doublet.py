import ppm
from ppm import pan3d

pnt1 = ppm.PanelVector3(-0.4, -0.5, 0)
pnt2 = ppm.PanelVector3(0.5, -0.2, 0)
pnt3 = ppm.PanelVector3(0.3, 0.5, 0)
pnt4 = ppm.PanelVector3(-0.5, 0.5, 0)
pnt5 = ppm.PanelVector3(0.0, 0, 1)
pnt6 = ppm.PanelVector3(0.3, 0.3, 0.3)

source = ppm.Panel3([pnt1, pnt2, pnt3, pnt4])


checklist = [pnt1, pnt2, pnt3, pnt4, pnt5, pnt6, source.center]
for trg_pnt in checklist:
    dip = pan3d.doublet_3_0_sphere(trg_pnt, source)
    print(trg_pnt, ": doublet:", dip)
