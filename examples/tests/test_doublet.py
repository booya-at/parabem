import paraBEM
from paraBEM import pan3d

pnt1 = paraBEM.PanelVector3(-0.4, -0.5, 0)
pnt2 = paraBEM.PanelVector3(0.5, -0.2, 0)
pnt3 = paraBEM.PanelVector3(0.3, 0.5, 0)
pnt4 = paraBEM.PanelVector3(-0.5, 0.5, 0)
pnt5 = paraBEM.PanelVector3(0.0, 0, 1)
pnt6 = paraBEM.PanelVector3(0.3, 0.3, 0.3)

source = paraBEM.Panel3([pnt1, pnt2, pnt3, pnt4])


checklist = [pnt1, pnt2, pnt3, pnt4, pnt5, pnt6, source.center]
for trg_pnt in checklist:
    dip = pan3d.doublet_3_0_sphere(trg_pnt, source)
    print(trg_pnt, ": doublet:", dip)
