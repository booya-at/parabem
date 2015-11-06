import unittest
import numpy as np
from prettytable import PrettyTable
from ppm.pan3d import *
from ppm import Panel3, Vector3, PanelVector3


def short(numbers):
    if hasattr(numbers, "__getitem__"):
        return tuple([short(number) for number in numbers] )
    return "%0.3f" % numbers

class test_elements_3(object):
    def setUp(self, name="2D Elements"):
        self.p1 = PanelVector3(-1, -1, 0)
        self.p2 = PanelVector3(1, -1, 0)
        self.p3 = PanelVector3(1, 1, 0)
        self.p4 = PanelVector3(-1, 1, 0)
        self.panel = Panel3([self.p1, self.p2, self.p3, self.p4])

    @property
    def testpoints(self):
        return [self.p1, self.p2, self.p3, self.p4, self.panel.center]

    @property
    def test_lines(self):
        values = np.linspace(-3, 3, 10) 
        return [
            [Vector3(i, 0, 0) for i in values],
            [Vector3(i, i, 0) for i in values],
            [Vector3(0, i, 0) for i in values],
            [Vector3(0, 0, i) for i in values],
            [Vector3(self.p1.x, self.p1.y, i) for i in values]
        ]

    def test_element_phi(self, target=Vector3(0, 0, 0)):
        return 0

    def test_element_v(self, target=Vector3(0, 0, 0)):
        return 0, 0

    @property
    def name(self):
        return "name"

    def test_values(self):            
        x = PrettyTable([self.__class__.__name__, "influence", "influence derivative"])
        for target in self.testpoints:
            x.add_row(short([
                target,
                self.test_element_phi(target),
                self.test_element_v(target)]))
        print(x)
        for targets in self.test_lines:
            x = PrettyTable([self.__class__.__name__, "influence", "influence derivative"])
            for target in targets:
                x.add_row(short([
                    target,
                    self.test_element_phi(target),
                    self.test_element_v(target)]))
            print(x)

class src_3_0_vsaero_test(test_elements_3, unittest.TestCase):
    def test_element_phi(self, target=Vector3(0, 0, 0)):
        return src_3_0_vsaero(target, self.panel)

    def test_element_v(self, target=Vector3(0, 0, 0)):
        return src_3_0_vsaero_v(target, self.panel)

class doublet_3_0_sphere_test(test_elements_3, unittest.TestCase):
    def test_element_phi(self, target=Vector3(0, 0, 0)):
        return doublet_3_0_sphere(target, self.panel)

    def test_element_v(self, target=Vector3(0, 0, 0)):
        return doublet_3_0_vsaero_v(target, self.panel)

class doublet_src_3_0_n0_test(test_elements_3, unittest.TestCase):
    def test_element_phi(self, target=Vector3(0, 0, 0)):
        return doublet_src_3_0_n0(target, self.panel)
        
if __name__ == "__main__":
    unittest.main()