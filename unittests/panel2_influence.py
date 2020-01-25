import unittest
import numpy as np
from prettytable import PrettyTable
from parabem.pan2d import *
from parabem import Panel2, Vector2, PanelVector2


def short(numbers):
    if hasattr(numbers, "__getitem__"):
        return tuple([short(number) for number in numbers] )
    return "%0.3f" % numbers

class test_elements_2(object):
    def setUp(self, name="2D Elements"):
        self.p1 = PanelVector2(-1, 0)
        self.p2 = PanelVector2(1, 0)
        self.panel = Panel2([self.p1, self.p2])

    @property
    def testpoints(self):
        return [self.p1, self.p2, self.panel.center]

    @property
    def test_lines(self):
        values = np.linspace(-3, 3, 10) 
        return [
            [Vector2(i, 0) for i in values],
            [Vector2(i, i) for i in values],
            [Vector2(0, i) for i in values]
        ]

    def test_element_phi(self, target=Vector2(0, 0)):
        return 0

    def test_element_v(self, target=Vector2(0, 0)):
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

    # matplotlib and vtk!

class source_2_0_test(test_elements_2, unittest.TestCase):
    def test_element_phi(self, target=Vector2(0, 0)):
        return source_2_0(target, self.panel)

    def test_element_v(self, target=Vector2(0, 0)):
        return source_2_0_v(target, self.panel)

class doublet_2_0_test(test_elements_2, unittest.TestCase):
    def test_element_phi(self, target=Vector2(0, 0)):
        return doublet_2_0(target, self.panel)

    def test_element_v(self, target=Vector2(0, 0)):
        return doublet_2_0_v(target, self.panel)

class doublet_2_1_test(test_elements_2, unittest.TestCase):
    def test_element_phi(self, target=Vector2(0, 0)):
        return doublet_2_1(target, self.panel, True)

    def test_element_v(self, target=Vector2(0, 0)):
        return doublet_2_1_v(target, self.panel, True)
        
if __name__ == "__main__":
    unittest.main()