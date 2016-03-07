#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <iostream>
#include "string.h"

#include <vector>

#include "vector.h"
#include "panel3.h"
#include "panel2.h"
#include "case3.h"
#include "case2.h"
#include "element_influence.h"
#include "lifting_line.h"

using std::cout;
using std::endl;

namespace py = pybind11;


PYBIND11_PLUGIN(_ppm) {
    py::module m("_ppm", "pybind11 example plugin");
    py::class_<PanelVector2, bases<Vector2>>("PanelVector2");
};

