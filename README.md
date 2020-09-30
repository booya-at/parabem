# parabem
[![Build Status](https://travis-ci.org/booya-at/parabem.svg?branch=master)](https://travis-ci.org/booya-at/parabem)

<img src="./doc/latex_doc/Abbildungen/png/14_2_wake-rollup.png" alt="result" width="400"/>


parabem is a python module that provides methods to compute the potential flow over 2D and 3D objects with boundary-elements (panel-method).

the c++ code is wrapped with pybind11 to python. This allows fast computation (eigen) with a high-level interface.

## Dependencies
* C++11:
    - [eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page "Eigen")
    - [pybind11](https://pybind11.readthedocs.org/en/latest/ "pybind11")
    - [open_mp](http://openmp.org/wp/ "open mp")
    - [cmake](http://www.cmake.org/ "cmake")
* python:
    - [openglider](https://github.com/hiaselhans/OpenGlider "OpenGlider")
    - [numpy](http://www.numpy.org/ "numpy")
    - [scipy](https://www.scipy.org/ "scipy")
    - [matplotlib](http://matplotlib.org/ "matplotlib")
    - [paraview](http://www.paraview.org/ "paraview")

## Install
### initialize
```bash
mkdir build && cd build
cmake ..
cmake -DPYTHON_EXECUTABLE=/usr/local/bin/python3 ..
```

use cmake-gui to fix not found packages or wrong versions (eg. python, pybind11, ...)

### build and install
if cmake doesn't complain, install the package. (the number at the end is the number of compile jobs)
```
sudo make install  -j2
```


## Using the package (with python)
```python
<<< import parabem
<<< help(parabem)
```


for more information look into the [tutorial](https://github.com/looooo/panel-methode/blob/master/doc/tutorial/tutorial.ipynb)
