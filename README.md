# paraBEM
[![Build Status](https://travis-ci.org/looooo/paraBEM.svg?branch=master)](https://travis-ci.org/looooo/paraBEM)

<img src="./doc/latex_doc/Abbildungen/png/14_2_wake-rollup.png" alt="result" width="400"/>


paraBEM is a python module which provides methods to calculate the potential flow over 2D and 3D meshes with boundary-elements (panel-method)

the c++ code is wrapped with pybind11 to python. This allows fast computation (eigen) with a high-level interface.

## Dependencies
* C++11:
    - [eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page "Eigen")
    - [pybind11](https://pybind11.readthedocs.org/en/latest/ "pybind11")
    - [open_mp](http://openmp.org/wp/ "open mp")
    - [cmake](http://www.cmake.org/ "cmake")
* python:
    - [openglider](https://github.com/hiaselhans/OpenGlider "OpenGlider")
    - [numpy](http://www.numpy.org/ "mumpy")
    - [matplotlib](http://matplotlib.org/ "matplotlib")
    - [paraview](http://www.paraview.org/ "paraview")

## Install
### initialize
```bash
mkdir build && cd build
cmake ..
use cmake-gui to fix not find packages or wrong versions (eg. pybind11)
```

### build and install
if cmake doesn't complain, install the package. (the number at the end is the number of compile jobs)
```
sudo make install  -j2
```


## Using the package (with python)
```python
<<< import paraBEM
<<< help(paraBEM)
```


for more information look into the [tutorial](https://github.com/looooo/panel-methode/blob/master/doc/tutorial/tutorial.ipynb)
