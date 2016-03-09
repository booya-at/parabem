#PM-PanelMethod
[![Build Status](https://travis-ci.org/looooo/panel-method.svg?branch=master)](https://travis-ci.org/looooo/panel-method)

PM is a package which provides some functionality to calculate the potentialflow with low order panel-methodes. The methodes are for 2d and 3d lifting and non lifting problems.


##PPM
PPM is the python binding for PM. it is used to for the tests and exdents the package with some extra functionality like airfoils, mesh import, export for post-processing.

##Dependencies
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

##Install
###initialize
```bash
mkdir build && cd build
cmake ..
use cmake-gui to fix not find packages or wrong versions (eg. boost-python)
```

### build and install
if cmake doesn't complain, install the package. (the number at the end is the number of compile jobs)
```
sudo make install  -j2
```


##Using the package (with python)
```python
<<< import PPM
<<< help(PPM)
```


for more information look into the [tutorial](https://github.com/looooo/panel-methode/blob/master/doc/tutorial/tutorial.ipynb)
