#PM-PanelMethode
PM is a package which provides some functionality to calculate the potentialflow with low order panel-methodes. The methodes are for 2d and 3d lifting and non lifting problems.

##PPM
PPM is the python binding for PM. it is used to for the tests and exdents the package with some extra functionality like airfoils, mesh import, export for post-processing.

##Dependencies
* C++11:
    - [eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page "Eigen")
    - [boost python](http://www.boost.org/doc/libs/1_58_0/libs/python/doc/ "boost python")
    - [open_mp](http://openmp.org/wp/ "open mp)
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
```
###Options for cmake
if another python version should be used (eg. version 3.4):
```
cmake .. -Dpy=3.4
```
if the boost-python library isn't found try to locate it and use the last part of the file as option input:
```
cmake .. -Dboost=python-py27
```
to loacte the directory:
```
locate libboost |grep python
```

### build and install
if cmake doesn't complain, install the package. (the number at the end is the number of compile jobs)
```
sudo make install  -j2
```


##Using the package
```python
<<< import PPM
<<< help(PPM)
```


for more information look into this [tutorial](https://bytebucket.org/lorenz_l/pm_0.2/raw/bf3725b4e1976004171f2a8e2e167dad7dc5db3c/docs/tutorial_1.pdf)
