from distutils.core import setup, Extension


extra_link_args = ['-lgomp']

include_dirs = ['/usr/include/eigen3',
                '/usr/local/include/eigen3',
                'src/headers']

src = ["src/python/paraBEM_ext.cpp",      "src/element_influence.cpp",
       "src/panel2.cpp",       "src/panel3.cpp",
       "src/case2.cpp",        "src/case3.cpp",
       "src/lifting_line.cpp"]

headers = ["src/headers/element_influence.h",
           "src/headers/vector.h",         "src/headers/panel2.h",
           "src/headers/panel3.h",         "src/headers/case2.h",
           "src/headers/case3.h",          "src/headers/lifting_line.h"]

files = ["paraBEM", "paraBEM.airfoil", "paraBEM.liftingline", "paraBEM.mesh",
         "paraBEM.pan2d", "paraBEM.pan3d", "paraBEM.utils", "paraBEM.vtk_export"]


setup(name='paraBEM._paraBEM',
      version='0.0.1',
      author='looooo',
      requires='eigen',
      author_email='sppedflyer@gmail.com',
      url="https://github.com/looooo/panelmethod",
      description='Wrap paraBEM using pybind11',
      packages=files,
      ext_modules=[Extension('paraBEM._paraBEM',
                   sources=src,
                   include_dirs=include_dirs,
                   extra_compile_args=['-std=c++11', '-fopenmp'],
                   extra_link_args=extra_link_args)])
