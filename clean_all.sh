#!/usr/bin/bash

BASEDIR=$(dirname $0)
cd $BASEDIR/parabem_tests

rm plots/results -r
rm vtk/results -r
rm openglider/results -r