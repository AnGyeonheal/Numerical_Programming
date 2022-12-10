// 21900416 ¾È°ßÈú 2022-09-25

#ifndef		_MY_NP_H		// use either (#pragma once) or  (#ifndef ...#endif)
#define		_MY_NP_H
#define		PI		3.14159265358979323846264338327950288419716939937510582

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// Nonlinear Solver

double newtonRaphson(double dfunc(double x), double d2func(double x), double _x0, double _tol);

#endif

