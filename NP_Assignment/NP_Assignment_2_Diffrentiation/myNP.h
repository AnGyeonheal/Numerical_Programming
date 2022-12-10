// 21900416 ¾È°ßÈú 2022-09-26

#ifndef		_MY_NP_H		// use either (#pragma once) or  (#ifndef ...#endif)
#define		_MY_NP_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void gradient1D(double _x[], double _y[], double dydx[], int m);
void gradientFunc(double func(const double x), double x[], double dydx[], int m);
void acceleration(double x[], double y[], double dy2dx[], int m);

#endif