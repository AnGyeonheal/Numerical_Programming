/*-------------------------------------------------------------------------------\
@ Numerical Methods by Young-Keun Kim - Handong Global University

Author           : [An Gyeonheal]
ID               : [21900416]
Created          : 02-11-2022
Modified         : 02-11-2022
Language/ver     : C++ in MSVS2019

Description      : myNP.h
-------------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double trapz(double x[], double y[], int m);
double simpson13(double x[], double y[], int m);
double integral(double func(const double x), double a, double b, int n);

