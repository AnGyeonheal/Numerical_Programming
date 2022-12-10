/*-------------------------------------------------------------------------------\
@ Numerical Methods by Young-Keun Kim - Handong Global University

Author           : [An Gyeonheal]
ID               : [21900416]
Created          : 02-11-2022
Modified         : 24-11-2022
Language/ver     : C++ in MSVS2019

Description      : myNP.h
-------------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define		PI		3.14159265358979323846264338327950288419716939937510582
#define		EU		0
#define		EM		1
#define		RK2		2
#define		RK3     3

/*===================================TEST1================================================*/

// Q1
double q1_integral(double func(const double x), double a, double b, int n);

// Q2
double diffrentiation(double x, double q2_func(double _x));

// Q3
double steffensen(double func(double x), double x0, double tol);

/*===================================TEST1================================================*/

// TaylorSeries
double sinTaylor(double _x);
double sinTaylor2(double _x);
double sindTaylor(double _x);
double cosTaylor(double _x);
double power(double _x, int N);
double factorial(int _x);

// Bisection and TU_Newtonrapson
double bisection(double func(double x), float _a0, float _b0, float _tol);
double TU_newtonRaphson(double func(double x), double dfunc(double x), double _x0, double _tol);

// Nonlinear
double newtonRaphson(double dfunc(double x), double d2func(double x), double _x0, double _tol);

// Diffrential
void gradient1D(double _x[], double _y[], double dydx[], int m);
void gradientFunc(double func(const double x), double x[], double dydx[], int m);
void acceleration(double x[], double y[], double dy2dx[], int m);

// Integration
double trapz(double x[], double y[], int m);
double simpson13(double x[], double y[], int m);
double integral(double func(const double x), double a, double b, int n);

// ode 1
void odeEU(double myfunc(const double t, const double y), double y[], double t0, double tf, double h);
void odeEM(double myfunc(const double t, const double y), double y[], double t0, double tf, double h);
void odeRK2(double myfunc(const double t, const double y), double y[], double t0, double tf, double h, double y0);
void odeRK3(double myfunc(const double t, const double y), double y[], double t0, double tf, double h, double y0);
void ode(double myfunc(const double t, const double y), double y[], double t0, double tf, double h, int method);

// ode 2
void sys2RK2(void func(const double t, const double Y[], double dYdt[]), double y1[], double y2[], double t0, double tf, double h, double y1_init, double y2_init);
