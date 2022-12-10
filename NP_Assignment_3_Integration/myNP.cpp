/*-------------------------------------------------------------------------------\
@ Numerical Methods by Young-Keun Kim - Handong Global University

Author           : [An Gyeonheal]
ID               : [21900416]
Created          : 02-11-2022
Modified         : 02-11-2022
Language/ver     : C++ in MSVS2019

Description      : myNp.cpp
-------------------------------------------------------------------------------*/


#include "myNP.h"

double trapz(double x[], double y[], int m) {
	double I_trapz = 0;
	int N = m - 1;
	for (int i = 0; i < N; i++)
		I_trapz += ((x[i + 1] - x[i]) / 2) * (y[i] + y[i + 1]);
	return I_trapz;
}

double simpson13(double x[], double y[], int m) {
	double I_simpson13 = 0;
	int N = m - 1;
	double h = x[1] - x[0];

	for (int i = 0; i < N / 2; i++) {
		I_simpson13 = I_simpson13 + 4 * y[N - (2 * i + 1)] + 2 * y[N - (2 * i + 2)];
	}
	I_simpson13 = (h / 3) * (I_simpson13 - y[0] + y[N]);

	return I_simpson13;
}

double integral(double func(const double x), double a, double b, int n) {
	double h = (b - a) / n;
	double I_integral = 0;

	for (int i = 0; i < n / 2; i++) {
		I_integral += 4 * func(a + h * (2 * i + 1)) + 2 * func(a + h * (2 * i + 2));
	}
	I_integral = (h / 3) * (I_integral - func(a) + func(b));

	return I_integral;
}
