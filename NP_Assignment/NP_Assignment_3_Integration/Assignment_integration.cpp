/*-------------------------------------------------------------------------------\
@ Numerical Methods by Young-Keun Kim - Handong Global University

Author           : [An Gyeonheal]
ID               : [21900416]
Created          : 02-11-2022
Modified         : 02-11-2022
Language/ver     : C++ in MSVS2019

Description      : Assignment_Integration.cpp
-------------------------------------------------------------------------------*/

#include "../../../Include/myNP.h"

// Integration using rectangular method for discrete data inputs
double IntegrateRect(double x[], double y[], int m);
double myFunc(const double x);

int main(int argc, char* argv[])
{
	// PART 1. Integration from Datasets, may not evenly distribted
	printf("\n**************************************************");
	printf("\n        Problem 1. Integration from Datasets         ");
	printf("\n**************************************************\n");

	double x[] = { 0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60 };
	double y[] = { 0, 3, 8, 20, 33, 42, 40, 48, 60, 12, 8, 4, 3 };
	int M = sizeof(x) / sizeof(x[0]);

	double I_rect = IntegrateRect(x, y, M);
	printf("I_rect  = %f\n", I_rect);
	
	double I_trapz = trapz(x, y, M);
	printf("I_trapz = %f\n\n", I_trapz);


	// PART 2. Integration from Datasets with evenly distributed
	printf("\n**************************************************");
	printf("\n        Problem 2. Integration from Datasets, with the same interval       ");
	printf("\n**************************************************\n");

	double x2[] = { -3, -2.25, -1.5, -0.75, 0, 0.75, 1.5, 2.25, 3};
	double y2[] = { 0, 2.1875, 3.75, 4.6875, 5, 4.6875, 3.75, 2.1875, 0};
	int M2 = sizeof(x2) / sizeof(x2[0]);

	double I_simpson13 = simpson13(x2, y2, M2);
	printf("I_simpson13  = %f\n\n", I_simpson13);


	// PART 3. Integration from a Function
	printf("\n**************************************************");
	printf("\n        Problem 3. Integration from a Function       ");
	printf("\n**************************************************\n");

	double I_integral = integral(myFunc,-1, 1, 12);
	printf("I_integral  = %f\n\n", I_integral);


	system("pause");
	return 0;
}

// Integration using rectangular method for discrete data inputs
double IntegrateRect(double x[], double y[], int m) {
	int N = m - 1;
	double I = 0;
	for (int i = 0; i < N; i++)
		I += y[i] * (x[i + 1] - x[i]);

	return I;
}

double myFunc(const double x) {
	double F = 0;
	F = sqrt(1 - x * x);
	return F;
}
