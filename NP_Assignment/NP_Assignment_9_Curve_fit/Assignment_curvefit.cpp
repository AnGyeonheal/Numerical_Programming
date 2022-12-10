/*-------------------------------------------------------------------------------\
@ Numerical Methods by Young-Keun Kim - Handong Global University

Author           : [An Gyeonheal]
ID               : [21900416]
Created          : 27-11-2022
Modified         : 27-11-2022
Language/ver     : C++ in MSVS2019

Description      : Assignment_curvefit.cpp
-------------------------------------------------------------------------------*/

#include "../../../Include/myNP.h"


int main(int argc, char* argv[]) {

	printf("===================PROBLEM 1===================\n\n");
	printf("	    < Linear_Regression	> \n\n");

	double xi[] = { 30, 40, 50, 60, 70, 80 };
	double yi[] = { 1.05, 1.07, 1.09, 1.14, 1.17, 1.21 };
	double z[2] = { 0 };
	int m = sizeof(xi) / sizeof(xi[0]);
	
	linearRegression(xi, yi, m, z);

	printf("a1 = z[0] = %f \na0 = z[1] = %f \n\nP = a0 * T + a1\n\n", z[0], z[1]);		// a1 = z[0], a0 = z[1]
	double T = 100;
	double Pr = 0;
	Pr = z[1] * T + z[0];							// a0*T + a1
	printf("Pressure  = %f\n\n", Pr);

	printf("===================PROBLEM 2===================\n\n");
	printf("		 < Polyfit >				   \n\n");

	double yt[] = { 0, 3, 4.5, 5.8, 5.9, 5.8, 6.2, 7.4, 9.6, 15.6, 20.7, 26.7, 31.1, 35.6, 39.3, 41.5 };
	int N = sizeof(yt) / sizeof(double);
	int n = 4;

	Matrix Y = createMat(1,N);
	Matrix X = zeros(1, N);
	Matrix Z = zeros(n + 1, 1);
	for (int i = 0; i < N; i++) {
		X.at[0][i] = 0.4 * i;
		Y.at[0][i] = yt[i];
	}

	polyfit(X, Y, Z, n);

	printf("polynomial : [ %d ]\n\n", n);
	printMat(Z, "z");

	system("pause");
	return 0;
}

