/*-------------------------------------------------------------------------------\
@ Numerical Methods by Young-Keun Kim - Handong Global University

Author           : [An Gyeonheal]
ID               : [21900416]
Created          : 02-11-2022
Modified         : 24-11-2022
Language/ver     : C++ in MSVS2019

Description      : Assignment_ode2.cpp
-------------------------------------------------------------------------------*/

#include "../../../Include/myNP.h"


void mckfunc(const double t, const double Y[], double dYdt[]);

int main(int argc, char* argv[]) {

	int t0 = 0;
	double tf = 1;
	double h = 0.01;
	double y1[200] = { 0 };
	double y2[200] = { 0 };
	double y1_init = 0;
	double y2_init = 0.2;
	double N = (tf - t0) / h + 1;
	double t = 0;

	sys2RK2(mckfunc, y1, y2, t0, tf, h, y1_init, y2_init);

	for (int k = 0; k < N; k++) {
		/*printf("t = %f \n y(t) = %f \n z(t) = %f \n\n", t + k * h, y1[k], y2[k]);*/
		printf("%f\n", y1[k]);
	}
	system("pause");
	return 0;
}


void mckfunc(const double t, const double Y[], double dYdt[]) {
	double m = 1;
	double c = 7;
	double k = 6.9;
	double f = 5;
	double A = 2;

	double Fin = A * cos(2 * PI * f * t);
	/*double Fin = A;*/
	/*double Fin = 0;		*/					
																		// Y[0]  = y(t) , Y[1] = z(t)
	dYdt[0] = Y[1];														// dydt == z(t)
	dYdt[1] = (-k * Y[0] - c * Y[1] + Fin) / m;							// dzdt
}