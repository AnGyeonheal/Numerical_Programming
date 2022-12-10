/*-------------------------------------------------------------------------------\
@ Numerical Methods by Young-Keun Kim - Handong Global University

Author           : [An Gyeonheal]
ID               : [21900416]
Created          : 02-11-2022
Modified         : 22-11-2022
Language/ver     : C++ in MSVS2019

Description      : Assignment_ode1.cpp
-------------------------------------------------------------------------------*/

#include "../../../Include/myNP.h"

double myfunc(const double t, const double y);
int main(int argc, char* argv[]) {

	double t0 = 0;
	double tf = 0.1;
	double h = 0.001;
	double y[101] = { 0 };

	ode(myfunc, y, t0, tf, h, EU);
	ode(myfunc, y, t0, tf, h, EM);
	ode(myfunc, y, t0, tf, h, RK2);
	ode(myfunc, y, t0, tf, h, RK3);

	system("pause");
	return 0;

}

double myfunc(const double t, const double y) {

	double tau = 1;
	double T = 1 / tau;
	double f = 10;
	double Vm = 1;
	double w = 2 * PI * f;
	
	double F = -T * y + T * Vm * cos(w * t);

	return F;
}



