/*------------------------------------------------------------------------------------------\
@ Numerical Programming by Young-Keun Kim - Handong Global University

Name	        : Your Name Goes Here
ID		        : Your ID  Goes Here
Created         : 07-03-2019

Language/ver    : C /  MSVS2017
Course          : Numerical Programming 2022
Description     : Assignment of Non-linear solver
/------------------------------------------------------------------------------------------*/



#include "stdio.h"
#include "stdlib.h"
#include <math.h>

#include "../../../Include/myNP.h"

double func(double x)
{
	double L = 4;
	double E = 70 * pow(10, 9);
	double I = 52.9e-6;
	double w0 = 20 * pow(10, 3);

	double F = 0;
	F = (w0 * x / (360 * L * E * I)) * (7 * pow(L, 4) - 10 * pow(L, 2) * pow(x, 2) + 3 * pow(x, 4));

	return F;
}


double dfunc(double x)
{
	double L = 4;
	double E = 70 * pow(10, 9);
	double I = 52.9 * pow(10, -6);
	double w0 = 20 * pow(10, 3);

	double dF = 0;
	dF = (w0 / (360 * L * E * I)) * (7 * pow(L, 4) - 30 * pow(L, 2) * pow(x, 2) + 15 * pow(x, 4));

	return dF;
}

double d2func(double x)
{
	double L = 4;
	double E = 70 * pow(10, 9);
	double I = 52.9 * pow(10, -6);
	double w0 = 20 * pow(10, 3);

	double d2F = 0;
	d2F = (w0 / (360 * L * E * I)) * (-60 * pow(L, 2) * x + 60 * pow(x, 3));

	return d2F;
}

void main() {

	float tol = 0.000001;

	/*==========================================================================*/
	/*               Assignment -     Newton Rhapson                            */
	/*==========================================================================*/

	/************      Variables declaration & initialization      ************/
	double sol_nr;
	double x0 = 3;

	printf("------------------------------------------------------------------------------------\n");
	printf("         Newton-Raphson Method Results             \n");
	printf("------------------------------------------------------------------------------------\n");

	/************      Solve  &		Show Output			           ************/
	printf("Newton-Raphson Method Result:\n");
	sol_nr = newtonRaphson(dfunc, d2func, x0, tol);

	printf("Final Solution: %f \t", sol_nr);
	printf("\n");


	system("pause");
}



