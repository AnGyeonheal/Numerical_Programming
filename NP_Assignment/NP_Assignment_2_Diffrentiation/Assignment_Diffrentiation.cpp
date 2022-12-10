/*-------------------------------------------------------------------------------\
@ Numerical Methods by Young-Keun Kim - Handong Global University
Language/ver     : C++ in MSVS2019
Description      : Assignment for Differentiation
-------------------------------------------------------------------------------*/

/*==========================================================================*/
/*               Assignment   -    Differentiation	                            */
/*==========================================================================*/

#include "stdio.h"
#include "stdlib.h"
#include "\Users\82107\source\repos\NP\Include\myNP.h"

double myFunc(const double x);
void printVec(double* _vec, int _row);

int main(int argc, char* argv[])
{
	/************			PROBLEM 1						     ************/

	printf("\n***********************************************************************************************");
	printf("\n|  1. Create a function for numerical 1st order differentiation from a set of discrete data.  |");
	printf("\n***********************************************************************************************\n");

	int k = 12;
	double t1[12] = { 0 };
	for (int i = 0; i < k; i++) t1[i] = 0.5 * i;
	double x[] = { -3.632, -0.3935, 1, 0.6487, -1.282, -4.518, -8.611, -12.82, -15.91, -15.88, -9.402, 9.017};
	double  dxdt[12] = { 0 };

	gradient1D(t1, x, dxdt, k);

	printVec(dxdt, k);

	/************			PROBLEM 2						     ************/

	printf("\n****************************************************************************************************");
	printf("\n|   2. Create a function for numerical 1st order differentiation from user defined math equation   |");
	printf("\n****************************************************************************************************\n");

	int m = 21;
	double t2[21] = {0};
	for (int i = 0; i < m; i++) t2[i] = 0.2 * i;

	double dydx[21] = { 0 };
	gradientFunc(myFunc, t2, dydx, m);
	printVec(dydx, m);

	/************			PROBLEM 3						     ************/

	printf("\n***********************************************************");
	printf("\n|   3. Create a function for 2nd order differentiation    |");
	printf("\n***********************************************************\n");

	int n = 21;
	double t3[21] = { 0 };
	double y1[21] = { 0 };
	for (int i = 0; i < n; i++) {
		t3[i] = 0.2 * i;
		y1[i] = myFunc(i*0.2);
	}
	double dy2dx[21] = { 0 };
	acceleration(t3, y1, dy2dx, n);
	printVec(dy2dx, n);

	system("pause");
	return 0;
}
/*==========================================================================*/
/*						Function Definitions								*/
/*==========================================================================*/

// User defined function:  example  y=x^3
double myFunc(const double x) {
	return (x * x) * x;
}

void printVec(double* _vec, int _row)
{
	for (int i = 0; i < _row; i++)
		printf("Vector[%d] = %f \n", i, _vec[i]);
	printf("\n");
}