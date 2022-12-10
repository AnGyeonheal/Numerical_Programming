// 21900416 ¾È°ßÈú 2022-09-25

#include "myNP.h"

// Nonlinear Solver
double newtonRaphson(double dfunc(double x), double d2func(double x), double _x0, double _tol)
{
	double xn = _x0;
	double ep = 1000;
	int Nmax = 1000;
	int k = 0;
	double h = 0;

	do {
		if (d2func(xn) == 0)
		{
			printf("[ERROR] d2F == 0 !!\n");
			break;
		}
		else
		{
			// get h=f/df @ x(k)
			h = -dfunc(xn) / d2func(xn);
			printf("%f\t%f\t", dfunc(xn), d2func(xn));

			// update x(k+1)=x(k)+h(k)
			xn = xn + h;

			// check tolerance
			ep = fabs(dfunc(xn));

			k++;

			printf("k:%d \t", k);
			printf("X(k): %f \t", xn);
			printf("Tol: %.10f\n", ep);

		}
	} while (k < Nmax && ep > _tol);

	return xn;
}
