/*-------------------------------------------------------------------------------\
@ Numerical Methods by Young-Keun Kim - Handong Global University

Author           : [An Gyeonheal]
ID               : [21900416]
Created          : 02-11-2022
Modified         : 24-11-2022
Language/ver     : C++ in MSVS2019

Description      : myNp.cpp
-------------------------------------------------------------------------------*/

#include "myNP.h"

/*===================================TEST1================================================*/

// QUESTION 1
double q1_integral(double func(const double x), double a, double b, int n) {
	double h = (b - a) / n;
	double I_integral = 0;

	for (int i = 0; i < n / 2; i++) {
		I_integral += 4 * func(a + h * (2 * i + 1)) + 2 * func(a + h * (2 * i + 2));
	}
	I_integral = (h / 3) * (I_integral - func(a) + func(b));

	return I_integral;
}
// QUESTION 2
double diffrentiation(double x, double q2_func(double _x))
{
	double df = 0;
	double h = 0.02;

	df = (q2_func(x + h) - q2_func(x - h)) / (2 * h);

	return df;
}
// QUESTION 3
double steffensen(double q3_func(double x), double x0, double tol)
{
	double N_max = pow(10, 100);
	double xn = x0;
	double ep = 0;
	double h = 0;
	int k = 0;

	do {

		h = -pow(q3_func(xn), 2) / (q3_func(xn + q3_func(xn)) - q3_func(xn));

		xn = xn + h;

		ep = fabs(q3_func(xn));

		k++;

	} while (k<N_max && ep>tol);

	return xn;
}

/*===================================TEST1================================================*/













/*-----------------------------Taylor Series------------------------------*/
double sinTaylor(double _x)
{
	int N_max = 10;
	double S_N = 0;

	for (int k = 0; k < N_max; k++)
		S_N = S_N + pow(-1, k) * pow(_x, 2 * k + 1) / factorial(2 * k + 1);

	return S_N;
}

double sindTaylor(double _x)
{
	return sinTaylor(_x * PI / 180);
}

double cosTaylor(double _x) {
	int N_max = 10;
	double C_N = 0;

	for (int k = 0; k < N_max; k++)
		C_N += pow(-1, k) * pow(_x, 2 * k) / factorial(2 * k);

	return C_N;
}

double factorial(int N)
{
	int y = 1;
	for (int k = 2; k <= N; k++)
		y = y * k;

	return y;
}

double power(double _x, int N)
{
	double y = 1;
	for (int k = 1; k <= N; k++) {
		y = y * _x;
	}
	return y;
}

double sinTaylor2(double _x)
{
	int N_max = 10;
	double S2_N = 0;

	for (int k = 0; k < N_max; k++)
		S2_N = S2_N + power(-1, k) * power(_x, 2 * k + 1) / factorial(2 * k + 1);

	return S2_N;
}

/*-----------------------------Bisection and TU_newtonRapson------------------------------*/
double bisection(double func(double x), float _a0, float _b0, float _tol)
{
	// Initialization
	int k = 0;
	int Nmax = 100;
	float a = _a0;
	float b = _b0;
	float xn = 0;
	float ep = 1000;

	// Bisection 
	while (k<Nmax && ep>_tol) {
		// Update xn as midpoint
		xn = (a + b) / 2;
		// Update range a, b
		if (((func(xn) * func(a)) < 0)) {
			b = xn;
		}
		else if (((func(xn) * func(a)) > 0)) {
			a = xn;
		}
		// Check tolerance
		ep = fabs(func(xn));

		k++;

		printf("k:%d \t", k);
		printf("Xn(k): %f \t", xn);
		printf("Tol: %.8f\n", ep);
	}

	return xn;
}

double TU_newtonRaphson(double func(double x), double dfunc(double x), double _x0, double _tol)
{
	double xn = _x0;
	double ep = 1000;
	int Nmax = 1000;
	int k = 0;
	double h = 0;

	do {
		if (func(xn) == 0)
		{
			printf("[ERROR] d2F == 0 !!\n");
			break;
		}
		else
		{
			// get h=f/df @ x(k)
			h = -func(xn) / dfunc(xn);
			printf("%f\t%f\t", func(xn), dfunc(xn));

			// update x(k+1)=x(k)+h(k)
			xn = xn + h;

			// check tolerance
			ep = fabs(func(xn));

			k++;

			printf("k:%d \t", k);
			printf("X(k): %f \t", xn);
			printf("Tol: %.10f\n", ep);

		}
	} while (k < Nmax && ep > _tol);

	return xn;
}



/*-----------------------------Nonlinear.cpp------------------------------*/
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



/*-----------------------------Diffrential.cpp------------------------------*/
// PROBLEM 1
void gradient1D(double x[], double y[], double dydx[], int m) {
	// Check if length x and y are equal	
	if (sizeof(x) != sizeof(y)) {
		printf("ERROR: length of x and y must be equal\n");
		return;
	}
	// Calculate h
	double h = x[1] - x[0];
	// Three-Point FWD  O(h^2). Need to use first 2 points
	dydx[0] = (-3 * y[0] + 4 * y[1] - y[2]) / (2 * h);
	// For mid points	
	for (int i = 1; i < m - 1; i++) {
		dydx[i] = (y[i + 1] - y[i - 1]) / (2 * h);
	}
	// For end point	
	// Three-Point BWD  O(h^2). Need to use last 2 points
	dydx[m - 1] = (y[m - 3] - 4 * y[m - 2] + 3 * y[m - 1]) / (2 * h);
}
// PROBLEM 2
void gradientFunc(double func(const double x), double x[], double dydx[], int m) {
	double* y;

	y = (double*)malloc(sizeof(double) * m);
	for (int i = 0; i < m; i++) {
		y[i] = func(x[i]);
	}
	gradient1D(x, y, dydx, m);

	free(y);
}
// PROBLEM 3
void acceleration(double x[], double y[], double dy2dx[], int m) {
	// Check if length x and y are equal	
	if (sizeof(x) != sizeof(y)) {
		printf("ERROR: length of x and y must be equal\n");
		return;
	}
	// calculate h
	double h = x[1] - x[0];
	// Four-point forward diffrence
	dy2dx[0] = (2 * y[0] - 5 * y[1] + 4 * y[2] - y[3]) / (h * h);
	// printf("%f\n", dy2dx[0]);
	// For mid points
	for (int i = 1; i < m - 1; i++) {
		dy2dx[i] = (y[i - 1] - 2 * y[i] + y[i + 1]) / (h * h);
	}
	// Four-point backward diffrence
	dy2dx[m - 1] = (-1 * y[m - 4] + 4 * y[m - 3] - 5 * y[m - 2] + 2 * y[m - 1]) / (h * h);
}



/*-----------------------------Integration.cpp------------------------------*/
// Trapezoidal method
double trapz(double x[], double y[], int m) {
	double I_trapz = 0;
	int N = m - 1;
	for (int i = 0; i < N; i++)
		I_trapz += ((x[i + 1] - x[i]) / 2) * (y[i] + y[i + 1]);
	return I_trapz;
}
// Simpson 1/3 method
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
// myfunc using Simson 1/3 method 
double integral(double func(const double x), double a, double b, int n) {
	double h = (b - a) / n;
	double I_integral = 0;

	for (int i = 0; i < n / 2; i++) {
		I_integral += 4 * func(a + h * (2 * i + 1)) + 2 * func(a + h * (2 * i + 2));
	}
	I_integral = (h / 3) * (I_integral - func(a) + func(b));

	return I_integral;
}

void odeEU(double myfunc(const double t, const double y), double y[], double t0, double tf, double h) {

	double t;
	double yi = 0;
	double slope = 0;
	int k = 0;
	y[k] = yi;
	
	printf("y[%d] = %f\n", k, y[k]);

	for (t = t0; t <= tf; t+=h) {
		slope = myfunc(t, yi);
		yi = yi + slope * h;

		k++;
		y[k] = yi;
		printf("%f\n",y[k]);
	}
}
void odeEM(double myfunc(const double t, const double y), double y[], double t0, double tf, double h) {

	double yi = 0;
	double yiE = 0;
	double slope1 = 0;
	double slope2 = 0;
	int k = 0;
	y[k] = yi;

	printf("y[%d] = %f\n", k, y[k]);

	for (double t = t0; t <= tf; t += h) {
		slope1 = myfunc(t, yi);
		yiE = yiE + slope1 * h;
		slope2 = myfunc(t + h, yiE);
		yi = yi + 0.5 * (slope1 + slope2) * h;

		k++;
		y[k] = yi;
		printf("%f\n", y[k]);
	}
}
void odeRK2(double myfunc(const double t, const double y), double y[], double t0, double tf, double h, double y0) {
	int alpha = 1;
	int beta = alpha;
	double C2 = 1 / 2 * alpha;
	double C1 = 1 - C2;
	double K1 = 0;
	double K2 = 0;
	double yi = y0;
	double yiE = y0;
	int k = 0;
	y[k] = y0;

	printf("y[%d] = %f\n", k, y[k]);

	for (double t = t0; t <= tf; t += h) {
		K1 = myfunc(t, yi);
		yiE = yiE + beta * K1 * h;
		K2 = myfunc(t + alpha * h, yiE);

		yi = yi + 0.5 * (K1 + K2) * h;

		k++;
		y[k] = yi;
		printf("%f\n", y[k]);
	}
		
}
void odeRK3(double myfunc(const double t, const double y), double y[], double t0, double tf, double h, double y0) {
	int alpha2 = 1;
	int alpha3 = 1;
	double beta21 = 0.5;
	int beta31 = -1;
	int beta32 = 2;
	double C1 = 1 / 6;
	double C2 = 4 / 6;
	double C3 = 1 / 6;
	double C4 = 1 / 2;

	double K1 = 0;
	double K2 = 0;
	double K3 = 0;

	double yi = y0;
	int k = 0;
	y[k] = y0;

	printf("y[%d] = %f\n", k, y[k]);
	
	for (int t = t0; t < tf; t+=h){

		K1 = myfunc(t, yi);
		K2 = myfunc(t + alpha2 * h, yi + beta21 * K1 * h);
		K3 = myfunc(t + alpha3 * h, yi + beta31 * K1 * h);
		
		yi = yi + (K1 + 4 * K2 + K3) * h / 6;

		k++;
		y[k] = yi;
		printf("%f\n", y[k]);
	}
}

void ode(double myfunc(const double t, const double y), double y[], double t0, double tf, double h, int method) {
	int y0 = 0;

	if (method == EU)
	{
		printf("\n================================================\n");
		printf("                     odeEU                        ");
		printf("\n================================================\n");
		odeEU(myfunc, y, t0, tf, h);
	}
	else if (method == EM)
	{
		printf("\n================================================\n");
		printf("                     odeEM                        ");
		printf("\n================================================\n");
		odeEM(myfunc, y, t0, tf, h);
	}
	else if (method ==RK2)
	{
		printf("\n================================================\n");
		printf("                     odeRK2                       ");
		printf("\n================================================\n");
		odeRK2(myfunc, y, t0, tf, h, y0);
	}
	else if (method == RK3)
	{
		printf("\n================================================\n");
		printf("                     odeRK3                       ");
		printf("\n================================================\n");
		odeRK2(myfunc, y, t0, tf, h, y0);
	}
	else
		printf("Error: Wrong method\n");
		system("pause");
}

void sys2RK2(void mckfunc(const double t, const double Y[], double dYdt[]), double y1[], double y2[], double t0, double tf, double h, double y1_init,
	double y2_init) {
	double K1[2] = { 0 };
	double K2[2] = { 0 };
	double tempY[2] = { 0 };
	double tempZ[2] = { 0 };
	double N = (tf - t0) / h + 1;
	double t = t0;
	double K1_0 = 0;							// K1 = { K1_0 ; K1_1 } 
	double K1_1 = 0;
	double K2_0 = 0;							// K2 = { K2_0 ; K2_0 } 
	double K2_1 = 0;

	y1[0] = y1_init;                          // y(t)					Y = { y(t) ; z(t) }
	y2[0] = y2_init;						  // z(t) 	

	for (int i = 0; i < N - 1; i++) {

		// slope 1 :							 f(t[i], y[i], z[i]) 
		tempY[0] = y1[i];						// Y[0] = y(t)
		tempY[1] = y2[i];						// Y[1] = z(t)

		mckfunc(t, tempY, K1);					// *** K1 = f(t[i], y[i], z[i]) ***

		// saving the calculated value (K1[] = dydt[])
		K1_0 = K1[0];							// K1_0 = dydt = z(t)
		K1_1 = K1[1];							// K1_1 = z'(t)

		// slope 2 							

		tempY[0] += K1_0 * h;					// y(t)
		tempY[1] += K1_1 * h;					// z(t)

		mckfunc(t + h, tempY, K2);				// *** K2 = f(t[i] + h, y[i] + K1 * h, z[i] + K1 * h) *** 

		// saving the calculated value (K2[] = dydt[])
		K2_0 = K2[0];							// z(t)
		K2_1 = K2[1];							// z'(t)

		y1[i + 1] = y1[i] + (0.5 * K1_0 + 0.5 * K2_0) * h;				// y(t)
		y2[i + 1] = y2[i] + (0.5 * K1_1 + 0.5 * K2_1) * h;				// z(t)

		t = t + h;								// t[i+1] = t[i] + h
	}

}
