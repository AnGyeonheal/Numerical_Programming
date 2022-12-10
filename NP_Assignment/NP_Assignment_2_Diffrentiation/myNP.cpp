// 21900416 ¾È°ßÈú 2022-09-26

#include "myNP.h"

// problem 1, problem 2, problem 3

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

// Truncation error should be O(h^2) 
void gradientFunc(double func(const double x), double x[], double dydx[], int m) {
	double* y;

	y = (double*)malloc(sizeof(double) * m);
	for (int i = 0; i < m; i++) {
		y[i] = func(x[i]);
	}
	gradient1D(x, y, dydx, m);

	free(y);
}

void acceleration(double x[], double y[], double dy2dx[], int m) {
	// Check if length x and y are equal	
	if (sizeof(x) != sizeof(y)) {
		printf("ERROR: length of x and y must be equal\n");
		return;
	}
	// calculate h
	double h = x[1] - x[0];
	// Four-point forward diffrence
	dy2dx[0] = (2 * dy2dx[0] - 5 * dy2dx[1] + 4 * dy2dx[2] - dy2dx[3]) / (h * h);
	// For mid points
	for (int i = 1; i < m - 1; i++) {
		dy2dx[i] = (y[i - 1] - 2 * y[i] + y[i + 1]) / (h * h);
	}
	// Four-point backward diffrence
	dy2dx[m - 1] = (-1 * y[m - 4] + 4 * y[m - 3] - 5 * y[m - 2] + 2 * y[m - 1]) / (h * h);
}