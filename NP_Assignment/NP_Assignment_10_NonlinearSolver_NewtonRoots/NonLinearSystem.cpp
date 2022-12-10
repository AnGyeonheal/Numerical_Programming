
/*------------------------------------------------------------------------------------------\
@ Numerical Methods by Young-Keun Kim - Handong Global University

Author          : Young-Keun Kim
Created         : 01-04-2019
Modified        : 12-08-2022
Language/ver	: C in MSVS2017
Course			: Numerical Programming

Description     : TU System of NonLinear
/------------------------------------------------------------------------------------------*/

#include "../../../Include/myNP.h"

void myjacobEX1(Matrix X, Matrix J);
void myfuncEX1(Matrix X, Matrix F);
void myjacobEX2(Matrix X, Matrix J);
void myfuncEX2(Matrix X, Matrix F);

int main(int argc, char* argv[]) {

	// problem 1

	printf("\n=========PROBLEM1===========\n");

	double xval[] = { 2.5, 2 };
	Matrix X1 = createMat(2, 1);
	X1 = arr2Mat(xval, 2, 1);
	double tol = 0.00001;
	
	newtonRoot(X1, myfuncEX1, myjacobEX1, tol);

	// problem 2
	
	printf("\n=========PROBLEM2===========\n");

	double x2val[] = { 10.0 / 180 * PI, 10, 10 };
	Matrix X2 = createMat(3, 1);
	X2 = arr2Mat(x2val, 3, 1);

	newtonRoot(X2, myfuncEX2, myjacobEX2, tol);


	system("pause");
	return 0;
}

void myjacobEX1(Matrix X, Matrix J) {
	// x[0] = x
	// x[1] = y
	J.at[0][0] = -0.25 *(exp(X.at[0][0] / 2) - exp(-X.at[0][0] / 2));
	J.at[0][1] = 1;
	J.at[1][0] = 18 * X.at[0][0];
	J.at[1][1] = 50 * X.at[1][0];
}

void myfuncEX1(Matrix X, Matrix F) {
	F.at[0][0] = X.at[1][0] - 0.5 * (exp(X.at[0][0] / 2) + exp(-X.at[0][0] / 2));
	F.at[1][0] = 9 * pow(X.at[0][0], 2) + 25 * pow(X.at[1][0], 2) - 225;
}

void myjacobEX2(Matrix X, Matrix J) {
	double x0 = 0;
	double x1 = 0;
	double y0 = 100;
	double y1 = -100;
	double th = 0;
	double dx = 0;
	double dy = 0;

	th = X.at[0][0];
	dx = X.at[1][0];
	dy = X.at[2][0];

	J.at[0][0] = -y0 * cos(th) - x0 * sin(th);
	J.at[0][1] = 1;
	J.at[0][2] = 0;
	J.at[1][0] = x0 * cos(th) - y0 * sin(th);
	J.at[1][1] = 0;
	J.at[1][2] = 1;
	J.at[2][0] = -y1 * cos(th) - x1 * sin(th);
	J.at[2][1] = 1;
	J.at[2][2] = 0;
}

void myfuncEX2(Matrix X, Matrix F) {
	double x0 = 0;
	double x1 = 0;
	double y0 = 100;
	double y1 = -100;
	double th = 0;
	double dx = 0;
	double dy = 0;
	double x0_new = 50;
	double y0_new = 186.6025;
	double x1_new = 150;
	double y1_new = 13.3975;

	th = X.at[0][0];
	dx = X.at[1][0];
	dy = X.at[2][0];

	F.at[0][0] = x0 * cos(th) - y0 * sin(th) + dx - x0_new;
	F.at[1][0] = x0 * sin(th) + y0 * cos(th) + dy - y0_new;
	F.at[2][0] = x1 * cos(th) - y1 * sin(th) + dx - x1_new;
}