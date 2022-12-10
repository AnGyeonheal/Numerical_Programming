/*-------------------------------------------------------------------------------\
@ Numerical Methods by Young-Keun Kim - Handong Global University

Author           : [An Gyeonheal]
ID               : [21900416]
Created          : 02-11-2022
Modified         : 27-11-2022
Language/ver     : C++ in MSVS2019

Description      : myNp.cpp
-------------------------------------------------------------------------------*/

#include "myNP.h"
/*===================================FUNCTION=============================================*/

// Taylorseries
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

// Nonlinear
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

// Diffrentiation
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
	dy2dx[0] = (2 * y[0] - 5 * y[1] + 4 * y[2] - y[3]) / (h * h);
	// printf("%f\n", dy2dx[0]);
	// For mid points
	for (int i = 1; i < m - 1; i++) {
		dy2dx[i] = (y[i - 1] - 2 * y[i] + y[i + 1]) / (h * h);
	}
	// Four-point backward diffrence
	dy2dx[m - 1] = (-1 * y[m - 4] + 4 * y[m - 3] - 5 * y[m - 2] + 2 * y[m - 1]) / (h * h);
}

// Integration
double trapz(double x[], double y[], int m) {
	double I_trapz = 0;
	int N = m - 1;
	for (int i = 0; i < N; i++)
		I_trapz += ((x[i + 1] - x[i]) / 2) * (y[i] + y[i + 1]);
	return I_trapz;
}
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
double integral(double func(const double x), double a, double b, int n) {
	double h = (b - a) / n;
	double I_integral = 0;

	for (int i = 0; i < n / 2; i++) {
		I_integral += 4 * func(a + h * (2 * i + 1)) + 2 * func(a + h * (2 * i + 2));
	}
	I_integral = (h / 3) * (I_integral - func(a) + func(b));

	return I_integral;
}

// ODE solver
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
	
	for (double t = t0; t < tf; t+=h){

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
void sys2RK2(void mckfunc(const double t, const double Y[], double dYdt[]), double y1[], double y2[], double t0, double tf, double h, double y1_init, double y2_init) {
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

// Curve fit
double linearRegression(double x[], double y[], int m, double z[]) {
	double sx = 0;
	double sxx = 0;
	double sxy = 0;
	double sy = 0;
	double S = 0;

	if (sizeof(x) != sizeof(y))
	{
		printf("ERROR: size of x and y is not same");
		system("pause");
		return 0;
	}
	for (int i = 0; i < m; i++) {
		sx += x[i];
		sxx += x[i] * x[i];
		sxy += x[i] * y[i];
		sy += y[i];
	}
	S = 1 / (m * sxx - sx * sx);
	z[0] = (sxx * sy - sxy * sx) * S;			// = a0
	z[1] = (m * sxy - sx * sy) * S;				// = a1

	return 0;
}
void polyfit(Matrix x, Matrix y, Matrix z, int n) {
	Matrix SX = zeros(1, 2 * n + 1);				// SX(1개, n+1개)
	Matrix S = zeros(n + 1, n + 1);					//  S(n+1개 , n+1개)
	Matrix b = zeros(n + 1, 1);
	Matrix invS = zeros(n + 1, n + 1);
	int m = x.cols;
	int cnt = 0;

	if (x.cols != y.cols) {
		printf("ERROR : The number of elements in x must be the same as in y\n");
		system("pause");
	}
	SX.at[0][0] = m;
	for (int l = 1; l <= n; l++) {
		double SXtemp1 = 0;
		double SXtemp2 = 0;
		double SXYtemp = 0;
		double SXYtemp0 = 0;
		for (int k = 0; k < m; k++) {
			SXtemp1 += pow(x.at[0][k], l);
			SXtemp2 += pow(x.at[0][k], 2 * n - l + 1);
			SXYtemp += y.at[0][k] * pow(x.at[0][k], n + 1 - l);
			SXYtemp0 += y.at[0][k] * pow(x.at[0][k], 0);
			cnt++;
		}
		b.at[l - 1][0] = SXYtemp;
		b.at[n][0] = SXYtemp0;
		SX.at[0][l] = SXtemp1;
		SX.at[0][2 * n - l + 1] = SXtemp2;
	}
	for (int i = 0; i < (double)(n + 1) / 2; i++) {
		S.at[i][i] = SX.at[0][2 * n - 2 * i];
		S.at[n - i][n - i] = SX.at[0][2 * i];
		S.at[i][n - i] = SX.at[0][n];
		S.at[n - i][i] = SX.at[0][n];
		for (int j = 0; j < i; j++) {
			S.at[i][j] = SX.at[0][n - i + n - j];
			S.at[j][i] = SX.at[0][n - i + n - j];		
			S.at[n - i][n - j] = SX.at[0][i + j];		
			S.at[n - j][n - i] = SX.at[0][i + j]; 
			S.at[j][n - i] = SX.at[0][n + i - j]; 
			S.at[n - i][j] = SX.at[0][n + i - j]; 
			S.at[i][n - j] = SX.at[0][n - i + j];
			S.at[n - j][i] = SX.at[0][n - i + j];
			cnt++;
		}
	}
	invMat(S, invS);
	copyVal(multiply(invS, b), z);
	printf("iteration = %d\n\n", cnt);
}

/*===================================MATRIX FUNCTION=============================================*/

Matrix	addMat(Matrix _A, Matrix _B)
{
	if (_A.rows != _B.rows || _A.cols != _B.cols) {
		printf("\n*************************************************");
		printf("\n  ERROR!!: dimension error at 'addMat' function");
		printf("\n*************************************************\n");
		return createMat(0, 0);
	}

	Matrix Out = createMat(_A.rows, _B.cols);
	for (int i = 0; i < _A.rows; i++)
		for (int j = 0; j < _B.cols; j++)
			Out.at[i][j] = _A.at[i][j] + _B.at[i][j];

	return Out;
}
Matrix	backSub(Matrix U, Matrix y, Matrix x)
{
	int i, j;

	for (i = U.rows - 1; i >= 0; i--)
	{
		double sub = 0;
		for (j = i + 1; j < U.rows; j++)
		{
			sub += U.at[i][j] * x.at[j][0];
		}
		x.at[i][0] = (y.at[i][0] - sub) / U.at[i][i];
	}
	return x;
}
Matrix fwdSub(Matrix L, Matrix b, Matrix y)
{
	int i, j;

	for (i = 0; i < L.rows; i++) {
		double sub = 0;
		for (j = 0; j < L.rows - 1; j++) {
			sub += L.at[i][j] * y.at[j][0];
		}
		y.at[i][0] = (b.at[i][0] - sub) / L.at[i][i];
	}
	return y;
}
double invMat(Matrix A, Matrix Ainv) {

	if (A.rows != A.cols) {
		printf("ERROR: IT IS NOT A N*N MATRIX");
	}
	for (int i = 0; i < A.cols; i++) {
		if (A.at[i][i] == 0) {
			printf("ERROR: IT IS NOT A FULL RANK MATRIX");

			return 0;
		}
	}
	Matrix L = zeros(A.rows, A.cols);
	Matrix U = zeros(A.rows, A.cols);
	Matrix P = eye(A.rows, A.cols);
	Matrix I = eye(A.rows, A.cols);
	Matrix x = zeros(A.rows, 1);
	Matrix AI = zeros(A.rows, 1);
	Matrix Uc = zeros(A.rows, A.cols);
	Matrix Lc = zeros(A.rows, A.cols);
	Matrix result = zeros(A.rows, A.cols);

	LUdecomp(A, L, U, P);

	for (int j = 0; j < A.cols; j++) {
		for (int i = 0; i < A.cols; i++) {
			AI.at[i][0] = I.at[i][j];
		}															// I를 열 벡터로 만들기
		x = backSub(U, AI, x);										// U * x = AI 를 back substitution으로 x 벡터 구하기
		for (int i = 0; i < A.cols; i++) {
			Uc.at[i][j] = x.at[i][0];
		}initMat(x, 0);												// 각 x벡터들을 Uc에 할당
	}
	for (int j = 0; j < A.cols; j++) {
		for (int i = 0; i < A.cols; i++) {
			AI.at[i][0] = I.at[i][j];
		}															// I를 열 벡터로 만들기
		x = fwdSub(L, AI, x);										// L * x = AI 를 back substitution으로 x 벡터 구하기
		for (int i = 0; i < A.cols; i++) {
			Lc.at[i][j] = x.at[i][0];
		}initMat(x, 0);												// 각 x벡터들을 Lc에 할당
	}
	copyVal(multiply(Uc, Lc), Ainv);								// copyVal을 통해 Ainv에 저장, Uc * Lc
	return 0;
}
// Gauss-Elimination with pivoting 
void gaussElim(Matrix A, Matrix b, Matrix U, Matrix d, Matrix P)
{
	if (A.rows != A.cols || b.rows != A.rows || A.rows != U.rows || A.cols != U.cols || A.rows != P.rows || A.cols != P.cols || d.rows != U.rows || b.rows != d.rows || b.cols != 1 || d.cols != 1) {
		printf("Square Matrix is not Used || Vector's rows are not same with Matrix's rows || Vector's cols are not 1");
		return;
	}															// 예외 처리 (행렬, 벡터 에러)

	int i, j, k;
	double max_a;
	double temp;
	double m;

	copyVal(eye(P.rows, P.cols), P);							// Permutation matrix
	copyVal(A, U);
	copyVal(b, d);

	for (k = 0; k < U.rows - 1; k++)
	{
		max_a = fabs(U.at[k][k]);
		for (i = k + 1; i < U.rows; i++)
		{
			if (max_a < fabs(U.at[i][k]))
			{
				max_a = fabs(U.at[i][k]);						// 가장 큰 행렬 찾기
				for (j = 0; j < U.rows; j++) {
					temp = U.at[i][j];
					U.at[i][j] = U.at[k][j];
					U.at[k][j] = temp;

					temp = P.at[i][j];
					P.at[i][j] = P.at[k][j];
					P.at[k][j] = temp;
				}
			}
		}														// 행렬 재정렬
		for (i = k + 1; i < U.rows; i++) {
			if (U.at[k][k] == 0) {
				printf("There's 0 in pivot element");			// 예외 처리 (0으로 나눌 때)
				return;
			}
			m = U.at[i][k] / U.at[k][k];
			for (j = k; j < U.cols; j++) {
				U.at[i][j] = U.at[i][j] - (m * U.at[k][j]);
			}
			d.at[i][0] = d.at[i][0] - m * d.at[k][0];

		}														// Gauss-Elimination
	}
}
// LU-Decomposition with pivoting
void LUdecomp(Matrix A, Matrix L, Matrix U, Matrix P) {

	Matrix sp = zeros(1, A.cols);
	double temp = 0;
	double max_sp = 0;
	double max_a = 0;
	int t = 0;
	double m = 0;
	copyVal(A, U);												// CopyVal을 쓰는 이유: 출력을 처음 입력과 같게 만든 후 출력 상에서 전개 후에 그대로 출력하기 위해서

	if (A.rows != A.cols) {
		printf("ERROR: IT IS NOT A N*N MATRIX");
		return;
	}

	for (int k = 0; k < U.rows - 1; k++) {
		/*for (int i = k; i < U.cols; i++) {
			for (int j = 0; j < U.cols; j++) {
				if (max_a < fabs(U.at[i][j])) {
					max_a = fabs(U.at[i][j]);					// 같은 행에서 max 값 찾기
				} 
			}
			sp.at[0][i] = fabs(U.at[i][k] / max_a);				// 각 행의 sp 값 구하기
			max_a = 0;											// max값 초기화
		}
		for (int i = 0; i < U.rows; i++){
			if (max_sp < sp.at[0][i]) {
				max_sp = sp.at[0][i];							// 각 행의 sp 값 중 가장 큰 sp 찾기
				t = i;											// 가장 큰 sp를 가진 행의 주소 값 찾기
			}
		}
		if (k != t) {											// 바꾸는 행과 바꿔지는 행의 sp값이 같을 때
			for (int i = 0; i < U.rows; i++) {
				temp = U.at[t][i];
				U.at[t][i] = U.at[k][i];
				U.at[k][i] = temp;								// 기존 행과 sp가 가장 큰 행 바꾸기

				temp = P.at[t][i];
				P.at[t][i] = P.at[k][i];
				P.at[k][i] = temp;								// P 바꾸기

				temp = L.at[t][i];
				L.at[t][i] = L.at[k][i];
				L.at[k][i] = temp;								// L 바꾸기
			}
		}*/
		for (int i = k + 1; i < U.rows; i++) {
			m = U.at[i][k] / U.at[k][k];
			L.at[i][k] = m;										// L 만들기
			for (int j = k; j < U.rows + 1; j++) {
				U.at[i][j] = U.at[i][j] - m * U.at[k][j];		// 가우스 소거법 --> U 만들기
			}
		}
		initMat(sp, 0);											// sp 행렬, sp 최대 값 초기화
		max_sp = 0;
	}
	for (int i = 0; i < L.rows; i++) {
		L.at[i][i] = 1;
	}
}
// Solving LU
void solveLU(Matrix L, Matrix U, Matrix P, Matrix b, Matrix x) {
	// PAx = Pb --> LUx = Pb
	Matrix y = zeros(b.rows, b.cols);							// LUx = b,   Ux = y
	Matrix ny = fwdSub(L, multiply(P, b), y);					// Ly=b (L,U 모두 scaled pivoting을 했으므로 vector b 또한 pivoting)
	x = backSub(U, ny, x);										// Ux=y

}
// eigen value and eigen vector
Matrix eig(Matrix A)
{
	if (A.rows != A.cols) {
		printf("Error: Please Input nXn Matrix\n");
		system("pause");
	}
	Matrix e = zeros(A.cols, 1);
	Matrix c = createMat(A.rows, 1);
	Matrix v = createMat(A.rows, 1);
	Matrix vt = createMat(1, A.cols);
	Matrix vvt = createMat(A.rows, A.cols);
	Matrix H = createMat(A.rows, A.cols);
	Matrix I = eye(H.rows, H.cols);
	Matrix R = createMat(A.rows, A.cols);
	Matrix eigenVal = zeros(1, A.cols);
	Matrix U = createMat(A.rows, A.cols);

	double normc = 0;
	double vtv = 0;
	int i, j, k, l, m;

	copyVal(A, R);
	for (m = 0; m < 10; m++) {
		Matrix Q = eye(A.rows, A.cols);										// Initializing Q = I
		// QR decomposition
		for (l = 0; l < A.rows - 1; l++) {
			for (k = 0; k < A.rows - 1; k++) {
				vtv = 0;													// initializing 
				initMat(e, 0);
				for (i = 0; i < A.cols; i++) {								// c 
					c.at[i][0] = R.at[i][k];
				}
				for (i = 0; i < k; i++) {
					c.at[i][0] = 0;
				}
				normc = volMat(c);											// ||c|| 

				if (c.at[k][0] >= 0)
					e.at[k][0] = 1;											// e
				else
					e.at[k][0] = -1;

				v = addMat(c, xMat(e, normc));								// v = c + |c|e
				vt = transpose(v);

				for (i = 0; i < A.rows; i++) {								// vtv = vt*v
					vtv += v.at[i][0] * vt.at[0][i];
					for (j = 0; j < A.cols; j++) {							// vvt = v*vt
						vvt.at[i][j] = v.at[i][0] * vt.at[0][j];
					}
				}
				if (vtv == 0) {												// 예외처리: 0으로 나눠질 때
					printf("Error: There is division by zero\n");
					system("pause");
				}
				copyVal(subMat(I, xMat(vvt, 2 / vtv)), H);					// H = I - 2*(v*vt)/(vt*v)
				copyVal(multiply(H, R), R); 								// R = H*R
				copyVal(multiply(Q, H), Q);									// Q = Q*H
			}
		}
		copyVal(multiply(R, Q), U);											// U = RQ
		copyVal(U, R);														// U-->R , QRdecomp(U)
	}
	printMat(R, "R");
	for (i = 0; i < A.rows; i++) {
		eigenVal.at[0][i] = R.at[A.rows - 1 - i][A.rows - 1 - i];				// eigenValue 
	}
	return eigenVal;
}
Matrix eigvec(Matrix A) {
	Matrix V = zeros(A.rows - 1, 1);
	Matrix I = eye(A.rows, A.cols);
	Matrix lamda = zeros(1, A.cols);
	Matrix B = zeros(A.rows, A.cols);
	Matrix subB = zeros(A.rows - 1, A.cols - 1);
	Matrix vecB = zeros(A.rows - 1, 1);
	Matrix eigenVec = eye(A.rows, A.cols);
	double normv = 0;
	int i, j, k;

	if (A.rows != 2 && A.rows != 3) {								// 예외처리:2by2, 3by3 행렬이 아닐 때
		printf("ERROR: Please Input 2X2 or 3X3 Matrix!\n");
		system("pause");
	}
	if (A.rows != A.cols) {											// 예외처리:n by n 행렬이 아닐 때
		printf("Error: Please Input nXn Matrix\n");
		system("pause");
	}

	for (k = 0; k < A.cols; k++) {
		lamda = eig(A);											// Initiallization
		Matrix invsubB = zeros(A.rows - 1, A.cols - 1);
		double x = 0;

		B = subMat(A, xMat(I, lamda.at[0][k]));
		if (k == 0) {											// V1
			for (i = 0; i < subB.rows; i++) {					// subB 1
				for (j = 0; j < subB.cols; j++) {
					subB.at[i][j] = B.at[i + 1][j + 1];
				}
			}
			for (i = 0; i < vecB.rows; i++) {
				vecB.at[i][0] = -B.at[i + 1][k];				// vecB 1
			}
			invMat(subB, invsubB);
			copyVal(multiply(invsubB, vecB), V);					// v1 = inv(subB)*vecB
			for (i = 0; i < V.rows; i++) {
				eigenVec.at[i + 1][k] = V.at[i][k];				// eigenVec에 v1 할당
			}
			for (i = 0; i < V.rows; i++) {
				x += pow(V.at[i][0], 2);
			}
			x += pow(1, 2);
			normv = sqrt(x);									// |v1| 구하기
			if (normv == 0) {
				printf("ERROR: There is division by 0");		// 예외처리: Division by 0
				system("pause");
			}
			for (i = 0; i < eigenVec.rows; i++) {
				eigenVec.at[i][k] = eigenVec.at[i][k] / normv;	// normalization
			}
		}
		else if (k > A.cols - A.cols && k < A.cols - 1) {		// V2
			subB.at[0][0] = B.at[0][0];							// subB 2
			subB.at[0][1] = B.at[0][2];
			subB.at[1][0] = B.at[2][0];
			subB.at[1][1] = B.at[2][2];
			vecB.at[0][0] = -B.at[0][1];						// vecB 2
			vecB.at[1][0] = -B.at[2][1];
			invMat(subB, invsubB);
			copyVal(multiply(invsubB, vecB), V);				// v2 = inv(subB)*vecB
			eigenVec.at[0][1] = V.at[0][0];
			eigenVec.at[2][1] = V.at[1][0];						// eigenVec에 V2 할당
			for (i = 0; i < V.rows; i++) {
				x += pow(V.at[i][0], 2);
			}
			x += pow(1, 2);										// |v2| 구하기
			normv = sqrt(x);
			if (normv == 0) {									// 예외처리: Division by 0
				printf("ERROR: There is division by 0");
				system("pause");
			}
			for (i = 0; i < eigenVec.rows; i++) {
				eigenVec.at[i][k] = eigenVec.at[i][k] / normv;	// normalization
			}
		}
		else if (k >= A.cols - 1) {								// V3
			for (i = 0; i < subB.rows; i++) {					// subB 3
				for (j = 0; j < subB.cols; j++) {
					subB.at[i][j] = B.at[i][j];
				}
			}
			for (i = 0; i < vecB.rows; i++) {
				vecB.at[i][0] = -B.at[i][k];					// vecB 3
			}
			invMat(subB, invsubB);
			copyVal(multiply(invsubB, vecB), V);				// v3 = inv(subB)*vecB
			for (i = 0; i < V.rows; i++) {
				eigenVec.at[i][k] = V.at[i][0];				    // eigenVec에 v3 할당
			}
			for (i = 0; i < V.rows; i++) {
				x += pow(V.at[i][0], 2);
			}
			x += pow(1, 2);
			normv = sqrt(x);									// |v3| 구하기
			if (normv == 0) {
				printf("ERROR: There is division by 0");		// 예외처리: Division by 0
				system("pause");
			}
			for (i = 0; i < eigenVec.rows; i++) {
				eigenVec.at[i][k] = eigenVec.at[i][k] / normv;	// normalization
			}
		}
	}
	return eigenVec;
}


/*===================================MID_TEST=============================================*/
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