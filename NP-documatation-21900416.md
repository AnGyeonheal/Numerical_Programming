# **<span style="color:blue">Numerical Programming_2022</span>**

# <span style="background-color:orange">.cpp Code of < System of NonLinear Equation ></span>

# <span style="background-color:yellow">myNP.cpp</span>

```c++
#include "myNP.h"
```

### sinTaylor()

Taylor series of sin function

```c++
double sinTaylor(double _x)
```

#### Parameters

+ _x : sin parameter (radian)

### sindTaylor()

Taylor series of sin function

```c++
double sindTaylor(double _x)
```

#### Parameters

+ _x : sin parameter (degree)

### cosTaylor()

Taylor series of cos function

```
double cosTaylor(double _x)
```

#### Parameters

+ _x : sin parameter (radian)

### factorial()

Function of getting factorial

```c++
double factorial(int N)
```

#### Parameters

+ N : number parameter

### power()

Function of getting squared number

```c++
double power(double _x, int N)
```

#### Parameters

+ N : number parameter (x^N)
+ _x : number parameter (x^N)

### sinTaylor2()

Taylor series using own function

```c++
double sinTaylor2(double _x)
```

#### Parameters

+ _x : number parameter

# bisection()

Bisection method to get root number

```c++
double bisection(double func(double x), float _a0, float _b0, float _tol)
```

**Parameters**

+ func(double x) : number parameter that in function
+ _a0 : initial number which is left side of root number
+ _b0 : initial number which is right side of root number
+ _tol : tolerance

# newtonRaphson()

NewtonRaphson method to get root number

```c++
double newtonRaphson(double dfunc(double x), double d2func(double x), double _x0, double _tol)
```

**Parameters**

+ dfunc(double x) : number parameter that in differential function
+ d2func(double x) : number parameter that in twice differential function
+ _x0 : root number
+ _tol : tolerance

# gradient1D()

Getting gradient of function that is series of number using Three-point FWD, mid-point and Three-point BWD.

```c++
void gradient1D(double x[], double y[], double dydx[], int m)
```

**Parameters**

+ x[] : number of X axis
+ y[] : number of Y axis
+ dydx[] : gradient number
+ m : times 

# gradientFunc()

Getting gradient from own function using gradient1D() function

```c++
void gradientFunc(double func(const double x), double x[], double dydx[], int m)
```

**Parameters**

+ func(const double x) : function number
+ double x[] : number of X axis
+ dydx[] : gradient number
+ m : times 

# acceleration()

Getting twice differentiated number.

```c++
void acceleration(double x[], double y[], double dy2dx[], int m)
```

**Parameters**

+ x[] : number of X axis
+ y[] : number of Y axis
+ dy2dx[] : twice differentiated number
+ m : times 

# trpaz()

Trapezoidal method for Integration

```c++
double trapz(double x[], double y[], int m)
```

**Parameters**

+ x[] : number of X axis
+ y[] : number of Y axis
+ m : times 

# simpson13()

Simpson method for Integration

```c++
double trapz(double x[], double y[], int m)
```

**Parameters**

+ x[] : number of X axis
+ y[] : number of Y axis
+ m : times 

# Integral()

Integral Using simpson13()

```c++
double integral(double func(const double x), double a, double b, int n)
```

**Parameters**

+ func(const double x) : function number(x)
+ a : initial number (getting start of Integration)
+ b : final number (closing Integration)
+ n : section



# <span style="background-color:orange">.cpp Code of < System of Linear Equation ></span>

# <span style="background-color:yellow">myMatrix.cpp</span>

```c++
#include "myMatrix.h"
```

# createMat()

Create Matrix

```c++
Matrix	createMat(int _rows, int _cols);
```

# arr2Mat()

Transfer array elements to Matrix

```c++
Matrix   arr2Mat(double* _1Darray, int _rows, int _cols)
```

**Parameters**

+ _1Darray : Array which want to move elements
+ _rows : Matrix 행의 개수
+ _cols : Matrix 열의 개수

example code :

```c++
	double a[] = { 1,2,3,4,5,6 };
	Matrix b = arr2Mat(a, 3, 2);
```

```
b =
       1.000000        2.000000
       3.000000        4.000000
       5.000000        6.000000
```

# freeMat()

Deallocate Memory

```c++
void	freeMat(Matrix _A);
```

**Parameters**

+ _A : Matrix form. (Should be nxn square.)

# txt2Mat()

Bring .txt File to Matrix

```c++
Matrix	txt2Mat(std::string _filePath, std::string _fileName);
```

# printMat()

Print the Matrix

```c++
void	printMat(Matrix _A, const char* _name);
```

**Parametersc**

+ _A : Matrix form.
+ char*_name: Explaining what this Matrix is.

# initMat()

Initialize Matrix elements

```c++
void	initMat(Matrix _A, double _val);
```

**Parameters**

+ _A : Matrix form. (Should be nxn square.)
+ _val : Initializing Matrix element with this value.

# zeros()

Create matrix of all zeros

```c++
Matrix	zeros(int _rows, int _cols);
```

**Parameters**

+ _rows : Set the number of rows
+ _cols : Set the number of cols

# eye()

Create identity Matrix

```c++
Matrix	eye(int _rows, int _cols);
```

**Parameters**

+ _rows : Set the number of rows
+ _cols : Set the number of cols

# transpose()

Transpose the Matrix

```c++
Matrix	transpose(Matrix _A);
```

**Parameters**

+ _A : Matrix form. (Should be nxn square.)

# copyMat()

Copy the Matrix

```c++
Matrix	copyMat(Matrix _A)
```

**Parameters**

+ _A : Matrix form. (Should be nxn square.)

# copyVal()

Copy the Matrix elements from _A to _B

```c++
void	copyVal(Matrix _A, Matrix _B);
```

**Parameters**

+ _A : Matrix which want to copy.  

+ _B : Matrix which want to paste.

  ** A  matrix's number of rows and columns should be same with _B*

# multiply()

Multiply Matrix _A and _B  

ex) P = AB 

```c++
Matrix multiply(Matrix _A, Matrix _B);
```

**Parameters**

+ _A : Matrix form. (Should be nxn Matrix.)

+ _B : Matrix form. (Should be nxn Matrix.)

  ** A  matrix's number of rows and columns should be same with _B*

# subMat()

Subtract Matrix _A - _B

ex) P = A - B

```c++
Matrix subMat(Matrix _A, Matrix _B);
```

**Parameters**

+ _A : Matrix form. (Should be nxn Matrix.)

+ _B : Matrix form. (Should be nxn Matrix.)

  ** A  matrix's number of rows and columns should be same with _B*

# xMat()

Multiply scalar with Matrix

ex) P = c * A

```c++
Matrix xMat(Matrix _A, double _c);
```

**Parameters**

+ _A : Matrix form.
+ _c : Scalar which want to multiply

# volMat()

Get volume of Matrix

```c++
double volMat(Matrix _A);
```

**Parameters**

+ _A : Matrix form.

# <span style="background-color:yellow">myNP.cpp</span>

```c++
#include "myNP.h"
```

# addMat()

Add Matrix _A and _B

ex) P = _A + _B

```c++
Matrix	addMat(Matrix _A, Matrix _B);
```

**Parameters**

+ _A : Matrix form.
+ _B: Matrix form.

# gaussElim()

Transform the linear system into an equivalent form with an upper and lower triangular matrix

Ay=b 인 시스템을 Ux = d 로 변환한다.

Gauss Elimination 과정에서 0으로 나눠지는 것을 방지하기 위해 Partial Pivoting 기법이 들어갔는데 이 때 P는 Pivoting된 위치이다.



```c++
void gaussElim(Matrix A, Matrix b, Matrix U, Matrix d, Matrix P);
```



**Parameters**

+ A : Input Matrix. (nxn)
+ b : Input Vector. (nx)
+ U : Output Matrix (upper triangular matrix) (nxn)
+ d : Output Vector (nx1)
+ P : Pivoting Matrix (Permutation Matrix)



사용 예제는 다음과 같다.

```c++
Matrix matA = txt2Mat(path, "matA");
Matrix vecb = txt2Mat(path, "vecb");
Matrix matU = zeros(matA.rows, matA.cols);
Matrix vecd = zeros(vecb.rows, vecb.cols);
Matrix P = eye(matA.rows, matA.cols);

gaussElim(matA, vecb, matU, vecd, P);
```
```
================= Input =================
matrix A =
       1.000000        3.000000       -2.000000        4.000000
       2.000000       -0.100000        3.000000       -1.000000
      -1.000000        7.000000       -4.000000        2.000000
       3.000000       -1.000000        4.000000        2.000000

vector b =
      12.000000
     -16.000000
      14.000000
      10.000000
      
================= Output =================
matrix U =
       3.000000       -1.000000        4.000000        2.000000
       0.000000        6.666667       -2.666667        2.666667
       0.000000        0.000000       -2.000000        2.000000
       0.000000        0.000000        0.000000       -2.000000

vector d =
      12.000000
     -20.000000
      28.000000
      11.540000

vector x =
      31.340000
      -8.600000
     -19.770000
      -5.770000

P =
       0.000000        0.000000        0.000000        1.000000
       0.000000        0.000000        1.000000        0.000000
       1.000000        0.000000        0.000000        0.000000
       0.000000        1.000000        0.000000        0.000000
```



# backSub()

Solve Upper Matrix 

gaussElim() 함수에서 반환 받은 U행렬과 d행렬을 각각 U와 y에 넣어 Ux=y 방정식에서 x 행렬을 반환한다.

또는 solveLU() 함수 안에서 U 행렬

U --> U

d --> y



```c++
Matrix	backSub(Matrix U, Matrix y, Matrix x);
```



**Parameters**

+ U : Upper Matrix form. (nxn)
+ y : Vector (nx1)
+ x : Answers (nx1)



사용 예제는 다음과 같다.

```c++
Matrix matA = txt2Mat(path, "matA");
Matrix vecb = txt2Mat(path, "vecb");
Matrix matU = zeros(matA.rows, matA.cols);
Matrix vecd = zeros(vecb.rows, vecb.cols);
Matrix P = eye(matA.rows, matA.cols);
Matrix x = createMat(vecd.rows, 1);       	// backSub 함수에서 x행렬을 받아오기 위해 행렬 생성

gaussElim(matA, vecb, matU, vecd, P);
	
Matrix vecx = backSub(matU, vecd, x);		// 새로운 행렬로 받아야 함
```

```
================= Input =================
matrix U =
       3.000000       -1.000000        4.000000        2.000000
       0.000000        6.666667       -2.666667        2.666667
       0.000000        0.000000       -2.000000        2.000000
       0.000000        0.000000        0.000000       -2.000000

vector d =
      12.000000
     -20.000000
      28.000000
      11.540000
      
================= Output =================
vector x =
      31.340000
      -8.600000
     -19.770000
      -5.770000
```



# fwdSub()

Solving Lower Matrix

LUdecomp() 함수에서 반환 받은 L, U행렬 solveLU() 함수 안에서 각각 L과 y에 넣어 Ly=b 방정식에서 y 행렬을 반환한다.

```c++
Matrix fwdSub(Matrix L, Matrix b, Matrix y);
```

**Parameters**

+ U : Upper Matrix form. (nxn)
+ b : Vector (nx1)
+ y : Answers  (nx1)



# LUdecomp()

Decompose Matrix A into L and U

PA = P(LU)

A 행렬을 L 행렬과 U 행렬로 나누는 과정에서 0으로 나뉘는 것을 방지하기 위해 scaled pivoting 기법이 사용되었다.

이때 변화한 주소값 P를 반환하여 solveLU() 함수에서 사용한다.

```c++
void LUdecomp(Matrix A, Matrix L, Matrix U, Matrix P)
```

**Parameters**

+ A : Input Matrix. (nxn)
+ U : Output Matrix (upper triangular matrix) (nxn)
+ L : Output Matrix (Lower triangular matrix) (nxn)
+ P : Pivoting Matrix (Permutation Matrix) (nxn)



사용 예제는 다음과 같다.

```c++
Matrix matA = txt2Mat(path, "matA");
Matrix matU = zeros(matA.rows, matA.cols);
Matrix matL = zeros(matA.rows, matA.cols);
Matrix matP = eye(matA.rows, matA.cols);

LUdecomp(matA, matL, matU, matP);
```

```
================= Input =================
matrix A =
       0.000000        4.000000        2.000000       -2.000000        8.000000
       1.000000        1.000000        2.000000        1.000000        3.000000
       1.000000        2.000000        1.000000        2.000000        2.000000
       2.000000        2.000000        1.000000       -1.000000        4.000000
       1.000000        2.000000        5.000000       -1.000000        4.000000
       
================= Output =================
matrix U =
       1.000000        2.000000        1.000000        2.000000        2.000000
       0.000000       -1.000000        1.000000       -1.000000        1.000000
       0.000000        0.000000       -3.000000       -3.000000       -2.000000
       0.000000        0.000000        0.000000      -12.000000        8.000000
       0.000000        0.000000        0.000000        0.000000       -5.333333

matrix L =
       1.000000        0.000000        0.000000        0.000000        0.000000
       1.000000        1.000000        0.000000        0.000000        0.000000
       2.000000        2.000000        1.000000        0.000000        0.000000
       0.000000       -4.000000       -2.000000        1.000000        0.000000
       1.000000       -0.000000       -1.333333        0.583333        1.000000

matrix P =
       0.000000        0.000000        1.000000        0.000000        0.000000
       0.000000        1.000000        0.000000        0.000000        0.000000
       0.000000        0.000000        0.000000        1.000000        0.000000
       1.000000        0.000000        0.000000        0.000000        0.000000
       0.000000        0.000000        0.000000        0.000000        1.000000
```

# solveLU()

Solve Matrix which is decomposed LU

원래 Ax = b 였던 방정식을 pivoting을 포함한 LU로 분해하였다.

LU에는 pivoting이 적용되었으므로 b에도 P를 곱해줘야 한다.

```c++
void solveLU(Matrix L, Matrix U, Matrix P, Matrix b, Matrix x)
```

**LUx = Pb**

**L(Ux) = Pb --> Ux = y**

**Ly = Pb 방정식에서 fwdSub() 함수를 이용하여 y를 구할 수 있다.**

**그 다음 Ux = y 방정식에서 backSub() 함수를 이용하여 최종적으로 x를 구할 수 있다.**

**Parameters**

+ L : Input Matrix (Lower triangular matrix)
+ U : Input Matrix (upper triangular matrix)
+ P : Pivoting Matrix (Permutation Matrix)
+ b : Vector for fwdSub()
+ x : Vector for backSub() , **FINAL OUTPUT**



사용 예제는 다음과 같다.

```c++
Matrix matA = txt2Mat(path, "matA");
Matrix matU = zeros(matA.rows, matA.cols);
Matrix matL = zeros(matA.rows, matA.cols);
Matrix matP = eye(matA.rows, matA.cols);
Matrix vecb = txt2Mat(path, "vecb");
Matrix vecx = createMat(vecb.rows, vecb.cols);

LUdecomp(matA, matL, matU, matP);
solveLU(matL, matU, matP, vecb, vecx);
```

```
================= Input =================
matrix A =
      75.000000      -20.000000        0.000000
     -20.000000       35.000000      -15.000000
       0.000000      -15.000000       15.000000

vector b =
      19.620001
      29.430000
      14.715000

matrix U =
      75.000000      -20.000000        0.000000
       0.000000       29.666667      -15.000000
       0.000000        0.000000        7.415730

matrix L =
       1.000000        0.000000        0.000000
      -0.266667        1.000000        0.000000
       0.000000       -0.505618        1.000000

================= Output =================
vector x =
       1.159364
       3.366614
       4.347614
```



# eig()

Get EigenValue from Matrix A with **QR decomposition**

QR decomposition를 통해 고유값(Eigen Value)를 보존하면서 A를 Q와 R로 나눈다.

이때 충분히 반복한 R 행렬의 대각 행렬에서 고유값을 반환한다.

```c++
Matrix eig(Matrix A)
```

**Parameters**

+ A : Input Matrix (nxn)



사용 예제는 다음과 같다.

```c++
Matrix matA = txt2Mat(path, "prob_matA");
Matrix eigen = createMat(matA.rows, matA.cols);

eigen = eig(matA);
```

```
------------------------------------------------------------------------------------------
                                        Input Matrix
------------------------------------------------------------------------------------------
A =
      10.000000       32.000000        4.000000
      15.000000        5.000000       42.000000
       1.000000       51.000000       43.000000

------------------------------------------------------------------------------------------
                                        Eigenvalue Results
------------------------------------------------------------------------------------------
EigenValue =
      13.952875      -33.019777       77.066902
```



# eigvec()

Get EigenVector from eig() 

eig() 함수에서 구한 고유값으로부터 고유 벡터를 구하는 함수이다.

```c++
Matrix eigvec(Matrix A)
```

**Parameters**

+ A : Input Matrix (nxn)



사용 예제는 다음과 같다.

```c++
Matrix matA = txt2Mat(path, "prob_matA");
Matrix eigenvec = createMat(matA.rows, matA.cols);

eigenvec = eigvec(matA);
```

```
------------------------------------------------------------------------------------
                                  Input Matrix
------------------------------------------------------------------------------------
A =
      10.000000       32.000000        4.000000
      15.000000        5.000000       42.000000
       1.000000       51.000000       43.000000
       
-------------------------------------------------------------------------------------
                               Eigenvector Results
-------------------------------------------------------------------------------------
Eigen Vector =
       0.967396       -0.494518        0.298550
      -0.242026        0.724781        0.526183
      -0.074624       -0.479734        0.796241
```



# ode()

주어진 미분방정식의 그래프를 그린다.

```c++
void ode(double myfunc(const double t, const double y), double y[], double t0, double tf, double h, double y0, int method);
```

**Parameters**

+ myfunc() : 주어진 미분 방정식
+ y : 받아올 y 좌표 배열
+ t0 : 시작 지점
+ tf : 종료 지점
+ h : t 간격
+ y0 : y의 초기값
+ method : 사용할 ode solver 



사용 예제는 다음과 같다.

![image](https://user-images.githubusercontent.com/118132313/206867372-d26eafd6-db7c-4084-a463-86fdaa755982.png)

```c++
double myfunc(const double t, const double y);
int main(int argc, char* argv[]) {

	double t0 = 0;
	double tf = 0.1;
	double h = 0.001;
	double y[101] = { 0 };
	double y0 = 0;

	ode(myfunc, y, t0, tf, h, y0, EU);
	ode(myfunc, y, t0, tf, h, y0, EM);
	ode(myfunc, y, t0, tf, h, y0, RK2);
	ode(myfunc, y, t0, tf, h, y0, RK3);

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
} 												// 예) 미분 방정식 
```



## odeEU()

![image](https://user-images.githubusercontent.com/118132313/206868304-1aa1f42f-ed0e-4dd3-b504-f716690605bf.png)

odeEU() is Euler method . One of the method that solves ODE problem

```c++
void odeEU(double myfunc(const double t, const double y), double y[], double t0, double tf, double h, double y0);
```

**Parameters**

+ myfunc() : 주어진 미분 방정식
+ y : 받아올 y 좌표 배열
+ t0 : 시작 지점
+ tf : 종료 지점
+ h : t 간격
+ y0 : y의 초기값



## odeEM()

![image](https://user-images.githubusercontent.com/118132313/206868688-d5fa8ec6-b837-4b20-a1a9-d4051dc89b1b.png)

odeEM() is Euler Modified method . One of the method that solves ODE problem

```c++
void odeEM(double myfunc(const double t, const double y), double y[], double t0, double tf, double h, double y0);
```

**Parameters**

+ myfunc() : 주어진 미분 방정식
+ y : 받아올 y 좌표 배열
+ t0 : 시작 지점
+ tf : 종료 지점
+ h : t 간격
+ y0 : y의 초기값



## odeRK2

![image](https://user-images.githubusercontent.com/118132313/206868856-23a8d794-028d-477a-807e-bc3417cfb2e7.png)

odeRK2() is Runge-Kutta 2nd order method . One of the method that solves ODE problem

```c++
void odeRK2(double myfunc(const double t, const double y), double y[], double t0, double tf, double h, double y0);
```

**Parameters**

+ myfunc() : 주어진 미분 방정식
+ y : 받아올 y 좌표 배열
+ t0 : 시작 지점
+ tf : 종료 지점
+ h : t 간격
+ y0 : y의 초기값



## odeRK3()

![image-20221211030026509](C:\Users\82107\AppData\Roaming\Typora\typora-user-images\image-20221211030026509.png)

odeRK3() is Runge-Kutta 3rd order method . One of the method that solves ODE problem

```c++
void odeRK3(double myfunc(const double t, const double y), double y[], double t0, double tf, double h, double y0);
```

**Parameters**

+ myfunc() : 주어진 미분 방정식
+ y : 받아올 y 좌표 배열
+ t0 : 시작 지점
+ tf : 종료 지점
+ h : t 간격
+ y0 : y의 초기값



# sys2RK2

![image](https://user-images.githubusercontent.com/118132313/206869202-5637f4e5-3ca3-4030-a473-f98328623499.png)

sys2RK2() is method for solving 2nd order ODE.

This function uses Runge-Kutta 2nd order method.

![image](https://user-images.githubusercontent.com/118132313/206869710-801a1c0b-aba6-4ab1-8fef-381ec6adce5e.png)

```c++
void sys2RK2(void mckfunc(const double t, const double Y[], double dYdt[]), double y1[], double y2[], double t0, double tf, double h, double y1_init, double y2_init);
```

**Parameters**

+ myfunc() : 주어진 미분 방정식
  + Y[0] = y(t)
  + Y[1] = y'(t)
+ y1 : 2계 미분방정식에서 y의 해 배열
+ y2 : 2계 미분방정식에서 y'의 해 배열
+ t0 : 시작 지점
+ tf : 종료 지점
+ h : t 간격
+ y1_init : y의 초기값 = y(0)
+ y2_init : y'의 초기값 = y'(0)



사용 예제는 다음과 같다.

```c++
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
		printf("t = %f \n y(t) = %f \n z(t) = %f \n\n", t + k * h, y1[k], y2[k]);
		//printf("%f\n", y1[k]);
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

	double Fin = A * cos(2 * PI * f * t);								// 예제
	/*double Fin = A;*/
	/*double Fin = 0;		*/					
																		// Y[0]  = y(t) , Y[1] = z(t)
	dYdt[0] = Y[1];														// dydt == z(t)
	dYdt[1] = (-k * Y[0] - c * Y[1] + Fin) / m;							// dzdt
}
```

```
=============OUTPUT=============
t = 1.000000
y(t) = 0.010953
z(t) = -0.000956
```



# linearRegression()

linearRegression()은 주어진 좌표를 이용하여 추세선을 그려 다음 값을 예측하는 함수이다.

이는 total residual error를 최소화 하는 기울기 a0와 y절편 a1을 구한다.

![image](https://user-images.githubusercontent.com/118132313/206870250-0f7ce064-3fbf-412b-8829-2b25cc0acb72.png)



```c++
double linearRegression(double x[], double y[], int m, double z[]);
```

**Parameters**

+ x : 주어진 x 좌표들의 배열 (x와 y의 배열 요소 수는 같다)
+ y : 주어진 y 좌표들의 배열 (x와 y의 배열 요소 수는 같다)
+ m : x와 y의 배열 수
+ z : 받아올 a1 과 a0
  + z[0] = a1
  + z[1] = a0

$$
y = a_0x + a_1
$$

사용 예제는 다음과 같다.

```c++
int main(int argc, char* argv[]) {

	printf("	    < Linear_Regression	> \n\n");

	double xi[] = { 30, 40, 50, 60, 70, 80 };
	double yi[] = { 1.05, 1.07, 1.09, 1.14, 1.17, 1.21 };
	double z[2] = { 0 };
	int m = sizeof(xi) / sizeof(xi[0]);
	
	linearRegression(xi, yi, m, z);

	printf("a1 = z[0] = %f \na0 = z[1] = %f \n\nP = a0 * T + a1\n\n", z[0], z[1]);		// a1 = z[0], a0 = z[1]
	double T = 100;
	double Pr = 0;
	Pr = z[1] * T + z[0];																// a0*T + a1
	printf("Pressure  = %f\n\n", Pr);
    
    system("pause");
	return 0;
}
```

```
< Linear_Regression >

a1 = z[0] = 0.940952
a0 = z[1] = 0.003286

P = a0 * T + a1

Pressure  = 1.269524
```



# polyfit()

polyfit() 또한 주어진 좌표를 이용하여 n차 방정식의 계수들을 구하는 함수이다.

![image](https://user-images.githubusercontent.com/118132313/206871251-7f203b06-b7d1-4cdc-add2-807662c8d2e0.png)

![image](https://user-images.githubusercontent.com/118132313/206871296-c61b75b2-06c8-44af-aa73-76f99ad5ea23.png)

```c++
void polyfit(Matrix x, Matrix y, Matrix z, int n);
```

**Parameters**

+ x : 주어진 x 좌표들
+ y : 주어진 y 좌표들
+ z : 받아올 n차 방정식의 계수들
+ n : 'n'차 방정식



사용 예제는 다음과 같다.

```c++
int main(int argc, char* argv[]) {
    
	printf("		 < Polyfit >				   \n\n");

	double yt[] = { 0, 3, 4.5, 5.8, 5.9, 5.8, 6.2, 7.4, 9.6, 15.6, 20.7, 26.7, 31.1, 35.6, 39.3, 41.5 };
	int N = sizeof(yt) / sizeof(double);
	int n = 4;

	Matrix Y = createMat(N, 1);
	Matrix X = zeros(N, 1);
	Matrix Z = zeros(n + 1, 1);
	for (int i = 0; i < N; i++) {
		X.at[i][0] = 0.4 * i;
		Y.at[i][0] = yt[i];
	}

	polyfit(X, Y, Z, n);

	printf("polynomial : [ %d ]\n\n", n);
	printMat(Z, "z");

	system("pause");
	return 0;
}
```

```
  < Polyfit >

iteration = 67

polynomial : [ 4 ]

z =
      -0.264389
       3.118549
     -10.192668
      12.877980
      -0.274607
```



# newtonRoot()

newtonRaphson 방법으로 행렬의 성질을 이용하여 해를 구하는 함수이다.

![image](https://user-images.githubusercontent.com/118132313/206871691-1bb5c117-30b5-427e-91e1-ea4c7b5e1ef6.png)

```c++
void newtonRoot(Matrix X, void myfunc(Matrix X, Matrix F), void myjacob(Matrix X, Matrix J), double tol);
```

**Parameters**

+ x : 주어진 초기값들
+ myfunc : 주어진 function
+ myjacob : myfunc의 jacobian function
+ tol : 오차 허용치



사용 예제는 다음과 같다.

```c++
void myjacob(Matrix X, Matrix J);
void myfunc(Matrix X, Matrix F);

int main(int argc, char* argv[]) {

	double xval[] = { 10.0 / 180 * PI, 10, 10 };
	Matrix X = createMat(3, 1);
	X = arr2Mat(xval, 3, 1);

	newtonRoot(X, myfunc, myjacob, tol);
	
	system("pause");
	return 0;
}

void myjacob(Matrix X, Matrix J) {
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

void myfunc(Matrix X, Matrix F) {
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
```

```
=========Input===========
X0 =
       0.174533
      10.000000
      10.000000
=========Output===========
X =
       0.523599
     100.000000
      99.999960
```

