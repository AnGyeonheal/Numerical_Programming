# **<span style="color:blue">Numerical Programming_2022</span>**

## <span style="background-color:orange">.cpp Code of < System of NonLinear Equation ></span>

## <span style="background-color:yellow">myNP.cpp</span>

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

### bisection()

Bisection method to get root number

```c++
double bisection(double func(double x), float _a0, float _b0, float _tol)
```

#### Parameters

+ func(double x) : number parameter that in function
+ _a0 : initial number which is left side of root number
+ _b0 : initial number which is right side of root number
+ _tol : tolerance

### newtonRaphson()

NewtonRaphson method to get root number

```c++
double newtonRaphson(double dfunc(double x), double d2func(double x), double _x0, double _tol)
```

#### Parameters

+ dfunc(double x) : number parameter that in differential function
+ d2func(double x) : number parameter that in twice differential function
+ _x0 : root number
+ _tol : tolerance

### gradient1D()

Getting gradient of function that is series of number using Three-point FWD, mid-point and Three-point BWD.

```c++
void gradient1D(double x[], double y[], double dydx[], int m)
```

#### Parameters

+ x[] : number of X axis
+ y[] : number of Y axis
+ dydx[] : gradient number
+ m : times 

### gradientFunc()

Getting gradient from own function using gradient1D() function

```c++
void gradientFunc(double func(const double x), double x[], double dydx[], int m)
```

#### Parameters

+ func(const double x) : function number
+ double x[] : number of X axis
+ dydx[] : gradient number
+ m : times 

### acceleration()

Getting twice differentiated number.

```c++
void acceleration(double x[], double y[], double dy2dx[], int m)
```

#### Parameters

+ x[] : number of X axis
+ y[] : number of Y axis
+ dy2dx[] : twice differentiated number
+ m : times 

### trpaz()

Trapezoidal method for Integration

```c++
double trapz(double x[], double y[], int m)
```

#### Parameters

+ x[] : number of X axis
+ y[] : number of Y axis
+ m : times 

### simpson13()

Simpson method for Integration

```c++
double trapz(double x[], double y[], int m)
```

#### Parameters

+ x[] : number of X axis
+ y[] : number of Y axis
+ m : times 

### Integral()

Integral Using simpson13()

```c++
double integral(double func(const double x), double a, double b, int n)
```

#### Parameters

+ func(const double x) : function number(x)
+ a : initial number (getting start of Integration)
+ b : final number (closing Integration)
+ n : section



## <span style="background-color:orange">.cpp Code of < System of Linear Equation ></span>

## <span style="background-color:yellow">myMatrix.cpp</span>

```c++
#include "myMatrix.h"
```



### createMat()

Create Matrix

```c++
Matrix	createMat(int _rows, int _cols);
```

### freeMat()

Deallocate Memory

```c++
void	freeMat(Matrix _A);
```

#### Parameters

+ _A : Matrix form. (Should be nxn square.)

### txt2Mat()

Bring .txt File to Matrix

```c++
Matrix	txt2Mat(std::string _filePath, std::string _fileName);
```

### printMat()

Print the Matrix

```c++
void	printMat(Matrix _A, const char* _name);
```

#### Parametersc

+ _A : Matrix form.
+ char*_name: Explaining what this Matrix is.

### initMat()

Initialize Matrix elements

```c++
void	initMat(Matrix _A, double _val);
```

#### Parameters

+ _A : Matrix form. (Should be nxn square.)
+ _val : Initializing Matrix element with this value.

### zeros()

Create matrix of all zeros

```c++
Matrix	zeros(int _rows, int _cols);
```

#### Parameters

+ _rows : Set the number of rows
+ _cols : Set the number of cols

### eye()

Create identity Matrix

```c++
Matrix	eye(int _rows, int _cols);
```

#### Parameters

+ _rows : Set the number of rows
+ _cols : Set the number of cols

### transpose()

Transpose the Matrix

```c++
Matrix	transpose(Matrix _A);
```

#### Parameters

+ _A : Matrix form. (Should be nxn square.)

### copyMat()

Copy the Matrix

```c++
Matrix	copyMat(Matrix _A)
```

#### Parameters

+ _A : Matrix form. (Should be nxn square.)

### copyVal()

Copy the Matrix elements from _A to _B

```c++
void	copyVal(Matrix _A, Matrix _B);
```

#### Parameters

+ _A : Matrix which want to copy.  

+ _B : Matrix which want to paste.

  ** A  matrix's number of rows and columns should be same with _B*

### multiply()

Multiply Matrix _A and _B  

ex) P = AB 

```c++
Matrix multiply(Matrix _A, Matrix _B);
```

#### Parameters

+ _A : Matrix form. (Should be nxn Matrix.)

+ _B : Matrix form. (Should be nxn Matrix.)

  ** A  matrix's number of rows and columns should be same with _B*

### subMat()

Subtract Matrix _A - _B

ex) P = A - B

```c++
Matrix subMat(Matrix _A, Matrix _B);
```

#### Parameters

+ _A : Matrix form. (Should be nxn Matrix.)

+ _B : Matrix form. (Should be nxn Matrix.)

  ** A  matrix's number of rows and columns should be same with _B*

### xMat()

Multiply scalar with Matrix

ex) P = c * A

```c++
Matrix xMat(Matrix _A, double _c);
```

#### Parameters

+ _A : Matrix form.
+ _c : Scalar which want to multiply

### volMat()

Get volume of Matrix

```c++
double volMat(Matrix _A);
```

#### Parameters

+ _A : Matrix form.



## <span style="background-color:yellow">myNP.cpp</span>

```c++
#include "myNP.h"
```

### addMat()

Add Matrix _A and _B

ex) P = _A + _B

```c++
Matrix	addMat(Matrix _A, Matrix _B);
```

#### Parameters

+ _A : Matrix form.
+ _B: Matrix form.

### backSub()

Solve Upper Matrix 

```c++
Matrix	backSub(Matrix U, Matrix y, Matrix x);
```

<img src="C:\Users\82107\AppData\Roaming\Typora\typora-user-images\image-20221114143331608.png" style="zoom:80%;" />

<img src="C:\Users\82107\AppData\Roaming\Typora\typora-user-images\image-20221114143409728.png" alt="image-20221114143409728" style="zoom: 80%;" />

#### Parameters

+ U : Upper Matrix form.
+ y : Vector (b1~b4)
+ x : Answers (x1~x4)

### fwdSub()

Solving Lower Matrix

```c++
Matrix fwdSub(Matrix L, Matrix b, Matrix y);
```

#### Parameters

+ U : Lower Matrix form.
+ y : Vector (b1~b4)
+ x : Answers (x1~x4)

### invMat()

Get inverse Matrix

```c++
double invMat(Matrix A, Matrix Ainv);
```

#### Parameters

+ A : Input Matrix which want to get inverse Matrix.
+ Ainv : Output Matrix (inversed Matrix of Matrix A)

### gaussElim()

Transform the linear system into an equivalent form with an upper and lower triangular matrix

```c++
void gaussElim(Matrix A, Matrix b, Matrix U, Matrix d, Matrix P);
```

![image-20221114144748246](C:\Users\82107\AppData\Roaming\Typora\typora-user-images\image-20221114144748246.png)

#### Parameters

+ A : Input Matrix.
+ b : Input Vector.
+ U : Output Matrix (upper triangular matrix)
+ d : Output Vector (d1~d3)
+ P : Pivoting Matrix (Permutation Matrix)

### LUdecomp()

Decompose Matrix A into L and U

ex) A = LU

```c++
void LUdecomp(Matrix A, Matrix L, Matrix U, Matrix P)
```

#### Parameters

+ A : Input Matrix.
+ U : Output Matrix (upper triangular matrix)
+ L : Output Matrix (Lower triangular matrix)
+ P : Pivoting Matrix (Permutation Matrix)

### solveLU()

Solve Matrix which is decomposed LU

```c++
void solveLU(Matrix L, Matrix U, Matrix P, Matrix b, Matrix x)
```

ex ) PA = LU

​		PAx = Pb

​		Ly = pb = d --> fwdSub (Ux = y)

 	   Ux = y --> backSub

<img src="C:\Users\82107\AppData\Roaming\Typora\typora-user-images\image-20221114150532597.png" alt="image-20221114150532597" style="zoom:80%;" />

#### Parameters

+ L : Output Matrix (Lower triangular matrix)
+ U : Output Matrix (upper triangular matrix)
+ P : Pivoting Matrix (Permutation Matrix)
+ b : Vector for fwdSub()
+ x : Vector for backSub()

### eig()

Get EigenValue from Matrix A with QR decomposition

```c++
Matrix eig(Matrix A)
```

#### Parameters

+ A : Input Matrx (Should be nxn Matrix)

### eigvec()

Get EigenVector from Matrix A with QR decomposition

```c++
Matrix eigvec(Matrix A)
```

#### Parameters

+ A : Input Matrx (Should be nxn Matrix)
