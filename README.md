# **<span style="color:blue">Numerical Programming_2022</span>**



## <span style="background-color:yellow">myMatrix.cpp</span>



### createMat()

Create Matrix

```
Matrix	createMat(int _rows, int _cols);
```

### freeMat()

Deallocate Memory

```
void	freeMat(Matrix _A);
```

#### Parameters

+ _A : Matrix form. (Should be nxn square.)

### txt2Mat()

Bring .txt File to Matrix

```
Matrix	txt2Mat(std::string _filePath, std::string _fileName);
```

### printMat()

Print the Matrix

```
void	printMat(Matrix _A, const char* _name);
```

#### Parameters

+ _A : Matrix form.
+ char*_name: Explaining what this Matrix is.

### initMat()

Initialize Matrix elements

```
void	initMat(Matrix _A, double _val);
```

#### Parameters

+ _A : Matrix form. (Should be nxn square.)
+ _val : Initializing Matrix element with this value.

### zeros()

Create matrix of all zeros

```
Matrix	zeros(int _rows, int _cols);
```

#### Parameters

+ _rows : Set the number of rows
+ _cols : Set the number of cols

### eye()

Create identity Matrix

```
Matrix	eye(int _rows, int _cols);
```

#### Parameters

+ _rows : Set the number of rows
+ _cols : Set the number of cols

### transpose()

Transpose the Matrix

```
Matrix	transpose(Matrix _A);
```

#### Parameters

+ _A : Matrix form. (Should be nxn square.)

### copyMat()

Copy the Matrix

```
Matrix	copyMat(Matrix _A)
```

#### Parameters

+ _A : Matrix form. (Should be nxn square.)

### copyVal()

Copy the Matrix elements from _A to _B

```
void	copyVal(Matrix _A, Matrix _B);
```

#### Parameters

+ _A : Matrix which want to copy.  

+ _B : Matrix which want to paste.

   ** A  matrix's number of rows and columns should be same with _B*

### multiply()

Multiply Matrix _A and _B  

ex) P = AB 

```
Matrix multiply(Matrix _A, Matrix _B);
```

#### Parameters

+ _A : Matrix form. (Should be nxn Matrix.)

+ _B : Matrix form. (Should be nxn Matrix.)

   ** A  matrix's number of rows and columns should be same with _B*

### subMat()

Subtract Matrix _A - _B

ex) P = A - B

```
Matrix subMat(Matrix _A, Matrix _B);
```

#### Parameters

+ _A : Matrix form. (Should be nxn Matrix.)

+ _B : Matrix form. (Should be nxn Matrix.)

   ** A  matrix's number of rows and columns should be same with _B*

### xMat()

Multiply scalar with Matrix

ex) P = c * A

```
Matrix xMat(Matrix _A, double _c);
```

#### Parameters

+ _A : Matrix form.
+ _c : Scalar which want to multiply

### volMat()

Get volume of Matrix

```
double volMat(Matrix _A);
```

#### Parameters

+ _A : Matrix form.
