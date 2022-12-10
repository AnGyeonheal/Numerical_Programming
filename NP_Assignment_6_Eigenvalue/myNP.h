/*----------------------------------------------------------------\
@ Numerical Programming by Young-Keun Kim - Handong Global University

Author           : [An Gyeonheal]
ID				 : 21900416
Created          : 26-03-2018
Modified         : 14-11-2022
Language/ver     : C++ in MSVS2019

Description      : myNP.h
----------------------------------------------------------------*/

#ifndef		_MY_NP_H		// use either (#pragma once) or  (#ifndef ...#endif)
#define		_MY_NP_H

#include "myMatrix.h"

// Matrix addition
extern	Matrix	addMat(Matrix _A, Matrix _B);

// Apply back-substitution
extern	Matrix	backSub(Matrix U, Matrix y, Matrix x);

// Apply forward-substitution
extern Matrix fwdSub(Matrix L, Matrix b, Matrix y);

//
extern double invMat(Matrix A, Matrix Ainv);

//// Gauss-Elimination with partial pivoting
extern void gaussElim(Matrix A, Matrix b, Matrix U, Matrix d, Matrix P);

// LU decomposition with scaled pivoting
extern void LUdecomp(Matrix A, Matrix L, Matrix U, Matrix P); 

extern void solveLU(Matrix L, Matrix U, Matrix P, Matrix b, Matrix x);

extern Matrix eig(Matrix A);

extern Matrix eigvec(Matrix A);
#endif