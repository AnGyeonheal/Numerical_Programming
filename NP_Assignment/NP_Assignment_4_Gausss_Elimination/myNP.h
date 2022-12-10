/*----------------------------------------------------------------\
@ Numerical Programming by Young-Keun Kim - Handong Global University

Author           : [An Gyeonheal]
ID				 : 21900416
Created          : 26-03-2018
Modified         : 31-10-2022
Language/ver     : C++ in MSVS2019

Description      : myNP.h
----------------------------------------------------------------*/

#ifndef		_MY_NP_H		// use either (#pragma once) or  (#ifndef ...#endif)
#define		_MY_NP_H

#include "myMatrix.h"

// Matrix addition
extern	Matrix	addMat(Matrix _A, Matrix _B);

// Apply back-substitution
extern	Matrix	backSub(Matrix U, Matrix d, Matrix x);

//// Gauss-Elimination with pivoting
void gaussElim(Matrix A, Matrix b, Matrix U, Matrix d, Matrix P);
#endif