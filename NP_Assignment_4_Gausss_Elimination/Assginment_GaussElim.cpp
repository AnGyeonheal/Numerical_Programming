/*-------------------------------------------------------------------------------\
@ Numerical Programming by Young-Keun Kim - Handong Global University
Author           : [An Gyeonheal]
ID				 : 21900416
Created          : 19-10-2022
Modified         : 31-10-2022
Language/ver     : C++ in MSVS2022
Description      : Assignment_gaussElim_21900416.cpp
-------------------------------------------------------------------------------*/

#define Assignment	4		// enter your assignment number
#define eval		0		// set 0

#include "../../../Include/myNP.h"
#include "../../../Include/myMatrix.h"

int main(int argc, char* argv[])
{
	/*	 [¡Ø DO NOT EDIT IT !!!]   Resources file path setting for evaluation	*/
	std::string path = "C:/NP_Matrix_Data/Assignment" + std::to_string(Assignment) + "/";

#if eval
	path += "eval/";
#endif

	/*==========================================================================*/
	/*					Variables declaration & initialization					*/
	/*--------------------------------------------------------------------------*/
	/*   - You can change the variable names									*/
	/*   - However, you must use the specified file name						*/

	/* ======== Q1.Analyze the current in each mesh in the given circuit ======== */
	Matrix matA1 = txt2Mat(path, "prob1_matA");
	Matrix vecb1 = txt2Mat(path, "prob1_vecb");

	/* =============== Q2.Determine the displacement (u1, u2, u3) =============== */
	Matrix matA2 = txt2Mat(path, "prob2_matA");
	Matrix vecb2 = txt2Mat(path, "prob2_vecb");


	/*==========================================================================*/
	/*					Apply your numerical method algorithm					*/
	/*==========================================================================*/

	/// =================================  Q1  ==================================
	Matrix matU1 = zeros(matA1.rows, matA1.cols);
	Matrix vecd1 = zeros(vecb1.rows, vecb1.cols);
	Matrix P1 = eye(matA1.rows, matA1.cols);
	Matrix x1 = createMat(vecd1.rows, 1);
	printf("\n================= Input of Q1 =================\n");

	printMat(matA1, "matrix A");
	printMat(vecb1, "vector b");

	// Gauss-Elimination with pivoting
	gaussElim(matA1, vecb1, matU1, vecd1, P1);

	Matrix vecx1 = backSub(matU1, vecd1, x1);

	/// =================================  Q2  ==================================
	Matrix matU2 = zeros(matA2.rows, matA2.cols);
	Matrix vecd2 = zeros(vecb2.rows, vecb2.cols);
	Matrix P2 = eye(matA2.rows, matA2.cols);
	Matrix x2 = createMat(vecd2.rows, 1);
	printf("\n================= Input of Q2 =================\n");

	printMat(matA2, "matrix A");
	printMat(vecb2, "vector b");

	// Gauss-Elimination with pivoting
	gaussElim(matA2, vecb2, matU2, vecd2, P2);

	Matrix vecx2 = backSub(matU2, vecd2, x2);

	/*==========================================================================*/
	/*							  Print your results							*/
	/*==========================================================================*/

	/// Q1

	printf("\n================= Output of Q1 =================\n");

	printMat(matU1, "matrix U");
	printMat(vecd1, "vector d");
	printMat(vecx1, "vector x");
	printMat(P1, "P");

	/// Q2

	printf("\n================= Output of Q2 =================\n");

	printMat(matU2, "matrix U");
	printMat(vecd2, "vector d");
	printMat(vecx2, "vector x");
	printMat(P2, "P");

	/*==========================================================================*/
	/*							  Deallocate memory 							*/
	/*==========================================================================*/

	freeMat(matA1);		freeMat(vecb1);
	freeMat(matU1);		freeMat(vecd1);		freeMat(vecx1);
	//freeMat(matA2);		freeMat(vecb2);
	//freeMat(matU2);		freeMat(vecd2);		freeMat(vecx2);

	system("pause");
	return 0;
}


