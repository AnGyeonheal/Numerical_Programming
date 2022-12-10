/*-------------------------------------------------------------------------------\
@ Numerical Programming by Young-Keun Kim - Handong Global University
Author           : [An Gyeonheal]
ID				 : 21900416
Created          : 19-10-2022
Modified         : 07-11-2022
Language/ver     : C++ in MSVS2022
Description      : Assignment_gaussElim_21900416.cpp
-------------------------------------------------------------------------------*/

#define Assignment	5		// enter your assignment number
#define eval		0		// set 0

#include "../../../include/myNP.h"
#include "../../../include/myMatrix.h"

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
	Matrix prob1_matA = txt2Mat(path, "prob1_matA");
	Matrix prob1_vecb = txt2Mat(path, "prob1_vecb");

	/* =============== Q2.Determine the displacement (u1, u2, u3) =============== */
	Matrix prob2_matA = txt2Mat(path, "prob2_matA");
	Matrix prob2_vecb = txt2Mat(path, "prob2_vecb");

	/* =============== Q3. Create a function that find the inverse of A =============== */
	Matrix prob3_matA = txt2Mat(path, "prob3_matA");

	/*==========================================================================*/
	/*					Apply your numerical method algorithm					*/
	/*==========================================================================*/

	/// ====================================  Q1  =====================================

	Matrix prob1_matU = zeros(prob1_matA.rows, prob1_matA.cols);
	Matrix prob1_matL = zeros(prob1_matA.rows, prob1_matA.cols);
	Matrix prob1_matP = eye(prob1_matA.rows, prob1_matA.cols);

	printf("\n================= Input of Q1 =================\n");

	printMat(prob1_matA, "matrix A");
	LUdecomp(prob1_matA, prob1_matL, prob1_matU, prob1_matP);



	/// ====================================  Q2  =====================================

	Matrix prob2_matU = zeros(prob2_matA.rows, prob2_matA.cols);
	Matrix prob2_matL = zeros(prob2_matA.rows, prob2_matA.cols);
	Matrix prob2_matP = eye(prob2_matA.rows, prob2_matA.cols);
	Matrix prob2_vecx = createMat(prob2_vecb.rows, prob2_vecb.cols);

	printf("\n================= Input of Q2 =================\n");

	printMat(prob2_matA, "matrix A");
	printMat(prob2_vecb, "vector b");

	LUdecomp(prob2_matA, prob2_matL, prob2_matU, prob2_matP);

	printMat(prob2_matU, "matrix U");
	printMat(prob2_matL, "matrix L");

	solveLU(prob2_matL, prob2_matU, prob2_matP, prob2_vecb, prob2_vecx);



	/// ====================================  Q3  =====================================

	Matrix Ainv = zeros(prob3_matA.rows, prob3_matA.cols);

	printf("\n================= Input of Q3 =================\n");

	printMat(prob2_matA, "matrix A");

	invMat(prob3_matA, Ainv);
	

	/*==========================================================================*/
	/*							  Print your results							*/
	/*==========================================================================*/

	/// Q1
	printf("\n================= Output of Q1 =================\n");
	
	printMat(prob1_matU, "matrix U");
	printMat(prob1_matL, "matrix L");
	printMat(prob1_matP, "matrix P");

	/// Q2
	printf("\n================= Output of Q2 =================\n");
	
	printMat(prob2_vecx, "vector x");
	
	/// Q3
	printf("\n================= Output of Q3 =================\n");

	printMat(Ainv, "Inverse A");

	/*==========================================================================*/
	/*							  Deallocate memory 							*/
	/*==========================================================================*/

	freeMat(prob1_matA);		
	freeMat(prob2_matA);

	system("pause");
	return 0;
}