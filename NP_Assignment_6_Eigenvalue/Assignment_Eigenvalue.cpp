/*------------------------------------------------------------------------------------------\
@ Numerical Methods by Young-Keun Kim - Handong Global University
Author          : SSSLAB
Created         : 01-04-2019
Modified        : 14-11-2022
Language/ver	: C in MSVS2022
Course			: Numerical Programming 2022
Description     : Assignment_eigenvalue_ID.cpp
/------------------------------------------------------------------------------------------*/

#define Assignment	6		// enter your assignment number
#define eval		0		// set 0

#include "../../../Include/myNP.h"
#include "../../../Include/myMatrix.h"

int main(int argc, char* argv[])
{
	/*	 [воик DO NOT EDIT IT !!!]   Resources file path setting for evaluation	*/
	std::string path = "C:/NP_Matrix_Data/Assignment" + std::to_string(Assignment) + "/";

#if eval
	path += "eval/";
#endif

	/*==========================================================================*/
	/*					Variables declaration & initialization					*/
	/*--------------------------------------------------------------------------*/
	/*   - You can change the variable names									*/
	/*   - However, you must use the specified file name						*/
	/*	   : For each assignment, the file name will be notified on HISNET		*/
	/*==========================================================================*/
	Matrix matA = txt2Mat(path, "prob_matA");
	/*==========================================================================*/
	/*					Apply your numerical method algorithm					*/
	/*==========================================================================*/
	Matrix eigen = createMat(matA.rows, matA.cols);
	eigen = eig(matA);
	Matrix eigenvec = createMat(matA.rows, matA.cols);
	eigenvec = eigvec(matA);
	printf("------------------------------------------------------------------------------------------\n");
	printf("					Input Matrix 						\n");
	printf("------------------------------------------------------------------------------------------\n");
	printMat(matA, "A");
	
	/*==========================================================================*/
	/*							  Print your results							*/
	/*==========================================================================*/
	printf("------------------------------------------------------------------------------------------\n");
	printf("					Eigenvalue Results					\n");
	printf("------------------------------------------------------------------------------------------\n");

	printMat(eigen, "EigenValue");

	printf("------------------------------------------------------------------------------------------\n");
	printf("			       Eigenvector Results					\n");
	printf("------------------------------------------------------------------------------------------\n");

	printMat(eigenvec, "Eigen Vector");

	/*==========================================================================*/
	/*							  Deallocate memory 							*/
	/*==========================================================================*/

	freeMat(matA);

	system("pause");
	return 0;
}