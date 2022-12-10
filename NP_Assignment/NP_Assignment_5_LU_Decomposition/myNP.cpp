/*----------------------------------------------------------------\
@ Numerical Programming by Young-Keun Kim - Handong Global University

Author           : [An Gyeonheal]
ID				 : 21900416
Created          : 26-03-2018
Modified         : 07-11-2022
Language/ver     : C++ in MSVS2019

Description      : myNP.cpp
----------------------------------------------------------------*/

#include "myNP.h"


// Matrix addition
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

 // Apply back-substitution
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
		}															// I�� �� ���ͷ� �����
		x = backSub(U, AI, x);										// U * x = AI �� back substitution���� x ���� ���ϱ�
		for (int i = 0; i < A.cols; i++) {
			Uc.at[i][j] = x.at[i][0];
		}initMat(x, 0);												// �� x���͵��� Uc�� �Ҵ�
	}
	for (int j = 0; j < A.cols; j++) {
		for (int i = 0; i < A.cols; i++) {
			AI.at[i][0] = I.at[i][j];
		}															// I�� �� ���ͷ� �����
		x = fwdSub(L, AI, x);										// L * x = AI �� back substitution���� x ���� ���ϱ�
		for (int i = 0; i < A.cols; i++) {
			Lc.at[i][j] = x.at[i][0];
		}initMat(x, 0);												// �� x���͵��� Lc�� �Ҵ�
	}
	copyVal(multiply(Uc, Lc), Ainv);								// copyVal�� ���� Ainv�� ����, Uc * Lc

	return 0;
}
// Gauss-Elimination with pivoting
void gaussElim(Matrix A, Matrix b, Matrix U, Matrix d, Matrix P)
{
	if (A.rows != A.cols || b.rows != A.rows || A.rows != U.rows || A.cols != U.cols || A.rows != P.rows || A.cols != P.cols || d.rows != U.rows || b.rows != d.rows || b.cols != 1 || d.cols != 1) {
		printf("Square Matrix is not Used || Vector's rows are not same with Matrix's rows || Vector's cols are not 1");
		return;
	}															// ���� ó�� (���, ���� ����)

	int i, j, k;
	double max_a;
	double temp;
	double m;

	copyVal(eye(P.rows, P.cols), P);							// Permutation matrix
	copyVal(A, U);
	copyVal(b, d);

	for (k = 0; k < U.rows-1; k++)
	{
		max_a = fabs(U.at[k][k]);
		for (i = k+1; i < U.rows; i++)
		{
			if (max_a < fabs(U.at[i][k]))
			{
				max_a = fabs(U.at[i][k]);						// ���� ū ��� ã��
				for (j = 0; j < U.rows; j++) {
					temp = U.at[i][j];
					U.at[i][j] = U.at[k][j];
					U.at[k][j] = temp;

					temp = P.at[i][j];
					P.at[i][j] = P.at[k][j];
					P.at[k][j] = temp;
				}
			}
		}														// ��� ������
		for (i = k + 1; i < U.rows; i++) {
			if (U.at[k][k] == 0) {
				printf("There's 0 in pivot element");			// ���� ó�� (0���� ���� ��)
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
	copyVal(A, U);												// CopyVal�� ���� ����: ����� ó�� �Է°� ���� ���� �� ��� �󿡼� ���� �Ŀ� �״�� ����ϱ� ���ؼ�

	if (A.rows != A.cols) {
		printf("ERROR: IT IS NOT A N*N MATRIX");
		return;
	}

	for (int k = 0; k < U.rows-1; k++) {
		for (int i = k; i < U.cols; i++) {
			for (int j = 0; j < U.cols; j++) {
				if (max_a < fabs(U.at[i][j])) {
					max_a = fabs(U.at[i][j]);					// ���� �࿡�� max �� ã��
				} 
			}
			sp.at[0][i] = fabs(U.at[i][k] / max_a);				// �� ���� sp �� ���ϱ�
			max_a = 0;											// max�� �ʱ�ȭ
		}
		for (int i = 0; i < U.rows; i++){
			if (max_sp < sp.at[0][i]) {
				max_sp = sp.at[0][i];							// �� ���� sp �� �� ���� ū sp ã��
				t = i;											// ���� ū sp�� ���� ���� �ּ� �� ã��
			}
		}
		if (k != t) {											// �ٲٴ� ��� �ٲ����� ���� sp���� ���� ��
			for (int i = 0; i < U.rows; i++) {
				temp = U.at[t][i];
				U.at[t][i] = U.at[k][i];
				U.at[k][i] = temp;								// ���� ��� sp�� ���� ū �� �ٲٱ�

				temp = P.at[t][i];
				P.at[t][i] = P.at[k][i];
				P.at[k][i] = temp;								// P �ٲٱ�

				temp = L.at[t][i];
				L.at[t][i] = L.at[k][i];
				L.at[k][i] = temp;								// L �ٲٱ�
			}
		}
		for (int i = k + 1; i < U.rows; i++) {
			m = U.at[i][k] / U.at[k][k];
			L.at[i][k] = m;										// L �����
			for (int j = k; j < U.rows + 1; j++) {
				U.at[i][j] = U.at[i][j] - m * U.at[k][j];		// ���콺 �ҰŹ� --> U �����
			}
		}
		initMat(sp, 0);											// sp ���, sp �ִ� �� �ʱ�ȭ
		max_sp = 0;
	}
	for (int i = 0; i < L.rows; i++) {
		L.at[i][i] = 1;
	}
}

void solveLU(Matrix L, Matrix U, Matrix P, Matrix b, Matrix x) {
																// PAx = Pb --> LUx = Pb
	Matrix y = zeros(b.rows, b.cols);							// LUx = b,   Ux = y
	Matrix ny = fwdSub(L, multiply(P, b), y);					// Ly=b (L,U ��� scaled pivoting�� �����Ƿ� vector b ���� pivoting)
	x = backSub(U, ny, x);										// Ux=y

}