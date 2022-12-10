/*----------------------------------------------------------------\
@ Numerical Programming by Young-Keun Kim - Handong Global University

Author           : [An Gyeonheal]
ID				 : 21900416
Created          : 26-03-2018
Modified         : 14-11-2022
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

	for (k = 0; k < U.rows-1; k++)
	{
		max_a = fabs(U.at[k][k]);
		for (i = k+1; i < U.rows; i++)
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

	for (int k = 0; k < U.rows-1; k++) {
		//for (int i = k; i < U.cols; i++) {
		//	for (int j = 0; j < U.cols; j++) {
		//		if (max_a < fabs(U.at[i][j])) {
		//			max_a = fabs(U.at[i][j]);					// 같은 행에서 max 값 찾기
		//		} 
		//	}
		//	sp.at[0][i] = fabs(U.at[i][k] / max_a);				// 각 행의 sp 값 구하기
		//	max_a = 0;											// max값 초기화
		//}
		//for (int i = 0; i < U.rows; i++){
		//	if (max_sp < sp.at[0][i]) {
		//		max_sp = sp.at[0][i];							// 각 행의 sp 값 중 가장 큰 sp 찾기
		//		t = i;											// 가장 큰 sp를 가진 행의 주소 값 찾기
		//	}
		//}
		//if (k != t) {											// 바꾸는 행과 바꿔지는 행의 sp값이 같을 때
		//	for (int i = 0; i < U.rows; i++) {
		//		temp = U.at[t][i];
		//		U.at[t][i] = U.at[k][i];
		//		U.at[k][i] = temp;								// 기존 행과 sp가 가장 큰 행 바꾸기

		//		temp = P.at[t][i];
		//		P.at[t][i] = P.at[k][i];
		//		P.at[k][i] = temp;								// P 바꾸기

		//		temp = L.at[t][i];
		//		L.at[t][i] = L.at[k][i];
		//		L.at[k][i] = temp;								// L 바꾸기
		//	}
		//}
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
// Solving Eigen
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
		for (l = 0; l < A.rows-1; l++) {
			for (k = 0; k < A.rows-1; k++) {
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
	for (i = 0; i < A.rows; i++) {
		eigenVal.at[0][i] = R.at[A.rows-1-i][A.rows - 1 - i];				// eigenValue 
	}
	return eigenVal;
}
// Solving Eigenvec
Matrix eigvec(Matrix A) {
	Matrix V = zeros(A.rows-1, 1);
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
					vecB.at[i][0] = - B.at[i + 1][k];				// vecB 1
				}
				invMat(subB, invsubB);
				copyVal(multiply(invsubB, vecB),V);					// v1 = inv(subB)*vecB
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
