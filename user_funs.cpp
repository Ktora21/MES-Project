#include <math.h>
#include"user_funs.h"
#include <iostream>
#include "user_structs.cpp"
using namespace std;

double valueInOnePoint(double x)
{
	return 2 * x * x + 3 * x - 8;
}
double valueInOnePoint1(double x)
{
	return x * x - 3 * x + 6;
}

double valueInTwoPoint(double x, double y)
{
	return -5 * x * x * y + 2 * x * y * y + 10;
};

double value(double pc[], double w[], int variable, int point)
{
	double value = 0;
	if (variable == 1) {
		for (int i = 0; i < point; i++)
			value += valueInOnePoint(pc[i]) * w[i];
		return value;
	}
	if (variable == 2) {
		for (int i = 0; i < point; i++) {
			for (int j = 0; j < point; j++)
				value += valueInTwoPoint(pc[i], pc[j]) * w[i] * w[j];
		}
		return value;
	}
}

double value1(double pc[], double w[], int point, double a, double b)  // obliczanie wartoœci dla lokalnego ukladu
{
	double value = 0;
	double detJ = (b - a) / 2;
	for (int i = 0; i < point; i++) {
		pc[i] = ((1 - pc[i]) / 2) * a + ((pc[i] + 1) / 2) * b;
		value += (valueInOnePoint1(pc[i]) * w[i]);
	}
	value = value * detJ;
	return value;
}

//obliczanie pochodnych po ksi lub eta (eta ==1 to matrixksi , N okresla ktora funkcja kszta³tu
double derivate(int eta, int N, double x) {
	if (eta == 1 && N == 1)
		return -1.0 / 4.0 * (1 - x);
	else if (eta == 1 && N == 2)
		return 1.0 / 4.0 * (1 - x);
	else if (eta == 1 && N == 3)
		return 1.0 / 4.0 * (1 + x);
	else if (eta == 1 && N == 4)
		return -1.0 / 4.0 * (1 + x);

	else if (eta == 0 && N == 1)
		return -1.0 / 4.0 * (1 - x);
	else if (eta == 0 && N == 2)
		return -1.0 / 4.0 * (1 + x);
	else if (eta == 0 && N == 3)
		return 1.0 / 4.0 * (1 + x);
	else if (eta == 0 && N == 4)
		return 1.0 / 4.0 * (1 - x);
	else
		return 0.0;
}
// obliczanie wartosci funkcji kszta³tu
double calculatingNFunction(int N, double eta, double ksi) {
	if (N == 1)
		return 1.0 / 4.0 * (1 - ksi) * (1 - eta);
	else if (N == 2)
		return 1.0 / 4.0 * (1 + ksi) * (1 - eta);
	else if (N == 3)
		return 1.0 / 4.0 * (1 + ksi) * (1 + eta);
	else if (N == 4)
		return 1.0 / 4.0 * (1 - ksi) * (1 + eta);
	else
		return 0.0;
}
// punkty ca³kowania dla 1D
double* pc2() {
	double pc2[2] = { (-1.0) / sqrt(3) , 1.0 / sqrt(3) };
	return pc2;
}
double* pc3() {
	double pc3[3] = { -sqrt(3.0 / 5.0), 0, sqrt(3.0 / 5.0) };
	return pc3;
}
double* pc4() {
	double pc4[4] = { -0.861136, -0.339981, 0.339981, 0.861136 };
	return pc4;
}
double* pc5() {
	double pc5[5] = { -0.906180, -0.538469, 0, 0.538469, 0.906180 };
	return pc5;
}

// punkty ca³kowania dla 2D - wartosci Eta i ksi
double* eta4() {
	double eta4[4] = { (-1.0) / sqrt(3), (-1.0) / sqrt(3), 1.0 / sqrt(3), 1.0 / sqrt(3) };
	return eta4;
}
double* ksi4() {
	double ksi4[4] = { (-1.0) / sqrt(3),  1.0 / sqrt(3), (-1.0) / sqrt(3), 1.0 / sqrt(3) };
	return ksi4;
}

double* eta9() {
	double eta9[9] = { -sqrt(3.0 / 5.0), -sqrt(3.0 / 5.0), -sqrt(3.0 / 5.0), 0, 0, 0,sqrt(3.0 / 5.0), sqrt(3.0 / 5.0), sqrt(3.0 / 5.0) };
	return eta9;
}
double* ksi9() {
	double ksi9[9] = { -sqrt(3.0 / 5.0), 0 , sqrt(3.0 / 5.0),  -sqrt(3.0 / 5.0), 0 , sqrt(3.0 / 5.0),  -sqrt(3.0 / 5.0), 0 , sqrt(3.0 / 5.0) };
	return ksi9;
}

double* eta16() {
	double eta16[16] = { -0.861136 ,-0.861136 ,-0.861136 ,-0.861136 , -0.339981 ,-0.339981 ,-0.339981 ,-0.339981 ,0.339981 ,0.339981 ,0.339981 ,0.339981, 0.861136,0.861136,0.861136,0.861136 };
	return eta16;
}
double* ksi16() {
	double ksi16[16] = { -0.861136, -0.339981, 0.339981, 0.861136,  -0.861136, -0.339981, 0.339981, 0.861136,  -0.861136, -0.339981, 0.339981, 0.861136,  -0.861136, -0.339981, 0.339981, 0.861136 };

	return ksi16;
}
double* eta25() {
	double eta25[25] = { -0.906180, -0.906180, -0.906180, -0.906180, -0.906180, -0.538469, -0.538469, -0.538469, -0.538469, -0.538469, 0, 0, 0, 0, 0, 0.538469, 0.538469, 0.538469, 0.538469, 0.538469, 0.906180, 0.906180, 0.906180, 0.906180, 0.906180 };
	return eta25;
}
double* ksi25() {
	double ksi25[25] = { -0.906180, -0.538469, 0, 0.538469, 0.906180, -0.906180, -0.538469, 0, 0.538469, 0.906180, -0.906180, -0.538469, 0, 0.538469, 0.906180, -0.906180, -0.538469, 0, 0.538469, 0.906180, -0.906180, -0.538469, 0, 0.538469, 0.906180 };
	return ksi25;
}

//wartoœci przemno¿onych przez siebie wag
double* multiplyWeights4()
{
	double weights[4] = { 1,1,1,1 };
	return weights;
}
double* multiplyWeights9()
{
	double weights[9] = { (5.0 / 9.0) * (5.0 / 9.0),(5.0 / 9.0) * (8.0 / 9.0),(5.0 / 9.0) * (5.0 / 9.0),(5.0 / 9.0) * (8.0 / 9.0),(8.0 / 9.0) * (8.0 / 9.0),(5.0 / 9.0) * (8.0 / 9.0),(5.0 / 9.0) * (5.0 / 9.0), (5.0 / 9.0) * (8.0 / 9.0), (5.0 / 9.0) * (5.0 / 9.0) };
	return weights;
}
double* multiplyWeights16()
{
	double weights[16] = { 0.347855 * 0.347855, 0.652145 * 0.347855,0.652145 * 0.347855, 0.347855 * 0.347855,
							0.347855 * 0.652145, 0.652145 * 0.652145 ,0.652145 * 0.652145 ,0.347855 * 0.652145,
							0.347855 * 0.652145,0.652145 * 0.652145 ,0.652145 * 0.652145 , 0.347855 * 0.652145,
							0.347855 * 0.347855, 0.652145 * 0.347855, 0.652145 * 0.347855, 0.347855 * 0.347855 };
	return weights;
}
double* multiplyWeights25() {
	double weights[25] = { 0.236927 * 0.236927, 0.478629 * 0.236927, 0.568889 * 0.236927, 0.478629 * 0.236927, 0.236927 * 0.236927,
						   0.236927 * 0.478629, 0.478629 * 0.478629, 0.568889 * 0.478629, 0.478629 * 0.478629, 0.236927 * 0.478629,
						   0.236927 * 0.568889, 0.478629 * 0.568889, 0.568889 * 0.568889, 0.478629 * 0.568889, 0.236927 * 0.568889,
						   0.236927 * 0.478629, 0.478629 * 0.478629, 0.568889 * 0.478629, 0.478629 * 0.478629, 0.236927 * 0.478629,
						   0.236927 * 0.236927, 0.478629 * 0.236927, 0.568889 * 0.236927, 0.478629 * 0.236927, 0.236927 * 0.236927 };
	return weights;
}

//Obliczanie pochodnych X po ksi, gdzie jako parametry jest tablica z dN/dKsi, tablica wartosci x i y dla jednego elementu
double dxdKsi(double** oneGrid, double** matrixKsi, int pc) {
	double result = 0;
	for (int i = 0; i < 4; i++) {
		result += matrixKsi[pc][i] * oneGrid[0][i];
	}
	return result;
}
//Obliczanie pochodnych X po eta, gdzie jako parametry jest tablica z dN/dKsi, tablica wartosci x i y dla jednego elementu
double dxdEta(double** oneGrid, double** matrixEta, int pc) {
	double result = 0;
	for (int i = 0; i < 4; i++) {
		result += matrixEta[pc][i] * oneGrid[0][i];
	}
	return result;
}

//Obliczanie pochodnych Y po ksi, gdzie jako parametry jest tablica z dN/dKsi, tablica wartosci x i y dla jednego elementu
double dydKsi(double** oneGrid, double** matrixKsi, int pc) {
	double result = 0;
	for (int i = 0; i < 4; i++) {
		result += matrixKsi[pc][i] * oneGrid[1][i];
	}
	return result;
}

//Obliczanie pochodnych Y po ksi, gdzie jako parametry jest tablica z dN/dKsi, tablica wartosci x i y dla jednego elementu
double dydEta(double** oneGrid, double** matrixEta, int pc) {
	double result = 0;
	for (int i = 0; i < 4; i++) {
		result += matrixEta[pc][i] * oneGrid[1][i];
	}
	return result;
};

//czêœciowe obliczenia dla Jakobianu przekszta³ceñ - obliczane s¹ ju¿ wartoœci po odwroceniu lecz nie pomno¿one przez 1/det
double*** matrixJacobi(double** oneGrid, double** matrixEta, double** matrixKsi, int pc) {
	double*** matrixs = new double** [pc];
	for (int i = 0; i < pc; i++) {
		matrixs[i] = new double* [2];
		for (int j = 0; j < 2; j++) {
			matrixs[i][j] = new double[2];
		}
		matrixs[i][0][0] = dydEta(oneGrid, matrixEta, i);
		matrixs[i][0][1] = (-1.0) * dydKsi(oneGrid, matrixKsi, i);
		matrixs[i][1][0] = (-1.0) * dxdEta(oneGrid, matrixEta, i);
		matrixs[i][1][1] = dxdKsi(oneGrid, matrixKsi, i);
	}
	return matrixs;
}
// przemno¿enie wartoœci przez odwrotnosc wyznacznika
double*** jacobian(double*** matrixJacobian, int pc) {
	double det = 0;
	double*** jacobian = new double** [pc];
	for (int i = 0; i < pc; i++) {
		det = matrixJacobian[i][0][0] * matrixJacobian[i][1][1] - (-1) * matrixJacobian[i][1][0] * (-1) * matrixJacobian[i][0][1];
		jacobian[i] = new double* [2];
		for (int j = 0; j < 2; j++) {
			jacobian[i][j] = new double[2];
		}
		jacobian[i][0][0] = 1.0 / det * matrixJacobian[i][0][0];
		jacobian[i][0][1] = 1.0 / det * matrixJacobian[i][0][1];
		jacobian[i][1][0] = 1.0 / det * matrixJacobian[i][1][0];
		jacobian[i][1][1] = 1.0 / det * matrixJacobian[i][1][1];
	}
	return jacobian;
}

// Obliczanie pochodniej dN po X gdzie podawany jest jakobian przekszta³cenia i macierze pochodnych kszta³tu po  ksi i eta
double** dNdX(double*** jacobian, double** matrixEta, double** matrixKsi, int pc) {
	double** X = new double* [pc];
	for (int i = 0; i < pc; i++) {
		X[i] = new double[4];
		for (int j = 0; j < 4; j++) {
			X[i][j] = jacobian[i][0][0] * matrixKsi[i][j] + jacobian[i][0][1] * matrixEta[i][j];
		}
	}

	return X;
}
// Obliczanie pochodniej dN po Y gdzie podawany jest jakobian przekszta³cenia i macierze pochodnych kszta³tu po  ksi i eta
double** dNdY(double*** jacobian, double** matrixEta, double** matrixKsi, int pc) {
	double** X = new double* [pc];
	for (int i = 0; i < pc; i++) {
		X[i] = new double[4];
		for (int j = 0; j < 4; j++) {
			X[i][j] = jacobian[i][1][0] * matrixKsi[i][j] + jacobian[i][1][1] * matrixEta[i][j];
		}
	}
	return X;
}

// Obliczenie macierzy H dla ka¿dego boku
double*** hPc(double*** jacobian, double** dNdX, double** dNdY, int pc, double kt) {
	double*** hPc = new double** [pc];
	double det = 0;
	for (int i = 0; i < pc; i++) {
		hPc[i] = new double* [4];
		for (int j = 0; j < 4; j++) {
			hPc[i][j] = new double[4];
		}
	}
	for (int i = 0; i < pc; i++) {
		for (int j = 0; j < 4; j++) {
			for (int k = 0; k < 4; k++) {
				det = jacobian[i][0][0] * jacobian[i][1][1] - jacobian[i][1][0] * jacobian[i][0][1];
				hPc[i][j][k] = kt * (dNdX[i][j] * dNdX[i][k] + dNdY[i][j] * dNdY[i][k]) * det;
			}
		}
	}
	return hPc;
}

//Obliczenie macierzy C dla ka¿dego z boków elementów
double*** MatrixCInPc(double*** jacobian, double** UniversalInOneNode, int pc, double specificHeat, double density) {
	double*** CInPc = new double** [pc];
	double det = 0;
	for (int i = 0; i < pc; i++) {
		CInPc[i] = new double* [4];
		for (int j = 0; j < 4; j++) {
			CInPc[i][j] = new double[4];
		}
	}
	for (int i = 0; i < pc; i++) {
		for (int j = 0; j < 4; j++) {
			for (int k = 0; k < 4; k++) {
				det = jacobian[i][0][0] * jacobian[i][1][1] - jacobian[i][1][0] * jacobian[i][0][1];
				CInPc[i][j][k] = specificHeat * density * (UniversalInOneNode[i][j] * UniversalInOneNode[i][k]) * det;
				//cout << CInPc[i][j][k] << " ";
			}
			//cout << endl;
		}
		//cout << endl;
	}
	return CInPc;
}

//Obliczenia dla fina³owych macierzy zarówno H jak i C dla ka¿dego z elementów
//tutaj nastpeuje dodanie do siebie macierzy dla boków do macierzy dla elementów i wymno¿enie przez wagi
double** finalMatrix(double*** hPc, double w[], int pc)
{
	double** H = new double* [4];
	for (int i = 0; i < 4; i++) H[i] = new double[4];
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++)
			H[i][j] = 0;
	}
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			for (int k = 0; k < pc; k++) {
				H[i][j] += hPc[k][i][j] * w[k];
			}
		}
	}
	return H;
}

//Metoda elmiminacji Gaussa
double* GaussElimination(double** HG, double* P, int nN)
{
	double** result = new double* [nN];
	for (int i = 0; i < nN; i++) result[i] = new double[nN + 1];
	for (int i = 0; i < nN; i++) {
		for (int j = 0; j < nN; j++)
			result[i][j] = HG[i][j];
	}
	for (int i = 0; i < nN; i++) {
		result[i][nN] = P[i];
	}
	//Wypisywanie HGAndP - polaczenie agregacji HG oraz agregacji wektora P
	//for (int i = 0; i < nN; i++) {
	//	for (int j = 0; j < nN + 1; j++)
	//		cout << HGAndP[i][j] << " ";
	//	cout << endl;
	//}

	double* Temp = new double[nN];
	for (int i = 0; i < nN; i++) Temp[i] = 0.0;
	double m, s;
	int i, j, k;

	for (i = 0; i < nN - 1; i++) {
		for (j = i + 1; j < nN; j++) {
			if (fabs(result[i][i]) < 1e-12) break;
			m = -result[j][i] / result[i][i];
			for (k = i + 1; k <= nN; k++)
				result[j][k] += m * result[i][k];
		}
	}
	for (i = nN - 1; i >= 0; i--) {
		s = result[i][nN];
		for (j = nN - 1; j >= i + 1; j--)
			s -= result[i][j] * Temp[j];
		if (fabs(result[i][i]) < 1e-12) break;
		Temp[i] = s / result[i][i];
	}
	return Temp;
}

//Funkcja szukaj¹ca minimum
double findMin(double* T, int size) {
	double min = T[0];
	for (int i = 0; i < size; i++)
	{
		if (min > T[i]) min = T[i];
	}
	return min;
}

//Funkcja szukaj¹ca minimum
double findMax(double* T, int size) {
	double max = T[0];
	for (int i = 0; i < size; i++)
	{
		if (max < T[i]) max = T[i];
	}
	return max;
}