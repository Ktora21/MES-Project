#include <iostream>
#include <fstream>
#include <math.h>
#include "user_structs.cpp"
#include "user_funs.h";
using namespace std;
// pobieranie z pliku tekstowego
void LoadData(string filename, Node& node, Element& element, Grid& grid, GlobalData& globalData)
{
	ifstream file;
	double number = 0.0;
	string word;
	double* data = new double[10];
	file.open(filename);
	if (file.fail())
	{
		cout << "B³¹d pliku" << endl;
	}
	int current_line = 0;
	for (current_line; current_line < 8; current_line++)
	{
		file >> word >> number;
		data[current_line] = number;
	};
	globalData.simulationTime = data[0];
	globalData.simulationStepTime = data[1];
	globalData.conductivity = data[2];
	globalData.alfa = data[3];
	globalData.tot = data[4];
	globalData.initialTemp = data[5];
	globalData.density = data[6];
	globalData.specificHeat = data[7];

	for (current_line; current_line < 10; current_line++)
	{
		file >> word >> word >> number;
		data[current_line] = number;
	};
	grid.nN = data[8];
	grid.nE = data[9];
	grid.ND = new Node[grid.nN];
	grid.El = new Element[grid.nE];

	file >> word;
	for (current_line = 0; current_line < grid.nN; current_line++)
	{
		file >> number >> word >> number >> word;
		grid.ND[current_line].x = number;
		file >> number;
		grid.ND[current_line].y = number;
	}
	file >> word >> word;
	for (current_line = 0; current_line < grid.nE; current_line++)
	{
		file >> number >> word >> number >> word;
		grid.El[current_line].ID[0] = number;
		file >> number >> word;
		grid.El[current_line].ID[1] = number;
		file >> number >> word;
		grid.El[current_line].ID[2] = number;
		file >> number;
		grid.El[current_line].ID[3] = number;
	};
	file >> word;
	int index;
	while (!file.eof())
	{
		file >> index >> word;
		grid.ND[index - 1].BC = 1;
	};
	node = grid.ND[1];
}

//Wywo³ywanie odpowiednich funkcji licz¹cych tabele pochodnych funkcji kszta³tu po ksi i eta
void MatrixEtaKsi(Universal& a, int n) {
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			if (n == 4) {
				a.matrixKsi[i][j] = derivate(1, j + 1, eta4()[i]);
				a.matrixEta[i][j] = derivate(0, j + 1, ksi4()[i]);
			}
			if (n == 9) {
				a.matrixKsi[i][j] = derivate(1, j + 1, eta9()[i]);
				a.matrixEta[i][j] = derivate(0, j + 1, ksi9()[i]);
			}
			if (n == 16) {
				a.matrixKsi[i][j] = derivate(1, j + 1, eta16()[i]);
				a.matrixEta[i][j] = derivate(0, j + 1, ksi16()[i]);
			}
			if (n == 25) {
				a.matrixKsi[i][j] = derivate(1, j + 1, eta25()[i]);
				a.matrixEta[i][j] = derivate(0, j + 1, ksi25()[i]);
			}
		}
	}
}

//Wywo³anie odpowiednich funkcji do obliczenia funkcji kszta³tu po ksi i eta dla poszczególnych punktow calkowania
void EtaKsiFunction(Universal& pkt, int n) {
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < 4; j++) {
			if (n == 4) {
				pkt.matrixEtaKsi[i][j] = calculatingNFunction(j + 1, eta4()[i], ksi4()[i]);
			}
			else if (n == 9) {
				pkt.matrixEtaKsi[i][j] = calculatingNFunction(j + 1, eta9()[i], ksi9()[i]);
			}
			else if (n == 16) {
				pkt.matrixEtaKsi[i][j] = calculatingNFunction(j + 1, eta16()[i], ksi16()[i]);
			}
			else if (n == 25) {
				pkt.matrixEtaKsi[i][j] = calculatingNFunction(j + 1, eta25()[i], ksi25()[i]);
			}
		}
	}
}

int main()
{
	string filename =  "test/Test3_31_31_Kwadrat.txt";//;
	Node node;
	Element element;
	Grid grid;
	GlobalData globalData;
	LoadData(filename, node, element, grid, globalData);

	double pc2[2] = { (-1.0) / sqrt(3) , 1.0 / sqrt(3) };
	double w2[2] = { 1,1 };
	double pc3[3] = { -sqrt(3.0 / 5.0), 0 , sqrt(3.0 / 5.0) };
	double w3[3] = { 5.0 / 9.0,8.0 / 9.0, 5.0 / 9.0 };
	double pc4[4] = { -0.861136, -0.339981, 0.339981, 0.861136 };
	double w4[4] = { 0.347855, 0.652145, 0.652145, 0.347855 };
	double pc5[5] = { -0.906180, -0.538469, 0, 0.538469, 0.906180 };
	double w5[5] = { 0.236927, 0.478629, 0.568889, 0.478629, 0.236927 };

	int a;
	cout << "Podaj schemat calkowania:" << endl;
	cin >> a;
	int pcSide;
	cout << "Podaj ilosc pkt calkowania dla boku:" << endl;
	cin >> pcSide;
	a = a * a;
	Universal pkt(a);
	for (int i = 0; i < a; i++) {
		pkt.matrixEta[i] = new double[4];
		pkt.matrixKsi[i] = new double[4];
		pkt.matrixEtaKsi[i] = new double[4];
	}

	MatrixEtaKsi(pkt, a);  //uzupe³nienie macierzy dN/dKsi oraz dN/dEta
	EtaKsiFunction(pkt, a); //uzupe³nienie macierzy sk³adaj¹cej siê z funkcji kszta³tu

	//---------------------------------------------WYPISYWANIE UNIWERSALNYCH TABLIC-----------------------------------
	for (int i = 0; i < a; i++) {
		for (int j = 0; j < 4; j++) {
			//cout << pkt.matrixEtaKsi[i][j] << " ";
			//cout << pkt.matrixEta[i][j] << " ";
			//cout << pkt.matrixKsi[i][j] << " ";
		}
		cout << endl;
	}

// Przypisanie odpowiednich wag pomno¿onych w zale¿noœci od wybranego schematu ca³kowania
	double* multipleWeights = new double[a];
	for (int i = 0; i < a; i++)
	{
		if (a == 4)
			multipleWeights[i] = multiplyWeights4()[i];
		else if (a == 9)
			multipleWeights[i] = multiplyWeights9()[i];
		else if (a == 16)
			multipleWeights[i] = multiplyWeights16()[i];
		else if (a == 25)
			multipleWeights[i] = multiplyWeights25()[i];
	}

	double*** resultH = new double** [grid.nE];  //tablica do rozwiazan koncowych macierzy H dla ka¿dego elementu
	for (int i = 0; i < grid.nE; i++) {
		resultH[i] = new double* [4];
		for (int j = 0; j < 4; j++) resultH[i][j] = new double[4];
	}
	double*** resultC = new double** [grid.nE]; //tablica do rozwiazan koncowych macierzy C dla ka¿dego elementu
	for (int i = 0; i < grid.nE; i++) {
		resultC[i] = new double* [4];
		for (int j = 0; j < 4; j++) resultC[i][j] = new double[4];
	}

	double*** grids = new double** [grid.nE];   // utworzenie tablicy dla wspolrzednych punktow nE X 2 X 4
	for (int i = 0; i < grid.nE; i++) {
		grids[i] = new double* [2];
		for (int j = 0; j < 2; j++) grids[i][j] = new double[a];
	}
	for (int i = 0; i < grid.nE; i++) {		// wpisywanie do tabeli danych z pliku dotycz¹cych wspolrzednych punktow
		for (int j = 0; j < 4; j++)
		{
			grids[i][0][j] = grid.ND[grid.El[i].ID[j] - 1].x;
			grids[i][1][j] = grid.ND[grid.El[i].ID[j] - 1].y;
		}
	}

	Side* side[4];
	for (int i = 0; i < 4; i++) {
		side[i] = new Side(pcSide);
	}
	for (int i = 0; i < pcSide; i++) {
		side[0]->pcCoordinates[i][1] = -1;
		if (pcSide == 2)
			side[0]->pcCoordinates[i][0] = pc2[i];
		else if (pcSide == 3)
			side[0]->pcCoordinates[i][0] = pc3[i];
		else if (pcSide == 4)
			side[0]->pcCoordinates[i][0] = pc4[i];
		else if (pcSide == 5)
			side[0]->pcCoordinates[i][0] = pc5[i];
	}
	for (int i = 0; i < pcSide; i++) {
		side[1]->pcCoordinates[i][0] = 1;
		if (pcSide == 2)
			side[1]->pcCoordinates[i][1] = pc2[i];
		else if (pcSide == 3)
			side[1]->pcCoordinates[i][1] = pc3[i];
		else if (pcSide == 4)
			side[1]->pcCoordinates[i][1] = pc4[i];
		else if (pcSide == 5)
			side[1]->pcCoordinates[i][1] = pc5[i];
	}
	for (int i = pcSide - 1; i >= 0; i--) {
		side[2]->pcCoordinates[i][1] = 1;
		if (pcSide == 2)
			side[2]->pcCoordinates[pcSide - 1 - i][0] = pc2[i];
		else if (pcSide == 3)
			side[2]->pcCoordinates[pcSide - 1 - i][0] = pc3[i];
		else if (pcSide == 4)
			side[2]->pcCoordinates[pcSide - 1 - i][0] = pc4[i];
		else if (pcSide == 5)
			side[2]->pcCoordinates[pcSide - 1 - i][0] = pc5[i];
	}
	for (int i = pcSide - 1; i >= 0; i--) {
		side[3]->pcCoordinates[i][0] = -1;
		if (pcSide == 2)
			side[3]->pcCoordinates[i][1] = pc2[i];
		else if (pcSide == 3)
			side[3]->pcCoordinates[i][1] = pc3[i];
		else if (pcSide == 4)
			side[3]->pcCoordinates[i][1] = pc4[i];
		else if (pcSide == 5)
			side[3]->pcCoordinates[i][1] = pc5[i];
	}
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < pcSide; j++) {
			side[i]->resultsN[j][0] = 1.0 / 4.0 * (1 - side[i]->pcCoordinates[j][0]) * (1 - side[i]->pcCoordinates[j][1]);
			side[i]->resultsN[j][1] = 1.0 / 4.0 * (1 + side[i]->pcCoordinates[j][0]) * (1 - side[i]->pcCoordinates[j][1]);
			side[i]->resultsN[j][2] = 1.0 / 4.0 * (1 + side[i]->pcCoordinates[j][0]) * (1 + side[i]->pcCoordinates[j][1]);
			side[i]->resultsN[j][3] = 1.0 / 4.0 * (1 - side[i]->pcCoordinates[j][0]) * (1 + side[i]->pcCoordinates[j][1]);
		}
	}
	//----------------------WYPISZ N W BOKACH----------------------------------------
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < pcSide; j++) {
			cout << side[i]->resultsN[j][0] << "\t" << side[i]->resultsN[j][1] << "\t" << side[i]->resultsN[j][2] << "\t" << side[i]->resultsN[j][3] << endl;
		}
	}

	// Utworzenie i wype³nieniee obliczenia wartosci funkcji kszta³ktu dla ka¿dego z punktow
	double*** UniversalEtaKsiForNodes = new double** [grid.nE];
	for (int i = 0; i < grid.nE; i++) {
		UniversalEtaKsiForNodes[i] = new double* [a];
		for (int j = 0; j < a; j++)
			UniversalEtaKsiForNodes[i][j] = new double[4];
	}
	for (int i = 0; i < grid.nE; i++) {
		for (int j = 0; j < a; j++)
			UniversalEtaKsiForNodes[i][j] = pkt.matrixEtaKsi[j];
	}

	for (int i = 0; i < grid.nE; i++) {
		double*** mPc = new double** [a];
		for (int i = 0; i < a; i++) {   // tablica l. pkt ca³kowania x 2 x 2   - przechowywanie wymno¿onych wspolrzednych z matrixKsi i Eta, potem jakobian
			mPc[i] = new double* [2];
			for (int j = 0; j < 2; j++) mPc[i][j] = new double[2];
		}
		mPc = matrixJacobi(grids[i], pkt.matrixEta, pkt.matrixKsi, a); // macierz odwotna bez 1/det
		mPc = jacobian(mPc, a);  // wymno¿ony przez wyznacznik -jakobian i  jakobian odwrotny

		/*for (int i = 0; i < a; i++) {
			for (int j = 0; j < 2; j++) {
				for (int k = 0; k < 2; k++)
					cout << mPc[i][j][k] << " ";
				cout << endl;
			}
			cout << endl;
		}*/
		double*** hbc = new double** [4];
		for (int i = 0; i < 4; i++) {
			hbc[i] = new double* [4];
			for (int j = 0; j < 4; j++) hbc[i][j] = new double[4];
		}
		for (int i = 0; i < 4; i++) {
			for (int j = 0; j < 4; j++) {
				for (int k = 0; k < 4; k++) {
					hbc[i][j][k] = 0;
				}
			}
		}
		//uzupelnienie macierzy hbc dla elementu - kazdego z bokow
		for (int ii = 0; ii < pcSide; ii++) {
			double det1 = sqrt(pow((grid.ND[grid.El[i].ID[0] - 1].x - grid.ND[grid.El[i].ID[1] - 1].x), 2) + pow((grid.ND[grid.El[i].ID[0] - 1].y - grid.ND[grid.El[i].ID[1] - 1].y), 2)) / 2;
			double det2 = sqrt(pow((grid.ND[grid.El[i].ID[1] - 1].x - grid.ND[grid.El[i].ID[2] - 1].x), 2) + pow((grid.ND[grid.El[i].ID[1] - 1].y - grid.ND[grid.El[i].ID[2] - 1].y), 2)) / 2;
			double det3 = sqrt(pow((grid.ND[grid.El[i].ID[2] - 1].x - grid.ND[grid.El[i].ID[3] - 1].x), 2) + pow((grid.ND[grid.El[i].ID[2] - 1].y - grid.ND[grid.El[i].ID[3] - 1].y), 2)) / 2;
			double det4 = sqrt(pow((grid.ND[grid.El[i].ID[0] - 1].x - grid.ND[grid.El[i].ID[3] - 1].x), 2) + pow((grid.ND[grid.El[i].ID[0] - 1].y - grid.ND[grid.El[i].ID[3] - 1].y), 2)) / 2;
			for (int j = 0; j < 4; j++) {
				for (int k = 0; k < 4; k++) {
					if (pcSide == 2) {
						hbc[0][j][k] += side[0]->resultsN[ii][j] * side[0]->resultsN[ii][k] * globalData.alfa * w2[ii] * det1;
						hbc[1][j][k] += side[1]->resultsN[ii][j] * side[1]->resultsN[ii][k] * globalData.alfa * w2[ii] * det2;
						hbc[2][j][k] += side[2]->resultsN[ii][j] * side[2]->resultsN[ii][k] * globalData.alfa * w2[ii] * det3;
						hbc[3][j][k] += side[3]->resultsN[ii][j] * side[3]->resultsN[ii][k] * globalData.alfa * w2[ii] * det4;
					}
					else if (pcSide == 3) {
						if (grid.ND[grid.El[i].ID[0] - 1].BC == 1 && grid.ND[grid.El[i].ID[1] - 1].BC == 1)
						hbc[0][j][k] += side[0]->resultsN[ii][j] * side[0]->resultsN[ii][k] * globalData.alfa * w3[ii] * det1;
						if (grid.ND[grid.El[i].ID[1] - 1].BC == 1 && grid.ND[grid.El[i].ID[2] - 1].BC == 1)
						hbc[1][j][k] += side[1]->resultsN[ii][j] * side[1]->resultsN[ii][k] * globalData.alfa * w3[ii] * det2;
						if (grid.ND[grid.El[i].ID[2] - 1].BC == 1 && grid.ND[grid.El[i].ID[3] - 1].BC == 1)
						hbc[2][j][k] += side[2]->resultsN[ii][j] * side[2]->resultsN[ii][k] * globalData.alfa * w3[ii] * det3;
						if (grid.ND[grid.El[i].ID[0] - 1].BC == 1 && grid.ND[grid.El[i].ID[3] - 1].BC == 1)
						hbc[3][j][k] += side[3]->resultsN[ii][j] * side[3]->resultsN[ii][k] * globalData.alfa * w3[ii] * det4;
					}
					else if (pcSide == 4) {
						hbc[0][j][k] += side[0]->resultsN[ii][j] * side[0]->resultsN[ii][k] * globalData.alfa * w4[ii] * det1;
						hbc[1][j][k] += side[1]->resultsN[ii][j] * side[1]->resultsN[ii][k] * globalData.alfa * w4[ii] * det2;
						hbc[2][j][k] += side[2]->resultsN[ii][j] * side[2]->resultsN[ii][k] * globalData.alfa * w4[ii] * det3;
						hbc[3][j][k] += side[3]->resultsN[ii][j] * side[3]->resultsN[ii][k] * globalData.alfa * w4[ii] * det4;
					}
					else if (pcSide == 5) {
						hbc[0][j][k] += side[0]->resultsN[ii][j] * side[0]->resultsN[ii][k] * globalData.alfa * w5[ii] * det1;
						hbc[1][j][k] += side[1]->resultsN[ii][j] * side[1]->resultsN[ii][k] * globalData.alfa * w5[ii] * det2;
						hbc[2][j][k] += side[2]->resultsN[ii][j] * side[2]->resultsN[ii][k] * globalData.alfa * w5[ii] * det3;
						hbc[3][j][k] += side[3]->resultsN[ii][j] * side[3]->resultsN[ii][k] * globalData.alfa * w5[ii] * det4;
					}
				}
			}
		}
		//------------------------------WYPISYWANIE HBC DLA KA¯DEGO BOKU I ELEMENTU ----------------------------------------
	/*	for (int i = 0; i < 4; i++) {
			for (int j = 0; j < 4; j++) {
				for (int k = 0; k < 4; k++) {
					cout << hbc[i][j][k] << "\t";
				}
				cout << endl;
			}
			cout << endl << endl;
		}*/

		//obliczanie macierzy H i C
		resultH[i] = finalMatrix(hPc(matrixJacobi(grids[i], pkt.matrixEta, pkt.matrixKsi, a), dNdX(mPc, pkt.matrixEta, pkt.matrixKsi, a), dNdY(mPc, pkt.matrixEta, pkt.matrixKsi, a), a, globalData.conductivity), multipleWeights, a);
		resultC[i] = finalMatrix(MatrixCInPc(matrixJacobi(grids[i], pkt.matrixEta, pkt.matrixKsi, a), UniversalEtaKsiForNodes[i], a, globalData.specificHeat, globalData.density), multipleWeights, a);

		//dodawanie H i HBC jesli spelnione warunki
		if (grid.ND[grid.El[i].ID[0] - 1].BC == 1 && grid.ND[grid.El[i].ID[1] - 1].BC == 1)
		{
			for (int j = 0; j < 4; j++) {
				for (int k = 0; k < 4; k++) {
					resultH[i][j][k] += hbc[0][j][k];
				}
			}
		}
		if (grid.ND[grid.El[i].ID[1] - 1].BC == 1 && grid.ND[grid.El[i].ID[2] - 1].BC == 1)
		{
			for (int j = 0; j < 4; j++) {
				for (int k = 0; k < 4; k++) {
					resultH[i][j][k] += hbc[1][j][k];
				}
			}
		}
		if (grid.ND[grid.El[i].ID[2] - 1].BC == 1 && grid.ND[grid.El[i].ID[3] - 1].BC == 1)
		{
			for (int j = 0; j < 4; j++) {
				for (int k = 0; k < 4; k++) {
					resultH[i][j][k] += hbc[2][j][k];
				}
			}
		}
		if (grid.ND[grid.El[i].ID[0] - 1].BC == 1 && grid.ND[grid.El[i].ID[3] - 1].BC == 1)
		{
			for (int j = 0; j < 4; j++) {
				for (int k = 0; k < 4; k++) {
					resultH[i][j][k] += hbc[3][j][k];
				}
			}
		}
	}
	//-----------------------------------WYPISZ POSZCZEGOLNE MACIERZE H + HBC dla ka¿dego node'a lub macierz C--------------------------------
	//cout.precision(15);
	//for (int i = 0; i < grid.nE; i++) {
	//	for (int j = 0; j < 4; j++) {
	//		for (int k = 0; k < 4; k++) {
	//			//cout << resultH[i][j][k] << " ";
	//			//cout << resultC[i][j][k] << " ";
	//		}
	//		cout << "\t" << endl;
	//	}
	//	cout << "\n\n\n";
	//}

	//------------------------------------------MACIERZ DO AGREGACJI MACIERZY H Z HBC -----------------------------------------------------------------------------------------
	double** HG = new double* [grid.nN];
	for (int i = 0; i < grid.nN; i++)
		HG[i] = new double[grid.nN];

	for (int j = 0; j < grid.nN; j++) {
		for (int k = 0; k < grid.nN; k++) {
			HG[j][k] = 0;
		}
	}
	//------------------------------------------MACIERZ DO AGREGACJI MACIERZY C-----------------------------------------------------------------------------------------
	double** CG = new double* [grid.nN];
	for (int i = 0; i < grid.nN; i++)
		CG[i] = new double[grid.nN];

	for (int j = 0; j < grid.nN; j++) {
		for (int k = 0; k < grid.nN; k++) {
			CG[j][k] = 0;
		}
	}

	//---------------------------------------------  AGREGACJA OBU MACIERZY-----------------------------------------------------------------
	for (int i = 0; i < grid.nE; i++) {
		for (int j = 0; j < 4; j++) {
			for (int k = 0; k < 4; k++) {
				HG[grid.El[i].ID[j] - 1][grid.El[i].ID[k] - 1] += resultH[i][j][k];
				CG[grid.El[i].ID[j] - 1][grid.El[i].ID[k] - 1] += resultC[i][j][k];
			}
		}
	}

	//------------------------------WYPISZ AGREGACJE MACIERZY H+ HBC LUB C--------------------------------------------
	//cout.precision(9);
	//for (int i = 0; i < grid.nN; i++) {
	//	for (int j = 0; j < grid.nN; j++) {
	//		cout << HG[i][j] << " ";
	//		//cout << CG[i][j] << " ";
	//	}
	//	cout << endl;
	//}

	double*** p = new double** [grid.nE];
	for (int i = 0; i < grid.nE; i++) {
		p[i] = new double* [4];
		for (int j = 0; j < 4; j++) p[i][j] = new double[4];
	}
	for (int i = 0; i < grid.nE; i++) {
		for (int j = 0; j < 4; j++) {
			for (int k = 0; k < 4; k++) {
				p[i][j][k] = 0;
			}
		}
	}
	for (int i = 0; i < grid.nE; i++) {
		for (int ii = 0; ii < pcSide; ii++) {
			double det1 = sqrt(pow((grid.ND[grid.El[i].ID[0] - 1].x - grid.ND[grid.El[i].ID[1] - 1].x), 2) + pow((grid.ND[grid.El[i].ID[0] - 1].y - grid.ND[grid.El[i].ID[1] - 1].y), 2)) / 2;
			double det2 = sqrt(pow((grid.ND[grid.El[i].ID[1] - 1].x - grid.ND[grid.El[i].ID[2] - 1].x), 2) + pow((grid.ND[grid.El[i].ID[1] - 1].y - grid.ND[grid.El[i].ID[2] - 1].y), 2)) / 2;
			double det3 = sqrt(pow((grid.ND[grid.El[i].ID[2] - 1].x - grid.ND[grid.El[i].ID[3] - 1].x), 2) + pow((grid.ND[grid.El[i].ID[2] - 1].y - grid.ND[grid.El[i].ID[3] - 1].y), 2)) / 2;
			double det4 = sqrt(pow((grid.ND[grid.El[i].ID[0] - 1].x - grid.ND[grid.El[i].ID[3] - 1].x), 2) + pow((grid.ND[grid.El[i].ID[0] - 1].y - grid.ND[grid.El[i].ID[3] - 1].y), 2)) / 2;
			for (int k = 0; k < 4; k++) {
				if (pcSide == 2) {
					if (grid.ND[grid.El[i].ID[0] - 1].BC == 1 && grid.ND[grid.El[i].ID[1] - 1].BC == 1)
						p[i][0][k] += side[0]->resultsN[ii][k] * globalData.alfa * globalData.tot * w2[ii] * det1;
					if (grid.ND[grid.El[i].ID[1] - 1].BC == 1 && grid.ND[grid.El[i].ID[2] - 1].BC == 1)
						p[i][1][k] += side[1]->resultsN[ii][k] * globalData.alfa * globalData.tot * w2[ii] * det2;
					if (grid.ND[grid.El[i].ID[2] - 1].BC == 1 && grid.ND[grid.El[i].ID[3] - 1].BC == 1)
						p[i][2][k] += side[2]->resultsN[ii][k] * globalData.alfa * globalData.tot * w2[ii] * det3;
					if (grid.ND[grid.El[i].ID[0] - 1].BC == 1 && grid.ND[grid.El[i].ID[3] - 1].BC == 1)
						p[i][3][k] += side[3]->resultsN[ii][k] * globalData.alfa * globalData.tot * w2[ii] * det4;
				}
				else if (pcSide == 3) {
					if (grid.ND[grid.El[i].ID[0] - 1].BC == 1 && grid.ND[grid.El[i].ID[1] - 1].BC == 1)
						p[i][0][k] += side[0]->resultsN[ii][k] * globalData.alfa * globalData.tot * w3[ii] * det1;
					if (grid.ND[grid.El[i].ID[1] - 1].BC == 1 && grid.ND[grid.El[i].ID[2] - 1].BC == 1)
						p[i][1][k] += side[1]->resultsN[ii][k] * globalData.alfa * globalData.tot * w3[ii] * det2;
					if (grid.ND[grid.El[i].ID[2] - 1].BC == 1 && grid.ND[grid.El[i].ID[3] - 1].BC == 1)
						p[i][2][k] += side[2]->resultsN[ii][k] * globalData.alfa * globalData.tot * w3[ii] * det3;
					if (grid.ND[grid.El[i].ID[0] - 1].BC == 1 && grid.ND[grid.El[i].ID[3] - 1].BC == 1)
						p[i][3][k] += side[3]->resultsN[ii][k] * globalData.alfa * globalData.tot * w3[ii] * det4;
				}
				else if (pcSide == 4) {
					if (grid.ND[grid.El[i].ID[0] - 1].BC == 1 && grid.ND[grid.El[i].ID[1] - 1].BC == 1)
						p[i][0][k] += side[0]->resultsN[ii][k] * globalData.alfa * globalData.tot * w4[ii] * det1;
					if (grid.ND[grid.El[i].ID[1] - 1].BC == 1 && grid.ND[grid.El[i].ID[2] - 1].BC == 1)
						p[i][1][k] += side[1]->resultsN[ii][k] * globalData.alfa * globalData.tot * w4[ii] * det2;
					if (grid.ND[grid.El[i].ID[2] - 1].BC == 1 && grid.ND[grid.El[i].ID[3] - 1].BC == 1)
						p[i][2][k] += side[2]->resultsN[ii][k] * globalData.alfa * globalData.tot * w4[ii] * det3;
					if (grid.ND[grid.El[i].ID[0] - 1].BC == 1 && grid.ND[grid.El[i].ID[3] - 1].BC == 1)
						p[i][3][k] += side[3]->resultsN[ii][k] * globalData.alfa * globalData.tot * w4[ii] * det4;
				}
				else if (pcSide == 5) {
					if (grid.ND[grid.El[i].ID[0] - 1].BC == 1 && grid.ND[grid.El[i].ID[1] - 1].BC == 1)
						p[i][0][k] += side[0]->resultsN[ii][k] * globalData.alfa * globalData.tot * w5[ii] * det1;
					if (grid.ND[grid.El[i].ID[1] - 1].BC == 1 && grid.ND[grid.El[i].ID[2] - 1].BC == 1)
						p[i][1][k] += side[1]->resultsN[ii][k] * globalData.alfa * globalData.tot * w5[ii] * det2;
					if (grid.ND[grid.El[i].ID[2] - 1].BC == 1 && grid.ND[grid.El[i].ID[3] - 1].BC == 1)
						p[i][2][k] += side[2]->resultsN[ii][k] * globalData.alfa * globalData.tot * w5[ii] * det3;
					if (grid.ND[grid.El[i].ID[0] - 1].BC == 1 && grid.ND[grid.El[i].ID[3] - 1].BC == 1)
						p[i][3][k] += side[3]->resultsN[ii][k] * globalData.alfa * globalData.tot * w5[ii] * det4;
				}
			}
		}
	}
	// ---------------------JEDEN WIERSZ TO JEST P DLA BOKU JEDNEGO - WYPISANIE P DLA KAZDEGO BOKU I  -----------------------------
	//for (int i = 0; i < grid.nE; i++) {
	//	for (int j = 0; j < 4; j++) {
	//		for (int k = 0; k < 4; k++) {
	//			cout << p[i][j][k] << "\t";
	//		}
	//		cout << endl;
	//	}
	//	cout << endl << endl;
	//}

	//----------------------------WEKTORY P DLA KAZDEGO ELEMENTU-------------------------------------------

	for (int i = 0; i < grid.nE; i++) {
		for (int j = 0; j < 4; j++) {
			for (int k = 1; k < 4; k++) {
				p[i][0][j] += p[i][k][j];
			}
		}
	}
	//--------------------WYPISANIE WEKTORA P DLA ELEMENTOW ------------------------------------------------
	//for (int i = 0; i < grid.nE; i++) {
	//	for (int j = 0; j < 4; j++) {
	//		cout << p[i][0][j] << "\t";
	//	}
	//	cout << endl << endl;
	//}
	//--------------------------------AGREGACJA WEKTORA P ------------------------------------------------------
	double* Pagreg = new double[grid.nN];
	for (int i = 0; i < grid.nN; i++) {
		Pagreg[i] = 0;
	}

	for (int i = 0; i < grid.nE; i++) {
		for (int k = 0; k < 4; k++) {
			Pagreg[grid.El[i].ID[k] - 1] += p[i][0][k];
		}
	}
	//-------------------------------------------------Wypisywanie P zagregowanego----------------------------------
	for (int i = 0; i < grid.nN; i++) {
		//cout << Pagreg[i] << endl;
	}

	// -----------------------------------------------OBLICZENIE METOD¥ ELIMINACJI GAUSSA JEDNO T-----------------------------------
	/*cout << endl;
	double* T = new double[grid.nN];
	T = GaussElimination(HG, Pagreg, grid.nN);
	for (int i = 0; i < grid.nN; i++) cout << T[i] << endl;*/

	// --------------------------------PRZYGOTOWANIE DO ROZWIAZYWANIA UKLADU ROWNAN-------------------------------------
	double* T0 = new double[grid.nN];
	for (int i = 0; i < grid.nN; i++) {
		T0[i] = globalData.initialTemp;
	}
	double** SumHAndC = new double* [grid.nN];
	for (int i = 0; i < grid.nN; i++) {
		SumHAndC[i] = new double[grid.nN];
	}
	for (int i = 0; i < grid.nN; i++) {
		for (int j = 0; j < grid.nN; j++)
			SumHAndC[i][j] = HG[i][j] + (CG[i][j] / globalData.simulationStepTime);
	}
	//-------------------------------WYPISANIE SUMY MACIERZY H I C------------------------------------------------------------
	//cout.precision(9);
	//for (int i = 0; i < grid.nN; i++) {
	//	for (int j = 0; j < grid.nN; j++) {
	//		cout << SumHAndC[i][j] << " ";
	//	}
	//	cout << endl;
	//}

	double* SumCAndP = new double[grid.nN];
	double minTemp = 0.0;
	double maxTemp = 0.0;
	//---------------------------- TEMPERATURY MAX I MIN W CZASIE CALEJ SYMULACJI-----------------------------
	cout.precision(6);
	cout << "Time[s] MinTemp[s]  MaxTemp[s]" << endl;
	for (int step = 0; step < globalData.simulationTime; step += globalData.simulationStepTime) {
		for (int i = 0; i < grid.nN; i++) {
			SumCAndP[i] = 0;
			for (int j = 0; j < grid.nN; j++) {
				SumCAndP[i] += CG[i][j] / globalData.simulationStepTime * T0[j];
			}
			SumCAndP[i] += Pagreg[i];
		}

		/*for (int i = 0; i < grid.nN; i++ ) {
			cout << SumCAndP[i] << endl;
		}
		cout << endl;*/
		T0 = GaussElimination(SumHAndC, SumCAndP, grid.nN);
		/*for (int i = 0; i < grid.nN; i++) cout << T0[i] << endl;
		cout << endl;*/
		minTemp = findMin(T0, grid.nN);
		maxTemp = findMax(T0, grid.nN);
		cout << step + globalData.simulationStepTime << "      " << minTemp << "       " << maxTemp << endl;
	}

	return 0;
}