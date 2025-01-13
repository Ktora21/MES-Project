struct Node {
	long double x, y; 
	int BC = 0;
};

struct Element {
	int* ID = new int[4];
};

struct Grid {
	int nE;
	int nN;
	Node* ND;
	Element* El;
};

struct GlobalData {
	double simulationTime, simulationStepTime, conductivity, alfa, tot, initialTemp, density, specificHeat;
};

struct Universal {
	double** matrixKsi;
	double** matrixEta;
	double** matrixEtaKsi;

	Universal(int n) {
		matrixKsi = new double* [n];
		matrixEta = new double* [n];
		matrixEtaKsi = new double* [n];
	}
};

struct Side {
	int pc;
	double** pcCoordinates;
	double** resultsN;

	Side() {}
	Side(int pc) {
		this->pc = pc;
		pcCoordinates = new double* [pc];
		for (int i = 0; i < pc; i++)
			pcCoordinates[i] = new double[2];
		resultsN = new double* [pc];
		for (int i = 0; i < pc; i++)
			resultsN[i] = new double[4];
	}
};