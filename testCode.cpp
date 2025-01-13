
//double pc2[2] = { (-1.0) / sqrt(3) , 1.0 / sqrt(3) };
//double w2[2] = { 1,1 };
//double pc3[3] = { -sqrt(3.0 / 5.0), 0 , sqrt(3.0 / 5.0) };
//double w3[3] = { 5.0 / 9.0,8.0 / 9.0, 5.0 / 9.0 };
//double pc4[4] = { -0.861136, -0.339981, 0.339981, 0.861136 };
//double w4[4] = { 0.347855, 0.652145, 0.652145, 0.347855 };
//double pc5[5] = { -0.906180, -0.538469, 0, 0.538469, 0.906180 };
//double w5[5] = { 0.236927, 0.478629, 0.568889, 0.478629, 0.236927 };

//cout.precision(15);
//cout << "2 pkt\t\t" << value1(pc2, w2, 2, 3, 8)<< endl;
//cout << "3 pkt\t\t" << value1(pc3, w3, 3, 3, 8) << endl;
//cout << "4 pkt\t\t" << value1(pc4, w4, 4, 3, 8) << endl;
//cout << "5 pkt\t\t" << value1(pc5, w5, 5, 3, 8) << endl;

//cout.precision(15);
//cout << "2 pkt, one variable:\t" << value(pc2, w2, 1, 2) << endl;
//cout << "2 pkt, two variable:\t" << value(pc2, w2, 2, 2) << endl;
//cout << "3 pkt, one variable:\t" << value(pc3, w3, 1, 3) << endl;
//cout << "3 pkt, two variable:\t" << value(pc3, w3, 2, 3) << endl;
//cout << "4 pkt, one variable:\t" << value(pc4, w4, 1, 4) << endl;
//cout << "4 pkt, two variable:\t" << value(pc4, w4, 2, 4) << endl;
//cout << "5 pkt, one variable:\t" << value(pc5, w5, 1, 5) << endl;
//cout << "5 pkt, two variable:\t" << value(pc5, w5, 2, 5) << endl;

//cout.precision(15);
//double** oneGrid = new double* [2];
//for (int i = 0; i < 2; i++) oneGrid[i] = new double[4];
//oneGrid[0][0] = 0;
//oneGrid[0][1] = 0.025;
//oneGrid[0][2] = 0.025;
//oneGrid[0][3] = 0;
//oneGrid[1][0] = 0;
//oneGrid[1][1] = 0;
//oneGrid[1][2] = 0.025;
//oneGrid[1][3] = 0.025;
//
//double*** mPc = new double** [a];
//for (int i = 0; i < a; i++) {
//	mPc[i] = new double* [2];
//	for (int j = 0; j < 2; j++)
//		mPc[i][j] = new double[2];
//}
//mPc = matrixJacobi(oneGrid, pkt.matrixEta, pkt.matrixKsi, a);
//mPc = jacobian(mPc, a);
//cout << " Macierze H PC\n";
//
//print3dMatrix(hPc(matrixJacobi(oneGrid, pkt.matrixEta, pkt.matrixKsi, a), dNdX(mPc, pkt.matrixEta, pkt.matrixKsi, a), dNdY(mPc, pkt.matrixEta, pkt.matrixKsi, a), a, 30), a, 4);
//cout << "\n\n\n\n\n\n";
//if(a == 4)
//	print2dMatrix(finalH(hPc(matrixJacobi(oneGrid, pkt.matrixEta, pkt.matrixKsi, a), dNdX(mPc, pkt.matrixEta, pkt.matrixKsi, a), dNdY(mPc, pkt.matrixEta, pkt.matrixKsi, a), a, 30), multipleWeights, a), a);
//else if(a == 9)
//	print2dMatrix(finalH(hPc(matrixJacobi(oneGrid, pkt.matrixEta, pkt.matrixKsi, a), dNdX(mPc, pkt.matrixEta, pkt.matrixKsi, a), dNdY(mPc, pkt.matrixEta, pkt.matrixKsi, a), a, 30),multipleWeights, a), a);
//else if (a == 16)
//	print2dMatrix(finalH(hPc(matrixJacobi(oneGrid, pkt.matrixEta, pkt.matrixKsi, a), dNdX(mPc, pkt.matrixEta, pkt.matrixKsi, a), dNdY(mPc, pkt.matrixEta, pkt.matrixKsi, a), a, 30), multipleWeights, a), a);
