#pragma once
#include "Element.h"
struct MatrixHbc {
	double Hb_czastkowe[4][4];
	double Hb_calkowite[4][4];
	double Hbleft[4][4] = { 0 };
	double Hbright[4][4] = { 0 };
	double N1;
	double N2;
	double N3;
	double N4;

	void calculateN(int punktyCalkowania,Element* element,Grid* grid, vector<double> ksi, vector<double> eta) {
		for (int i = 0; i < punktyCalkowania; i++) {
			for (int j = 0; j < grid->boundaryConditions.size(); j++) {
				if (element->obecne[i] == grid->boundaryConditions[j]) {
					N1 = 0.25 * (1 - ksi[i]) * (1 - eta[i]);
					N2 = 0.25 * (1 + ksi[i]) * (1 - eta[i]);
					N3 = 0.25 * (1 + ksi[i]) * (1 + eta[i]);
					N4 = 0.25 * (1 - ksi[i]) * (1 + eta[i]);
					cout<< "Punkt Bc" << element->obecne[i] << " " << N1 << " " << N2 << " " << N3 << " " << N4 << endl;
				}
			}
		}
	}
};