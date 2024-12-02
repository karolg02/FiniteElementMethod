#pragma once
#include "Element.h"
struct MatrixHbc {
	double Hb_czastkowe[4][4] = { 0 };
	double Hb_calkowite[4][4] = { 0 };
	double Hbleft[4][4] = { 0 };
	double Hbright[4][4] = { 0 };
	double detJ;

	void calculateN(int punktyCalkowania,Element* element,Grid* grid, vector<double> ksiBc, vector<double> etaBc, double alfa) {
		detJ = 0.0125;
		alfa = 25;
		for (int i = 0; i < ksiBc.size(); i++) {
			double N1 = 0.25 * (1 - ksiBc[i]) * (1 - etaBc[i]); // Wêze³ 1
			double N2 = 0.25 * (1 + ksiBc[i]) * (1 - etaBc[i]); // Wêze³ 2
			double N3 = 0.25 * (1 + ksiBc[i]) * (1 + etaBc[i]); // Wêze³ 3
			double N4 = 0.25 * (1 - ksiBc[i]) * (1 + etaBc[i]); // Wêze³ 4

			double N[4] = { N1, N2, N3, N4 };

			cout << "Punkt brzegowy: " << i + 1 << endl;
			cout << "ksiBc: " << ksiBc[i] << " etaBc: " << etaBc[i] << endl;
			cout << "N1: " << N1 << " N2: " << N2 << " N3: " << N3 << " N4: " << N4 << endl;
			
			for (int j = 0; j < 4; j++) {
                for (int k = 0; k < 4; k++) {
                    Hb_czastkowe[j][k] = alfa * (1 * N[j] * N[k]) * detJ;
                }
            }

			/*for (int j = 0; j < 4; j++) {
				for (int k = 0; k < 4; k++) {
					cout << Hb_czastkowe[j][k] << ",";
				}
				cout << endl;
			}*/

			// Dodanie wk³adu Hb_czastkowe do Hb_calkowite
			/*for (int j = 0; j < 4; j++) {
				for (int k = 0; k < 4; k++) {
					Hb_calkowite[j][k] += Hb_czastkowe[j][k];
				}
			}*/

		}
		/*cout << "Macierz Hbc:" << endl;
		for (int i = 0; i < 4; i++) {
			for (int j = 0; j < 4; j++) {
				cout << Hb_calkowite[i][j] << " ";
			}
			cout << endl;
		}*/
	}
};