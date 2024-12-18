#include "MatrixC.h"

void MatrixC::calculateC(int punktyCalkowania, Element* element, Node* node, Grid* grid, GlobalData* globaldata, Jakobian* jakobian)
{
	for (int i = 0; i < punktyCalkowania; i++) {
		double N1 = 0.25 * (1 - grid->ksi[i]) * (1 - grid->eta[i]);
		double N2 = 0.25 * (1 + grid->ksi[i]) * (1 - grid->eta[i]);
		double N3 = 0.25 * (1 + grid->ksi[i]) * (1 + grid->eta[i]);
		double N4 = 0.25 * (1 - grid->ksi[i]) * (1 + grid->eta[i]);

		cout << N1 << ", " << N2 << ", " << N3 << ", " << N4 << endl;
		double N[4] = { N1, N2, N3, N4 };

		jakobian->calculateMatrix(i, element);
		jakobian->calculateDetJ();

		for (int j = 0; j < 4; ++j) {
			for (int k = 0; k < 4; ++k) {
				C[j][k] += globaldata->Density* globaldata->SpecificHeat * N[j] * N[k] * jakobian->detJ * node->weights[node->weightsPathX[i]] * node->weights[node->weightsPathY[i]];
			}
		}
	}

	for (int i = 0; i < 4; ++i) {
		for (int j = 0; j < 4; ++j) {
			node->C_GLOBAL[element->obecne[i] - 1][element->obecne[j] - 1] += C[i][j];
		}
	}

}
