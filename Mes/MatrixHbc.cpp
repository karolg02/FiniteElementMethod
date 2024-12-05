#include "MatrixHbc.h"

void MatrixHbc::calculateN(int punktyCalkowania, Element* element, Grid* grid, GlobalData* globaldata, Node* node)
{
    // dla kazdej krawedzi, w elemencie
    for (int krawedz = 0; krawedz < 4; ++krawedz) {
        int node1 = element->obecne[krawedz];
        int node2 = element->obecne[(krawedz + 1) % 4];

        bool isBoundaryEdge = false;
        for (int i = 0; i < grid->boundaryConditions.size(); i++) {
            if (node1 == grid->boundaryConditions[i]) {
                for (int j = 0; j < grid->boundaryConditions.size(); j++) {
                    if (node2 == grid->boundaryConditions[j]) {
                        isBoundaryEdge = true;
                        break;
                    }
                }
            }
            if (isBoundaryEdge) break;
        }

        if (!isBoundaryEdge) continue;


        /*cout << "x1= " << element->pointsX[krawedz] << " x2= " << element->pointsX[(krawedz + 1) % 4] << endl;
        cout << "y1= " << element->pointsY[krawedz] << " y2= " << element->pointsY[(krawedz + 1) % 4] << endl;*/
        detJ = sqrt(pow((element->pointsX[krawedz] - element->pointsX[(krawedz + 1) % 4]), 2)
            + pow((element->pointsY[krawedz] - element->pointsY[(krawedz + 1) % 4]), 2)) / 2;

        //cout << "detJ " << detJ << endl;

        int x = sqrt(punktyCalkowania);

        for (int punkt = 0; punkt < x; ++punkt) {

            double N1 = 0.25 * (1 - grid->ksiBc[krawedz * x + punkt]) * (1 - grid->etaBc[krawedz * x + punkt]);
            double N2 = 0.25 * (1 + grid->ksiBc[krawedz * x + punkt]) * (1 - grid->etaBc[krawedz * x + punkt]);
            double N3 = 0.25 * (1 + grid->ksiBc[krawedz * x + punkt]) * (1 + grid->etaBc[krawedz * x + punkt]);
            double N4 = 0.25 * (1 - grid->ksiBc[krawedz * x + punkt]) * (1 + grid->etaBc[krawedz * x + punkt]);

            double N[4] = { N1, N2, N3, N4 };

            //cout << "N1: " << N1 << " " << "N2: " << N2 << " " << "N3: " << N3 << " " << "N4: " << N4 << endl;

            for (int i = 0; i < 4; ++i) {
                P[i] += globaldata->Alfa * globaldata->Tot * N[i] * detJ * node->weights[punkt];
                for (int j = 0; j < 4; ++j) {
                    Hb_calkowite[i][j] += globaldata->Alfa * N[i] * N[j] * detJ * node->weights[punkt];
                }
            }
        }
    }

    cout << "Macierz Hbc:" << endl;
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            cout << Hb_calkowite[i][j] << "  ";
        }
        cout << endl;
    }
    cout << "\n\n";

    cout << "Wektor P:" << endl;
    for (int i = 0; i < 4; ++i) {
        cout << P[i] << " ";
    }
    cout << "\n\n";

    for (int i = 0; i < 4; ++i) {
        node->P_GLOBAL[element->obecne[i] - 1] += P[i];
        for (int j = 0; j < 4; ++j) {
            node->H_GLOBAL[element->obecne[i] - 1][element->obecne[j] - 1] += Hb_calkowite[i][j];
        }
    }
}
