#pragma once
#include "Element.h"

struct MatrixHbc {
    double Hb_czastkowe[4][4] = { 0.0 };
    double Hb_calkowite[4][4] = { 0.0 };
    double detJ;


    void calculateN(int punktyCalkowania, Element* element, Grid* grid, vector<double> ksiBc, vector<double> etaBc, double alfa, Node* node) {
        detJ = 0.016666666;

        memset(Hb_czastkowe, 0, sizeof(Hb_czastkowe));
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

            for (int punkt = 0; punkt < sqrt(punktyCalkowania); ++punkt) {
                double N1 = 0.25 * (1 - ksiBc[krawedz * 2]) * (1 - etaBc[krawedz * 2]);
                double N2 = 0.25 * (1 + ksiBc[krawedz * 2]) * (1 - etaBc[krawedz * 2]);
                double N3 = 0.25 * (1 + ksiBc[krawedz * 2]) * (1 + etaBc[krawedz * 2]);
                double N4 = 0.25 * (1 - ksiBc[krawedz * 2]) * (1 + etaBc[krawedz * 2]);

                double N[4] = { N1, N2, N3, N4 };

                for (int i = 0; i < 4; ++i) {
                    for (int j = 0; j < 4; ++j) {
                        Hb_calkowite[i][j] += alfa * N[i] * N[j] * detJ;
                    }
                }
            }
        }

        cout << "Macierz Hbc:" << endl;
        for (int i = 0; i < 4; ++i) {
            for (int j = 0; j < 4; ++j) {
                cout << Hb_calkowite[i][j] << " ";
            }
            cout << endl;
        }
        cout << "\n\n";

        for (int i = 0; i < 4; ++i) {
            for (int j = 0; j < 4; ++j) {
                node->H_GLOBAL[element->obecne[i] - 1][element->obecne[j] - 1] += Hb_calkowite[i][j];
            }
        }
    }
};