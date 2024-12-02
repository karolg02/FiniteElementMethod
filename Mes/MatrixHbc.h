#pragma once
#include "Element.h"

struct MatrixHbc {
    double Hb_czastkowe[4][4] = { 0.0 };
    double Hb_calkowite[4][4] = { 0.0 };
    double detJ;
    

    void calculateN(int punktyCalkowania, Element* element, Grid* grid, vector<double> ksiBc, vector<double> etaBc, double alfa, Node* node) {
        detJ = 0.01666666;
        for (int krawedz = 0; krawedz < 4; ++krawedz) {
            bool isBoundaryEdge = false; // flaga dla granicy
            for (int i = 0; i < grid->boundaryConditions.size(); i++) {
                if (element->obecne[krawedz] == grid->boundaryConditions[i]) {
                    isBoundaryEdge = true;
                    break;
                }
            }

            if (!isBoundaryEdge) {
                continue; // jeœli krawêdŸ nie jest graniczn¹, przechodzimy do nastêpnej krawêdzi
            }

            // iteracja po punktach ca³kowania na krawêdzi
            for (int punkt = 0; punkt < punktyCalkowania; ++punkt) {
                double N1 = 0.25 * (1 - ksiBc[punkt]) * (1 - etaBc[punkt]); // Wêze³ 1
                double N2 = 0.25 * (1 + ksiBc[punkt]) * (1 - etaBc[punkt]); // Wêze³ 2
                double N3 = 0.25 * (1 + ksiBc[punkt]) * (1 + etaBc[punkt]); // Wêze³ 3
                double N4 = 0.25 * (1 - ksiBc[punkt]) * (1 + etaBc[punkt]); // Wêze³ 4

                double N[4] = { N1, N2, N3, N4 };

                for (int i = 0; i < 4; ++i) {
                    for (int j = 0; j < 4; ++j) {
                        Hb_czastkowe[i][j] = alfa * N[i] * N[j] * detJ;
                    }
                }
                for (int i = 0; i < 4; ++i) {
                    for (int j = 0; j < 4; ++j) {
                        Hb_calkowite[i][j] += Hb_czastkowe[i][j];
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

        for (int i = 0; i < 4; ++i) {
            for (int j = 0; j < 4; ++j) {
                node->H_GLOBAL[element->obecne[i] - 1][element->obecne[j] - 1] += Hb_calkowite[i][j];
            }
        }
    }
};
