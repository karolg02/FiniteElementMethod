#include "MatrixH.h"

void MatrixH::calculateDerivetivedN(int punktyCalkowania, vector<vector<double>>& dNdx, vector<vector<double>>& dNdy, Element* element, Jakobian* jakobian)
{
    for (int i = 0; i < punktyCalkowania; i++) {
        for (int j = 0; j < 4; j++) {
            //dla x
            dNdx[i][j] = jakobian->macierzOdwrotna[0] * element->Ksi[i][j] + jakobian->macierzOdwrotna[1] * element->Eta[i][j];
            // dla y
            dNdy[i][j] = jakobian->macierzOdwrotna[2] * element->Ksi[i][j] + jakobian->macierzOdwrotna[3] * element->Eta[i][j];
        }
    }
}

void MatrixH::calculateH(int punktyCalkowania, vector<vector<double>>& dNdx, vector<vector<double>>& dNdy, Jakobian* jakobian, Node* node, double k, Element* element)
{
    for (int i = 0; i < punktyCalkowania; i++) {
        for (int j = 0; j < 4; j++) {
            Hleft[j][0] = dNdx[i][j] * dNdx[i][0];
            Hleft[j][1] = dNdx[i][j] * dNdx[i][1];
            Hleft[j][2] = dNdx[i][j] * dNdx[i][2];
            Hleft[j][3] = dNdx[i][j] * dNdx[i][3];
        }
        //right
        for (int j = 0; j < 4; j++) {
            Hright[j][0] = dNdy[i][j] * dNdy[i][0];
            Hright[j][1] = dNdy[i][j] * dNdy[i][1];
            Hright[j][2] = dNdy[i][j] * dNdy[i][2];
            Hright[j][3] = dNdy[i][j] * dNdy[i][3];
        }
        if (punktyCalkowania == 4) {
            node->weightsPathX = { 0, 1, 1, 1 };
            node->weightsPathY = { 0, 0, 1, 1 };

            for (int j = 0; j < 4; j++) {
                int x = node->weightsPathX[i];
                int y = node->weightsPathY[i];

                H[j][0] = k * jakobian->detJ * (Hleft[j][0] + Hright[j][0]) * node->weights[x] * node->weights[y];
                H[j][1] = k * jakobian->detJ * (Hleft[j][1] + Hright[j][1]) * node->weights[x] * node->weights[y];
                H[j][2] = k * jakobian->detJ * (Hleft[j][2] + Hright[j][2]) * node->weights[x] * node->weights[y];
                H[j][3] = k * jakobian->detJ * (Hleft[j][3] + Hright[j][3]) * node->weights[x] * node->weights[y];
            }
        }
        else if (punktyCalkowania == 9) {
            node->weightsPathX = { 0, 1, 2, 2, 1, 0, 0, 1, 2 };
            node->weightsPathY = { 0, 0, 0, 1, 1, 1, 2, 2, 2 };

            for (int j = 0; j < 4; j++) {
                int x = node->weightsPathX[i];
                int y = node->weightsPathY[i];

                H[j][0] = k * jakobian->detJ * (Hleft[j][0] + Hright[j][0]) * node->weights[x] * node->weights[y];
                H[j][1] = k * jakobian->detJ * (Hleft[j][1] + Hright[j][1]) * node->weights[x] * node->weights[y];
                H[j][2] = k * jakobian->detJ * (Hleft[j][2] + Hright[j][2]) * node->weights[x] * node->weights[y];
                H[j][3] = k * jakobian->detJ * (Hleft[j][3] + Hright[j][3]) * node->weights[x] * node->weights[y];
            }
        }
        else if (punktyCalkowania == 16) {
            node->weightsPathX = { 0, 1, 2, 3,
                                  0, 1, 2, 3,
                                  0, 1, 2, 3,
                                  0, 1, 2, 3 };
            node->weightsPathY = { 0, 0, 0, 0,
                                  1, 1, 1, 1,
                                  2, 2, 2, 2,
                                  3, 3, 3, 3 };

            for (int j = 0; j < 4; j++) {
                int x = node->weightsPathX[i];
                int y = node->weightsPathY[i];

                H[j][0] = k * jakobian->detJ * (Hleft[j][0] + Hright[j][0]) * node->weights[x] * node->weights[y];
                H[j][1] = k * jakobian->detJ * (Hleft[j][1] + Hright[j][1]) * node->weights[x] * node->weights[y];
                H[j][2] = k * jakobian->detJ * (Hleft[j][2] + Hright[j][2]) * node->weights[x] * node->weights[y];
                H[j][3] = k * jakobian->detJ * (Hleft[j][3] + Hright[j][3]) * node->weights[x] * node->weights[y];
            }
        }
        //sumowanie h czastkowych
        for (int j = 0; j < 4; j++) {
            Hc[j][0] = Hc[j][0] + H[j][0];
            Hc[j][1] = Hc[j][1] + H[j][1];
            Hc[j][2] = Hc[j][2] + H[j][2];
            Hc[j][3] = Hc[j][3] + H[j][3];
        }
    }
    //sumowanie do globalnej macierzy H
    for (int j = 0; j < 4; ++j) {
        for (int k = 0; k < 4; ++k) {
            node->H_GLOBAL[element->obecne[j] - 1][element->obecne[k] - 1] += Hc[j][k];
        }
    }
}
