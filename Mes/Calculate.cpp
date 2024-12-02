#include "Calculate.h"
#include "Grid.h"
#include "Element.h"
#include "Node.h"
#include "Jakobian.h"
#include "MatrixH.h"
#include "MatrixHbc.h"

void calculate(int iloscPc, Grid* grid, GlobalData* globaldata, Node* node) 
{
	Element element;
    for (int i = 0; i < globaldata->nE; i++) {                          
        element.addPointsX(                       
            grid->node[grid->element[i].ID[0] - 1].x,
            grid->node[grid->element[i].ID[1] - 1].x,
            grid->node[grid->element[i].ID[2] - 1].x,
            grid->node[grid->element[i].ID[3] - 1].x
        );
        element.addPointsY(
            grid->node[grid->element[i].ID[0] - 1].y,
            grid->node[grid->element[i].ID[1] - 1].y,
            grid->node[grid->element[i].ID[2] - 1].y,
            grid->node[grid->element[i].ID[3] - 1].y
        );
        element.obecne[0] = grid->element[i].ID[0];
        element.obecne[1] = grid->element[i].ID[1];
        element.obecne[2] = grid->element[i].ID[2];
        element.obecne[3] = grid->element[i].ID[3];

        startCalculating(iloscPc,&element, node, globaldata, grid);
    }
}

void startCalculating(int punktyCalkowania, Element* element, Node* node, GlobalData* globaldata, Grid* grid)
{
    vector<double> ksi;
    vector<double> eta;
    double k = globaldata->Conductivity;

    //kolejnosc znakow etc dla punktow
    if (punktyCalkowania == 4) {
        ksi = { node->points[0], node->points[1], node->points[1], node->points[0] };
        eta = { node->points[0], node->points[0], node->points[1], node->points[1] };
    }
    else if (punktyCalkowania == 9) {
        ksi = { node->points[0], node->points[1], node->points[2], node->points[2], node->points[1], node->points[0], node->points[0], node->points[1], node->points[2] };
        eta = { node->points[0], node->points[0], node->points[0], node->points[1], node->points[1], node->points[1], node->points[2], node->points[2], node->points[2] };
    }
    else if (punktyCalkowania == 16) {
        ksi = { node->points[0], node->points[1], node->points[2], node->points[3],
        node->points[0], node->points[1], node->points[2], node->points[3],
        node->points[0], node->points[1], node->points[2], node->points[3],
        node->points[0], node->points[1], node->points[2], node->points[3] };

        eta = { node->points[0], node->points[0], node->points[0], node->points[0],
                node->points[1], node->points[1], node->points[1], node->points[1],
                node->points[2], node->points[2], node->points[2], node->points[2],
                node->points[3], node->points[3], node->points[3], node->points[3] };
    }
    else {
        std::cerr << "Nieprawid³owa liczba punktów calkowania!" << std::endl;
        return;
    }

    for (int i = 0; i < punktyCalkowania; i++) {
        // KSI
        element->Ksi[i][0] = -0.25 * (1 - eta[i]);
        element->Ksi[i][1] = 0.25 * (1 - eta[i]);
        element->Ksi[i][2] = 0.25 * (1 + eta[i]);
        element->Ksi[i][3] = -0.25 * (1 + eta[i]);

        // ETA
        element->Eta[i][0] = -0.25 * (1 - ksi[i]);
        element->Eta[i][1] = -0.25 * (1 + ksi[i]);
        element->Eta[i][2] = 0.25 * (1 + ksi[i]);
        element->Eta[i][3] = 0.25 * (1 - ksi[i]);

    }
    //tu dodaje funkcje ktore dzialaja na obecnych node'ach
    Jakobian* jakobian = new Jakobian;
    calculateMatrixJakobiego(punktyCalkowania,element, jakobian);
    calculateMatrixH(punktyCalkowania, element, jakobian, node, k);
    calculateMatrixHbc(punktyCalkowania, element, node, grid, ksi, eta);
}

void calculateMatrixJakobiego(int punktyCalkowania, Element* element, Jakobian* jakobian)
{
    for (int i = 0; i < punktyCalkowania; i++) {
        jakobian->calculateMatrix(i,element);
        jakobian->calculateDetJ();
        jakobian->calculateMatrixReversed();
    }
}

void calculateMatrixH(int punktyCalkowania, Element* element, Jakobian* jakobian, Node* node, double k)
{
    vector<vector<double>> dNdx(punktyCalkowania, vector<double>(4));           //Zbiory na pochodne funkcji ksztaltu
    vector<vector<double>> dNdy(punktyCalkowania, vector<double>(4));          //po x i y
    MatrixH matrixH;
    matrixH.calculateDerivetivedN(punktyCalkowania, dNdx, dNdy, element, jakobian);
    matrixH.calculateH(punktyCalkowania, dNdx, dNdy, jakobian, node, k, element);
}

void calculateMatrixHbc(int, Element* element, Node* node, Grid* grid, vector<double> ksi, vector<double> eta)
{
    MatrixHbc matrixHbc;
    matrixHbc.calculateN();
}
