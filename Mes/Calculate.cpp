#include "Calculate.h"
#include "Grid.h"
#include "Element.h"
#include "Node.h"
#include "Jakobian.h"
#include "MatrixH.h"
#include "MatrixHbc.h"

void calculate(int punktyCalkowania, Element* element, Grid* grid, GlobalData* globaldata, Node* node)
{
    Jakobian* jakobian = new Jakobian;
    // ustawia ksi i eta w zaleznosci od punktow calkowania
    grid->settings(punktyCalkowania, node);

    // ustawia wartosci ksi i eta
    calculateKsiEta(punktyCalkowania, element, node, globaldata, grid);

    // ten fragment wykonuje sie dla kazdego elementu
    for (int i = 0; i < globaldata->nE; i++) {    

        element->addPointsX(                       
            grid->node[grid->element[i].ID[0] - 1].x,
            grid->node[grid->element[i].ID[1] - 1].x,
            grid->node[grid->element[i].ID[2] - 1].x,
            grid->node[grid->element[i].ID[3] - 1].x
        );
        element->addPointsY(
            grid->node[grid->element[i].ID[0] - 1].y,
            grid->node[grid->element[i].ID[1] - 1].y,
            grid->node[grid->element[i].ID[2] - 1].y,
            grid->node[grid->element[i].ID[3] - 1].y
        );
        element->obecne[0] = grid->element[i].ID[0];
        element->obecne[1] = grid->element[i].ID[1];
        element->obecne[2] = grid->element[i].ID[2];
        element->obecne[3] = grid->element[i].ID[3];
        // tutaj dodaje funckje
        calculateMatrixH(punktyCalkowania, element, jakobian, node, globaldata->Conductivity);
        calculateMatrixHbc(punktyCalkowania, element, node, grid, globaldata);
    }
    node->printGlobalH();
    node->printGlobalP();
}

void calculateKsiEta(int punktyCalkowania, Element* element, Node* node, GlobalData* globaldata, Grid* grid)
{
    for (int i = 0; i < punktyCalkowania; i++) {
        element->Ksi[i][0] = -0.25 * (1 - grid->eta[i]);
        element->Ksi[i][1] = 0.25 * (1 - grid->eta[i]);
        element->Ksi[i][2] = 0.25 * (1 + grid->eta[i]);
        element->Ksi[i][3] = -0.25 * (1 + grid->eta[i]);

        element->Eta[i][0] = -0.25 * (1 - grid->ksi[i]);
        element->Eta[i][1] = -0.25 * (1 + grid->ksi[i]);
        element->Eta[i][2] = 0.25 * (1 + grid->ksi[i]);
        element->Eta[i][3] = 0.25 * (1 - grid->ksi[i]);
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

void calculateMatrixHbc(int punktyCalkowania, Element* element, Node* node, Grid* grid, GlobalData* globaldata)
{
    MatrixHbc matrixHbc;
    matrixHbc.calculateN(punktyCalkowania, element, grid, globaldata, node);
}
