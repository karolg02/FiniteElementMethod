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
//kod wykonujacy sie dla kazdego elementu
void startCalculating(int punktyCalkowania, Element* element, Node* node, GlobalData* globaldata, Grid* grid)
{
    vector<double> ksi;
    vector<double> eta;
    vector<double> ksiBc;
    vector<double> etaBc;
    double k = globaldata->Conductivity;

    //kolejnosc znakow etc dla punktow
    if (punktyCalkowania == 4) {
        ksi = { node->points[0], node->points[1], node->points[1], node->points[0] }; //points[0] - 0,577 -- points[1] + 0,577
        eta = { node->points[0], node->points[0], node->points[1], node->points[1] };

        ksiBc = { node->points[0], node->points[1], 1., 1., node->points[1] ,node->points[0], -1., -1.};
        etaBc = { -1., -1., node->points[0], node->points[1], 1., 1., node->points[1], node->points[0]};
    }
    else if (punktyCalkowania == 9) {
        ksi = { node->points[0], node->points[1], node->points[2], node->points[2], node->points[1], node->points[0], node->points[0], node->points[1], node->points[2] };
        eta = { node->points[0], node->points[0], node->points[0], node->points[1], node->points[1], node->points[1], node->points[2], node->points[2], node->points[2] };

        ksiBc = {
            node->points[0], node->points[1], node->points[2],  // dolna krawêdŸ
            1., 1., 1.,                                        // prawa krawêdŸ
            node->points[2], node->points[1], node->points[0],  // górna krawêdŸ
            -1., -1., -1.                                      // lewa krawêdŸ
        };

        etaBc = {
            -1., -1., -1.,                                     // dolna krawêdŸ
            node->points[0], node->points[1], node->points[2], // prawa krawêdŸ
            1., 1., 1.,                                        // górna krawêdŸ
            node->points[2], node->points[1], node->points[0]  // lewa krawêdŸ
        };
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

        ksiBc = {
            node->points[0], node->points[1], node->points[2], node->points[3],  // dolna krawêdŸ
            1., 1., 1., 1.,                                                     // prawa krawêdŸ
            node->points[3], node->points[2], node->points[1], node->points[0],  // górna krawêdŸ
            -1., -1., -1., -1.                                                  // lewa krawêdŸ
        };

        etaBc = {
            -1., -1., -1., -1.,                                                 // dolna krawêdŸ
            node->points[0], node->points[1], node->points[2], node->points[3], // prawa krawêdŸ
            1., 1., 1., 1.,                                                     // górna krawêdŸ
            node->points[3], node->points[2], node->points[1], node->points[0]  // lewa krawêdŸ
        };
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
    calculateMatrixHbc(punktyCalkowania, element, node, grid, etaBc, ksiBc, globaldata->Alfa);
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

void calculateMatrixHbc(int punktyCalkowania, Element* element, Node* node, Grid* grid, vector<double> ksi, vector<double> eta, double alfa)
{
    MatrixHbc matrixHbc;
    matrixHbc.calculateN(punktyCalkowania, element, grid,eta, ksi, alfa, node);
}
