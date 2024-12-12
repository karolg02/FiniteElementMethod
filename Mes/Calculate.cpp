﻿#include "Calculate.h"
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
    matrixH.calculateH(punktyCalkowania, dNdx, dNdy, jakobian, node, k, element);
}

void calculateMatrixHbc(int punktyCalkowania, Element* element, Node* node, Grid* grid, GlobalData* globaldata)
{
    MatrixHbc matrixHbc;
    matrixHbc.calculateHbc(punktyCalkowania, element, grid, globaldata, node);
}

void GaussElimination(const vector<vector<double>>& H_GLOBAL, const vector<double>& P_GLOBAL) {
    int size = P_GLOBAL.size();
    vector<vector<double>> tablica(size, vector<double>(size + 1, 0.0));

    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            tablica[i][j] = H_GLOBAL[i][j];
        }
        tablica[i][size] = P_GLOBAL[i];
    }

    double wspolczynnik;
    for (int z = 0; z < size - 1; z++) {
        for (int i = 0; i < size - z - 1; i++) {
            wspolczynnik = tablica[z + i + 1][z] / tablica[z][z];
            for (int j = z; j < size + 1; j++) {
                tablica[z + i + 1][j] -= wspolczynnik * tablica[z][j];
            }
        }
    }

    vector<double> wyniki(size);
    wyniki[size - 1] = tablica[size - 1][size] / tablica[size - 1][size - 1];
    int licznik = 1;
    double suma = 0;
    for (int i = size - 2; i >= 0; i--) {
        suma = 0;
        for (int j = size; j > size - licznik; j--) {
            suma += wyniki[j - 1] * tablica[i][j - 1];
        }
        licznik++;
        wyniki[i] = (tablica[i][size] - suma) / tablica[i][i];
    }
    cout << "\nWyniki temperatury otoczenia\n\n";
    for (double val : wyniki) {
        cout << val << " ";
    }
    cout << endl;
}