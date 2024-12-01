#pragma once
#include "Calculate.h"
#include "Grid.h"
#include "Element.h"
#include "Node.h"
#include "Jakobian.h"
struct MatrixH {
    double H[4][4];                 //macierz czastkowa
    double Hc[4][4] = { 0 };       //macierz ostateczna H dla danych PC
    double Hleft[4][4] = { 0 };
    double Hright[4][4] = { 0 };
    int size = 0;

    void calculateDerivetivedN(int punktyCalkowania, vector<vector<double>>& dNdx, vector<vector<double>>& dNdy, Element* element, Jakobian* jakobian);

    void calculateH(int punktyCalkowania, vector<vector<double>>& dNdx, vector<vector<double>>& dNdy, Jakobian* jakobian, Node* node, double k, Element* element);
};