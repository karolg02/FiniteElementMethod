#pragma once
#include "Element.h"
#include "Grid.h"
using namespace std;

struct MatrixHbc {
    double Hb_czastkowe[4][4] = { 0.0 };
    double Hb_calkowite[4][4] = { 0.0 };
    double detJ;

    void calculateN(int punktyCalkowania, Element* element, Grid* grid, double alfa, Node* node);
};