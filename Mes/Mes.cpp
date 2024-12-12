#include <iostream>
#include "Element.h"
#include "GlobalData.h"
#include "LoadGlobalData.h"
#include "Grid.h"
#include "Node.h"
#include "Calculate.h"

void static run(int liczbaWezlow) {
    GlobalData* globaldata = new GlobalData;
    Grid* grid = new Grid;
    Element* element = new Element;
    loadGlobalData(globaldata, grid);
    Node* node = new Node(liczbaWezlow, globaldata);
    calculate(liczbaWezlow*liczbaWezlow, element, grid, globaldata, node);
    GaussElimination(node->H_GLOBAL, node->P_GLOBAL);
}

int main()
{
    run(4);
}