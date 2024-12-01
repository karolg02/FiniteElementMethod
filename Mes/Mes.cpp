#include <iostream>
#include "Element.h"
#include "GlobalData.h"
#include "LoadGlobalData.h"
#include "Grid.h"
#include "Node.h"
#include "Calculate.h"

void run(int liczbaWezlow) {
    GlobalData* globaldata = new GlobalData;
    Grid* grid = new Grid;
    loadGlobalData(globaldata, grid);
    grid->printNodes();
    grid->printElements();
    grid->printBC();
    Node* node = new Node(liczbaWezlow, globaldata);
    calculate(liczbaWezlow*liczbaWezlow, grid, globaldata, node);
    node->printGlobalH();
}

int main()
{
    run(2);
}