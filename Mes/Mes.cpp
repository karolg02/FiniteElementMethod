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
    //grid->printNodes();
    //grid->printElements();
    //grid->printBC();
    calculate(liczbaWezlow*liczbaWezlow, element, grid, globaldata, node);
}

int main()
{
    run(2);
}