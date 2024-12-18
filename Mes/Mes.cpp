#include <iostream>
#include "Element.h"
#include "GlobalData.h"
#include "LoadGlobalData.h"
#include "Grid.h"
#include "Node.h"
#include "Calculate.h"
using namespace std;

void static run(int liczbaWezlow) {
    GlobalData* globaldata = new GlobalData;
    Grid* grid = new Grid;
    Element* element = new Element;
    loadGlobalData(globaldata, grid);
    Node* node = new Node(liczbaWezlow, globaldata);
    calculate(liczbaWezlow*liczbaWezlow, element, grid, globaldata, node);

    delete element, grid;

    finalCalculation(node, globaldata);
}

int main()
{
    run(4);
}