﻿#include <iostream>
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
    //GaussElimination(node->H_GLOBAL, node->P_GLOBAL);
    finalCalculation(node, globaldata);
}

int main()
{
    cout << "Ilosc punktow calkowania: " << 16 << endl;
    run(4);         //2, 3 lub 4 punktowy schemat calkowania
}