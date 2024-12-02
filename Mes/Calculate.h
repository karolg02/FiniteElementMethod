#pragma once
#include "Grid.h"
#include "GlobalData.h"
#include "Jakobian.h"
#include "Element.h"

void calculate(int , Grid* grid, GlobalData* globaldata, Node* node);
void startCalculating(int, Element* element, Node* node, GlobalData* globaldata, Grid* grid);
void calculateMatrixJakobiego(int, Element* element, Jakobian* jakobian);
void calculateMatrixH(int, Element* element, Jakobian* jakobian, Node* node, double);
void calculateMatrixHbc(int, Element* element, Node* node, Grid* grid, vector<double>, vector<double>, double);
