#pragma once
#include "Grid.h"
#include "GlobalData.h"
#include "Jakobian.h"
#include "Element.h"

void calculate(int ,Element* element, Grid* grid, GlobalData* globaldata, Node* node);
void calculateKsiEta(int, Element* element, Node* node, GlobalData* globaldata, Grid* grid);
void calculateMatrixH(int, Element* element, Jakobian* jakobian, Node* node, double);
void calculateMatrixHbc(int, Element* element, Node* node, Grid* grid, GlobalData* globaldata);
void GaussElimination(const vector<vector<double>>& H_GLOBAL, const vector<double>& P_GLOBAL);
