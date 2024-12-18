#pragma once
#include "Element.h"
#include "Node.h"
#include "Grid.h"
#include "Jakobian.h"
struct MatrixC {
	double C[4][4] = { 0 };
	void calculateC(int punktyCalkowania, Element* element, Node* node, Grid* grid, GlobalData* globaldata, Jakobian* jakobian);
};