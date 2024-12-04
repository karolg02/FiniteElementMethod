#pragma once
#include <iomanip>
#include <iostream>
#include "Element.h"
#include "GlobalData.h"
#include <vector>
struct Node {
    double x = 0.0;
    double y = 0.0;
    vector<double> weights;
    vector<double> points;
    vector<int> weightsPathX;
    vector<int> weightsPathY;
    vector<vector < double >> H_GLOBAL;
    GlobalData* globalDataRef = nullptr;

    Node() {}

    Node(int iloscWezlow, GlobalData* globaldata) : globalDataRef(globaldata) {
        H_GLOBAL.resize(globaldata->nN, vector<double>(globaldata->nN, 0));
        if (iloscWezlow == 2) {
            addTwoNode();
        }
        else if (iloscWezlow == 3) {
            addThreeNode();
        }
        else if (iloscWezlow == 4) {
            addFourNode();
        }
    }

    void addTwoNode();

    void addThreeNode();

    void addFourNode();

    void printGlobalH();
};