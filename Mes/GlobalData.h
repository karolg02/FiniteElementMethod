#pragma once
#include <iostream>
using namespace std;
struct GlobalData {
    double SimulationTime = 500.0;
    double SimulationStepTime = 50.0;
    double Conductivity = 25.0;
    double Alfa = 300;
    double Tot = 1200;
    double InitialTemp = 100;
    double Density = 7800;
    double SpecificHeat = 700;

    int nH = 4;
    int nW = 4;

    int nN = nW * nH;
    int nE = (nH - 1) * (nW - 1);

    void getGlobalData();
};