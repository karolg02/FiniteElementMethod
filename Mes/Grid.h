#pragma once
#include <vector>
#include "Element.h"
#include "Node.h"
using namespace std;
struct Grid {
    int nN = 0;
    int nE = 0;

    vector<Node> node;
    vector<Element> element;
    vector<int> boundaryConditions;
    vector<double> ksi;
    vector<double> eta;
    vector<double> ksiBc;
    vector<double> etaBc;

    void printNodes();
    void printElements();
    void printBC();
    void settings(int punktyCalkowania, Node* node);
};