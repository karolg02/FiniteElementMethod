#pragma once
#include <vector>
#include "Element.h"
#include "Node.h"
using namespace std;
struct Grid {
    int nN = 0;
    int nE = 0;

    std::vector<Node> node;
    std::vector<Element> element;
    std::vector<int> boundaryConditions;

    void printNodes();
    void printElements();
    void printBC();
};