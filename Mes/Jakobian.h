#pragma once
#include "Element.h"
struct Jakobian {
    double scoreX = 0.0;
    double scoreY = 0.0;
    double scoreXY = 0.0;
    double scoreYX = 0.0;
    double macierzJakobiego[4] = { 0 };
    double detJ = 0.0;
    double macierzOdwrotna[4] = { 0 };

    void calculateMatrix(int i, Element* element);

    void calculateDetJ();

    void calculateMatrixReversed();
};