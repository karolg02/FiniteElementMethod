#include "Jakobian.h"
#include "Element.h"

void Jakobian::calculateMatrix(int i, Element* element)
{
    scoreX = scoreXY = scoreY = scoreYX = 0.0;
    for (int j = 0; j < 4; j++) {
        scoreX += element->Ksi[i][j] * element->pointsX[j];
        scoreXY += element->Ksi[i][j] * element->pointsY[j];
        scoreY += element->Eta[i][j] * element->pointsY[j];
        scoreYX += element->Eta[i][j] * element->pointsX[j];
    }

    macierzJakobiego[0] = scoreX;
    macierzJakobiego[1] = scoreXY;
    macierzJakobiego[2] = scoreYX;
    macierzJakobiego[3] = scoreY;
}

void Jakobian::calculateDetJ()
{
    detJ = scoreX * scoreY - scoreXY * scoreYX;
}

void Jakobian::calculateMatrixReversed()
{
    macierzOdwrotna[0] = (1 / detJ) * scoreY;
    macierzOdwrotna[1] = (-1 / detJ) * scoreXY;
    macierzOdwrotna[2] = (-1 / detJ) * scoreYX;
    macierzOdwrotna[3] = (1 / detJ) * scoreX;
}
