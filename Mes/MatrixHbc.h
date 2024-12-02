#pragma once
#include "Element.h"
struct MatrixHbc {
	double Hb_czastkowe[4][4];
	double Hb_calkowite[4][4];
	double Hbleft[4][4] = { 0 };
	double Hbright[4][4] = { 0 };
	double N1;
	double N2;
	double N3;
	double N4;

	void calculateN(Element* element) {
		N1 = 0.25 * (1 - ksi[i]) * (1 - eta[i]);
		N2 = 0.25 * (1 + ksi[i]) * (1 - eta[i]);
		N3 = 0.25 * (1 + ksi[i]) * (1 + eta[i]);
		N4 = 0.25 * (1 - ksi[i]) * (1 + eta[i]);
	}
};