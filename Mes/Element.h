#pragma once
struct Element {
    int ID[4];
    double Ksi[16][4];
    double Eta[16][4];
    double pointsX[4];
    double pointsY[4];
    int obecne[4];
    void addPointsX(double x1, double x2, double x3, double x4);
    void addPointsY(double y1, double y2, double y3, double y4);
};