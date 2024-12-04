#include "Grid.h"
#include <iostream>
#include <iomanip>
void Grid::printNodes() {
    cout << "Nodes:" << endl;
    for (int i = 0; i < nN; i++) {
        cout << setw(3) << "[" << i + 1 << "] " << setw(10) << node[i].x << " " << setw(10) << node[i].y << endl;
    }
    cout << "\n";
}
void Grid::printElements() {
    cout << "Elements:" << endl;
    for (int i = 0; i < nE; i++) {
        cout << setw(3) << "[" << i + 1 << "] " << setw(3) << element[i].ID[0] << ", " << setw(3) << element[i].ID[1] << ", " << setw(3) << element[i].ID[2] << ", " << setw(3) << element[i].ID[3] << endl;
    }
    cout << "\n";
}
void Grid::printBC() {
    cout << "BC:" << endl;
    for (int i : boundaryConditions) {
        std::cout << i << " ";
    }
    cout << "\n";
}

void Grid::settings(int punktyCalkowania, Node* node)
{
    {
        if (punktyCalkowania == 4) {
            ksi = { node->points[0], node->points[1], node->points[1], node->points[0] }; //points[0] - 0,577 -- points[1] + 0,577
            eta = { node->points[0], node->points[0], node->points[1], node->points[1] };

            ksiBc = { node->points[0], node->points[1], 1., 1., node->points[1] ,node->points[0], -1., -1. };
            etaBc = { -1., -1., node->points[0], node->points[1], 1., 1., node->points[1], node->points[0] };
        }
        else if (punktyCalkowania == 9) {
            ksi = { node->points[0], node->points[1], node->points[2], node->points[2], node->points[1], node->points[0], node->points[0], node->points[1], node->points[2] };
            eta = { node->points[0], node->points[0], node->points[0], node->points[1], node->points[1], node->points[1], node->points[2], node->points[2], node->points[2] };

            ksiBc = {
                node->points[0], node->points[1], node->points[2],  // dolna kraw�d�
                1., 1., 1.,                                        // prawa kraw�d�
                node->points[2], node->points[1], node->points[0],  // g�rna kraw�d�
                -1., -1., -1.                                      // lewa kraw�d�
            };

            etaBc = {
                -1., -1., -1.,                                     // dolna kraw�d�
                node->points[0], node->points[1], node->points[2], // prawa kraw�d�
                1., 1., 1.,                                        // g�rna kraw�d�
                node->points[2], node->points[1], node->points[0]  // lewa kraw�d�
            };
        }
        else if (punktyCalkowania == 16) {
            ksi = { node->points[0], node->points[1], node->points[2], node->points[3],
            node->points[0], node->points[1], node->points[2], node->points[3],
            node->points[0], node->points[1], node->points[2], node->points[3],
            node->points[0], node->points[1], node->points[2], node->points[3] };

            eta = { node->points[0], node->points[0], node->points[0], node->points[0],
                    node->points[1], node->points[1], node->points[1], node->points[1],
                    node->points[2], node->points[2], node->points[2], node->points[2],
                    node->points[3], node->points[3], node->points[3], node->points[3] };

            ksiBc = {
                node->points[0], node->points[1], node->points[2], node->points[3],  // dolna kraw�d�
                1., 1., 1., 1.,                                                     // prawa kraw�d�
                node->points[3], node->points[2], node->points[1], node->points[0],  // g�rna kraw�d�
                -1., -1., -1., -1.                                                  // lewa kraw�d�
            };

            etaBc = {
                -1., -1., -1., -1.,                                                 // dolna kraw�d�
                node->points[0], node->points[1], node->points[2], node->points[3], // prawa kraw�d�
                1., 1., 1., 1.,                                                     // g�rna kraw�d�
                node->points[3], node->points[2], node->points[1], node->points[0]  // lewa kraw�d�
            };
        }
        else {
            std::cerr << "Nieprawid�owa liczba punkt�w calkowania!" << std::endl;
            return;
        }
    }
}
