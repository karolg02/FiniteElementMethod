#include "Node.h"

void Node::addTwoNode() {
    points.push_back(-0.57735);
    points.push_back(0.57735);

    weights.push_back(1.0);
    weights.push_back(1.0);
}
void Node::addThreeNode() {
    points.push_back((std::sqrt(3) / std::sqrt(5)) * (-1));
    points.push_back(0);
    points.push_back(std::sqrt(3) / std::sqrt(5));

    weights.push_back(5.0 / 9.0);
    weights.push_back(8.0 / 9.0);
    weights.push_back(5.0 / 9.0);
}
void Node::addFourNode() {
    points.push_back(-0.861136);
    points.push_back(-0.339981);
    points.push_back(0.339981);
    points.push_back(0.861136);

    weights.push_back(0.347855);
    weights.push_back(0.652145);
    weights.push_back(0.652145);
    weights.push_back(0.347855);
}
void Node::printGlobalH() {
    cout << "\nH GLOBALNE\n" << endl;
    for (int j = 0; j < H_GLOBAL.size(); ++j) {
        for (int k = 0; k < H_GLOBAL.size(); ++k) {
            cout << H_GLOBAL[j][k] << " ";
        }
        cout << endl;
    };
}
