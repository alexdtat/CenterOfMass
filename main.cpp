#include <vector>
#include "MolecularDynamics/MolecularDynamics.h"

int main() {
    const double sigma = 3.54;
    const double epsilon = 0.00801;


    int nSteps = 10;
    auto simulation = MolecularDynamics(10, 6.6335209e-26, 1000.0, sigma, epsilon, 293.0, 10.0e-12);
    simulation.startSimulation(nSteps);

    auto xs = std::vector<double>(nSteps + 1, 0.0);
    auto ys = std::vector<double>(nSteps + 1, 0.0);
    auto zs = std::vector<double>(nSteps + 1, 0.0);

    for (int i = 0; i <= nSteps; i++) {
        xs[i] = simulation.statesHistory[i][0].coordinates[0];
        ys[i] = simulation.statesHistory[i][0].coordinates[1];
        zs[i] = simulation.statesHistory[i][0].coordinates[2];
    }
    auto particleTrajectory = std::tuple<std::vector<double>, std::vector<double>, std::vector<double>> {xs, ys, zs};

    simulation.showTrajectory(0);

    return 0;
}
