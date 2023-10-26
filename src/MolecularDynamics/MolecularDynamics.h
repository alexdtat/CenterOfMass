//
// Created by Alexey on 22.10.2023.
//

#ifndef CENTEROFMASS_MOLECULARDYNAMICS_H
#define CENTEROFMASS_MOLECULARDYNAMICS_H


#include <array>
#include <vector>

struct Particle {
    std::array<double, 3> coordinates;
    std::array<double, 3> velocity;
    std::array<double, 3> acceleration;
};

class MolecularDynamics {
    double sigma = 3.54;
    double epsilon = 0.00801;
    double time = 0.0;
    double dt = 0.0;
    double m = 1.0;
    int step = 0;
    int nSteps = 0;
    double temperature = 0.0; // absolute temperature in Kelvins
    const double kB = 8.617333262e-5; // eV/K

    void makeSystemDimensionless();
    std::array<double, 3> getForce(int particleIndA, int particleIndB);
    void updateAcceleration(int particleInd);
    void updateVelocity(int particleInd);
    void updateCoordinates(int particleInd);
    std::array<double, 3> calculateCenterOfMass();
    void updateSystem();
    void vrand_MB(std::array<double, 3>& v, double *u) const;
    void vrand_MB_eq(std::array<double,3>& v) const;

    public:
        std::array<double, 3> centerOfMass = {0.0, 0.0, 0.0};
        std::vector<Particle> particles;
        std::vector<std::vector<Particle>> statesHistory;
        std::vector<std::array<double, 3>> centerOfMassHistory;
        double cubeSize = 1.0;
        double halfSize = 0.5;

        MolecularDynamics(int nParticlesNew, double mNew, double cubeSizeNew, double sigmaNew, double epsilonNew, double temperatureNew, double dtNew);
        void startSimulation (int nStepsNew);
        void showTrajectory (int particleInd);
        void showCenterOfMassChanges();
};


#endif //CENTEROFMASS_MOLECULARDYNAMICS_H
