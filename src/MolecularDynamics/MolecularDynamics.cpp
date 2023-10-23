#include <valarray>
#include <random>
#include <fstream>
#include "MolecularDynamics.h"

// Variables containing current state of the generator
// Values here correspond to the default state
unsigned int m_u = 521288629, m_v = 362436069;

void SetSeed(unsigned int u, unsigned int v = unsigned(362436069)) {
    m_u = u;
    m_v = v;
}

unsigned int IntUniform(unsigned int &u, unsigned int &v) {
    v = 36969 * (v & 65535) + (v >> 16);
    u = 18000 * (u & 65535) + (u >> 16);
    return (v << 16) + u;
}

double brng() {
    unsigned int z = IntUniform(m_u, m_v);
    return z * 2.328306435996595e-10;
}

void MolecularDynamics::makeSystemDimensionless() {
    double tau = sigma * sqrt(m / epsilon);
    time /= tau;
    dt /= tau;
    temperature *=
            kB / epsilon; // dimensionless temperature = kB * T / epsilon with T the absolute temperature in Kelvins
    cubeSize /= sigma;
    halfSize /= sigma;
    for (auto particle: particles) {
        for (int j = 0; j < 3; j++) {
            particle.coordinates[j] /= sigma;
            particle.velocity[j] *= tau / sigma;
            particle.acceleration[j] *= tau * tau / sigma;
        }
    }
}

// Dimensionless
std::array<double, 3> MolecularDynamics::getForce(int particle_ind_a, int particle_ind_b) {
    Particle particleA = particles[particle_ind_a];
    Particle particleB = particles[particle_ind_b];
    double difX = particleA.coordinates[0] - particleB.coordinates[0];
    double difY = particleA.coordinates[1] - particleB.coordinates[1];
    double difZ = particleA.coordinates[2] - particleB.coordinates[2];
    double r = sqrt(difX * difX + difY * difY + difZ * difZ);
    // double sr = sigma / r;
    // double cf = -24.0 * epsilon * (pow(sr, 7.0) - 2 * pow(sr, 13.0)) / r;
    double cf = -24.0 * (pow(r, -7.0) - 2 * pow(r, -13.0)) / r;
    auto f = std::array<double, 3>();
    for (int j = 0; j < 3; j++) {
        f[j] = cf * (particleA.coordinates[j] - particleB.coordinates[j]);
    }

    return f;
}

void MolecularDynamics::updateAcceleration(int particle_ind) {
    for (int i = 0; i < particles.size(); i++) {
        if (i != particle_ind) {
            auto fi = getForce(particle_ind, i);
            for (int j = 0; j < 3; j++) {
                // f[j] += fi[j];
                particles[particle_ind].acceleration[j] += fi[j];
            }
        }
    }
}

void MolecularDynamics::updateVelocity(int particle_ind) {
    for (int j = 0; j < 3; j++) {
        particles[particle_ind].velocity[j] += particles[particle_ind].acceleration[j] * dt;
    }
}

void MolecularDynamics::updateCoordinates(int particle_ind) {
    // auto bufferParticle = bufferParticles[particle_ind];
    for (int j = 0; j < 3; j++) {
        particles[particle_ind].coordinates[j] +=
                particles[particle_ind].velocity[j] * dt + particles[particle_ind].acceleration[j] * dt * dt / 2.0;
        double diff = fabs(particles[particle_ind].coordinates[j]) - halfSize;
        if (diff > 0.0) {
            if (particles[particle_ind].coordinates[j] > 0.0) {
                particles[particle_ind].coordinates[j] = halfSize - diff;
            } else {
                particles[particle_ind].coordinates[j] = -halfSize + diff;
            }
            particles[particle_ind].velocity[j] *= -1;
        }
    }
}

std::array<double, 3> MolecularDynamics::calculateCenterOfMass() {
    std::array<double, 3> centerOfMass = {0.0, 0.0, 0.0};
    double totalMass = 0.0;

    for (Particle particle: particles) {
        totalMass += 1.0;
    }

    for (Particle particle: particles) {
        for (int j = 0; j < 3; j++) {
            centerOfMass[j] += particle.coordinates[j];
        }
    }

    for (int j = 0; j < 3; j++) {
        centerOfMass[j] /= totalMass;
    }

    return centerOfMass;
}

void MolecularDynamics::updateSystem() {
    for (int i = 0; i < particles.size(); i++) {
        updateCoordinates(i);
    }

    for (int i = 0; i < particles.size(); i++) {
        updateVelocity(i);
    }

    for (int i = 0; i < particles.size(); i++) {
        updateAcceleration(i);
    }
}

static double frand_Gaussian(double E, double V) // Here V is the variance
{
    return E + sqrt(-2.0 * V * log(brng())) * cos(2.0 * M_PI * brng());
}

void MolecularDynamics::vrand_MB(std::array<double, 3> &v, double *u) const {
    // double RT = kB * T / m; -- Variance for dimensions case
    v[0] = frand_Gaussian(u[0], temperature);
    v[1] = frand_Gaussian(u[1], temperature);
    v[2] = frand_Gaussian(u[2], temperature);
}

void MolecularDynamics::vrand_MB_eq(std::array<double, 3> &v) const {
    // double RT = kB * T / m; -- Variance for dimensions case
    v[0] = frand_Gaussian(0.0, temperature);
    v[1] = frand_Gaussian(0.0, temperature);
    v[2] = frand_Gaussian(0.0, temperature);
}


MolecularDynamics::MolecularDynamics(int nParticlesNew, double mNew, double cubeSizeNew, double sigmaNew,
                                     double epsilonNew, double temperatureNew, double dtNew) {
    cubeSize = cubeSizeNew;
    halfSize = cubeSizeNew / 2.0;

    sigma = sigmaNew;
    epsilon = epsilonNew;
    m = mNew;
    dt = dtNew;
    temperature = temperatureNew;

    particles = std::vector<Particle>(nParticlesNew, {{0.0, 0.0, 0.0},
                                                      {0.0, 0.0, 0.0},
                                                      {0.0, 0.0, 0.0}});
    makeSystemDimensionless();
    static std::random_device rd;
    static std::mt19937 gen(rd());
    static std::uniform_real_distribution dis(-halfSize, halfSize);

    for (int i = 0; i < nParticlesNew; i++) {
        for (int j = 0; j < nParticlesNew; j++) {
            particles[i].coordinates[j] = dis(gen);
        }
        vrand_MB_eq(particles[i].velocity);
    }
    for (int i = 0; i < nParticlesNew; i++) {
        updateAcceleration(i);
    }
}

void MolecularDynamics::startSimulation(int nStepsNew) {
    nSteps = nStepsNew;
    statesHistory.push_back(particles);
    for (step = 0; step < nSteps; step++) {
        updateSystem();
        statesHistory.push_back(particles);
        if (step % 5 == 0) {
            estimatedCenterOfMass = calculateCenterOfMass();
        }
        time += dt;
    }
}

void MolecularDynamics::showTrajectory(int particleInd) {
    std::ofstream particlePathFile;
    particlePathFile.open("../src/MolecularDynamics/particlePathFile.txt");

    particlePathFile << nSteps + 1 << "\n";
    for (int i = 0; i <= nSteps; i++) {
        particlePathFile << statesHistory[i][particleInd].coordinates[0] << " "
                         << statesHistory[i][particleInd].coordinates[1] << " "
                         << statesHistory[i][particleInd].coordinates[2] << "\n";
    }

    particlePathFile.close();

    std::string command = "python drawer.py -f ./particlePathFile.txt";
    std::system(command.c_str());
}
