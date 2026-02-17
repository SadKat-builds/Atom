#include <cmath>
#include <iomanip>
#include <iostream>
#include <random>
#include <string>
#include <utility>
#include <vector>

struct Electron {
    int shell;
    double angle;
    double angularVelocity;
};

class AtomSimulation {
public:
    AtomSimulation(std::string elementName,
                   int protons,
                   int neutrons,
                   const std::vector<int>& shellOccupancy,
                   unsigned int seed = 42)
        : elementName_(std::move(elementName)),
          protons_(protons),
          neutrons_(neutrons),
          shellOccupancy_(shellOccupancy),
          rng_(seed) {
        initializeElectrons();
    }

    void step(double dt) {
        for (auto& electron : electrons_) {
            electron.angle = std::fmod(electron.angle + electron.angularVelocity * dt, 2.0 * kPi);
            if (electron.angle < 0.0) {
                electron.angle += 2.0 * kPi;
            }
        }
        time_ += dt;
    }

    void printState() const {
        std::cout << "Time: " << std::fixed << std::setprecision(2) << time_ << " s\n";
        std::cout << "Element: " << elementName_ << " (p=" << protons_ << ", n=" << neutrons_
                  << ", e=" << electrons_.size() << ")\n";
        std::cout << "Electron positions (x, y):\n";

        for (size_t i = 0; i < electrons_.size(); ++i) {
            const Electron& e = electrons_[i];
            const double radius = shellRadius(e.shell);
            const double x = radius * std::cos(e.angle);
            const double y = radius * std::sin(e.angle);

            std::cout << "  e" << (i + 1) << " [shell " << e.shell << "]"
                      << " -> (" << std::setw(7) << std::setprecision(3) << x << ", " << std::setw(7)
                      << std::setprecision(3) << y << ")\n";
        }
        std::cout << "----------------------------------------\n";
    }

private:
    static constexpr double kPi = 3.14159265358979323846;

    void initializeElectrons() {
        std::uniform_real_distribution<double> angleDist(0.0, 2.0 * kPi);

        for (size_t shellIndex = 0; shellIndex < shellOccupancy_.size(); ++shellIndex) {
            const int shell = static_cast<int>(shellIndex) + 1;
            const int count = shellOccupancy_[shellIndex];

            // Inner shells orbit faster in this simple toy model.
            const double baseVelocity = 3.0 / static_cast<double>(shell);

            for (int i = 0; i < count; ++i) {
                const double jitter = (static_cast<double>(i) / std::max(1, count)) * 0.2;
                electrons_.push_back(Electron{
                    shell,
                    angleDist(rng_),
                    baseVelocity + jitter,
                });
            }
        }
    }

    static double shellRadius(int shell) {
        return 0.8 * static_cast<double>(shell);
    }

    std::string elementName_;
    int protons_;
    int neutrons_;
    std::vector<int> shellOccupancy_;
    std::vector<Electron> electrons_;
    double time_ = 0.0;
    std::mt19937 rng_;
};

int main() {
    // Example: Carbon atom with Bohr-like shell occupancy [2, 4].
    AtomSimulation carbon("Carbon", 6, 6, {2, 4});

    constexpr double dt = 0.15;
    constexpr int steps = 12;

    for (int i = 0; i < steps; ++i) {
        carbon.printState();
        carbon.step(dt);
    }

    return 0;
}
