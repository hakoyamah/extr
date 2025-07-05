#include <iostream>
#include <random>
#include <cmath>
#include <fstream>

/* --- Simulation parameters --- */
constexpr double DT    = 0.001;                 // Time step
constexpr int N_REP    = 10000;                 // Number of trials

/* --- Model parameters --- */
constexpr double MU    = 0.0023556;             // Log growth rate
const double SIGMA    = std::sqrt(0.01087);     // Environmental standard deviation
constexpr double NQ    = 47.0;                  // Latest population size
constexpr double NE    = 10.0;                  // Extinction threshold
constexpr double TIME_LIMIT = 600.0;         // Time horizon for extinction probability (years)
const double XD = std::log(NQ / NE);            // Initial log ratio

int main() {
    std::mt19937 mt(std::random_device{}());
    std::uniform_real_distribution<double> dist(0.0, 1.0);
    
    std::ofstream ofs("extinction_times.csv");
    if (!ofs) {
        std::cerr << "Error opening extinction_times.csv for writing\n";
        return 1;
    }
    ofs << "extinction_time\n";
    
    int extinct_within_limit = 0;
    double te_sum = 0.0;
    double te2_sum = 0.0;
    
    for (int i = 0; i < N_REP; ++i) {
        double x = std::log(NQ);
        int step = 0;
        bool extinct = false;
        
        while (x > std::log(NE) && DT * step < TIME_LIMIT) {
            double noise = std::sqrt(3.0 * DT) * (2.0 * dist(mt) - 1.0);
            x += MU * DT + SIGMA * noise;
            ++step;
            
            if (x < std::log(NE)) {
                extinct = true;
                break;
            }
        }
        
        if (extinct) {
            double extinction_time = DT * step;
            te_sum += extinction_time;
            te2_sum += extinction_time * extinction_time;
            ++extinct_within_limit;
            ofs << extinction_time << "\n";
        }
    }
    
    ofs.close();
    
    double mean_te = te_sum / extinct_within_limit;
    double var_te = (te2_sum / extinct_within_limit) - (mean_te * mean_te);
    double sd_te = std::sqrt(var_te);
    
    // Analytical E[T] and Var[T] (Dennis式 18, 19)
    double theory_mean = XD / std::abs(MU);
    double theory_var = XD * SIGMA * SIGMA / std::pow(std::abs(MU), 3);
    double theory_sd = std::sqrt(theory_var);
    
    std::cout << "=== Simulation result (time limit = " << TIME_LIMIT << " years) ===\n";
std::cout << "Extinction count: " << extinct_within_limit << " / " << N_REP << "\n";
std::cout << "Mean extinction time (only extinct runs): " << mean_te << "\n";
std::cout << "Variance of extinction time (only extinct runs): " << var_te << "\n";
std::cout << "SD of extinction time (only extinct runs): " << sd_te << "\n\n";

std::cout << "=== Analytical result ===\n";
std::cout << "Theoretical E[T]: " << theory_mean << "\n";
std::cout << "Theoretical Var[T]: " << theory_var << "\n";
std::cout << "Theoretical SD[T]: " << theory_sd << "\n";

return 0;
}
