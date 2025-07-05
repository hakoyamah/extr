#include <iostream>
#include <random>
#include <cmath>
#include <fstream>

/* --- Simulation parameters --- */
constexpr double DT    = 0.001;                 // Time step
constexpr int N_REP    = 1000000;                // Number of trials

/* --- Model parameters --- */
constexpr double MU    = 0.0023556;           // Log growth rate
const double SIGMA    = std::sqrt(0.01087);   // Environmental standard deviation
constexpr double NQ    = 47.0;                   // Latest population size
constexpr double NE    = 10.0;                   // Extinction threshold
constexpr double TIME_LIMIT = 100.0;             // Time horizon for extinction probability (years)
const double XD = std::log(NQ / NE);             // Initial log ratio

// Normal cumulative distribution function (equivalent to pnorm)
double pnorm(double x) {
    return 0.5 * std::erfc(-x / std::sqrt(2.0));
}

// w function
double w(double mu, double xd, double s, double tt) {
    return (mu * tt + xd) / std::sqrt(s * tt);
}

// z function
double z(double mu, double xd, double s, double tt) {
    return (-mu * tt + xd) / std::sqrt(s * tt);
}

// Extinction probability function pr (cumulative inverse Gaussian distribution)
double pr(double mu, double xd, double s, double tt) {
    double w_val = w(mu, xd, s, tt);
    double z_val = z(mu, xd, s, tt);
    
    if (z_val < 35.0) {
        return pnorm(-w_val) + std::exp((z_val * z_val - w_val * w_val) / 2.0) * pnorm(-z_val);
    } else {
        double z_inv = 1.0 / z_val;
        double series = z_inv
                - z_inv * z_inv * z_inv / 3.0
                + 3.0 * std::pow(z_inv, 5)
                - 15.0 * std::pow(z_inv, 7)
                + 105.0 * std::pow(z_inv, 9)
                - 945.0 * std::pow(z_inv, 11)
                + 10395.0 * std::pow(z_inv, 13);
        return pnorm(-w_val) + std::exp(-0.5 * w_val * w_val) * std::sqrt(2.0) / (2.0 * std::sqrt(M_PI)) * series;
    }
}

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
    double te_sum = 0.0, te2_sum = 0.0;
    
    for (int i = 0; i < N_REP; ++i) {
        double x = XD;
        int step = 0;
        bool extinct = false;
        
        while ((DT * step) < TIME_LIMIT) {
            if (x <= 0.0) {
                extinct = true;
                break;
            }
            
            double noise = (2.0 * dist(mt) - 1.0) * std::sqrt(3.0);
            x += MU * DT + SIGMA * std::sqrt(DT) * noise;
            ++step;
        }
        
        double extinction_time = DT * step;
        if (extinct) {
            ++extinct_within_limit;
            te_sum += extinction_time;
            te2_sum += extinction_time * extinction_time;
            ofs << extinction_time << "\n";  // Record for reading in R
        }
    }
    
    ofs.close();
    
    double extinction_prob_sim = static_cast<double>(extinct_within_limit) / N_REP;
  double mean_te = (extinct_within_limit > 0) ? te_sum / extinct_within_limit : 0.0;
    double var_te = 0.0;
    if (extinct_within_limit > 1) {
        var_te = (te2_sum / extinct_within_limit) - (mean_te * mean_te);
    }
    
    std::cout << "=== Simulation result (time limit = " << TIME_LIMIT << " years) ===\n";
std::cout << "Extinction count: " << extinct_within_limit << " / " << N_REP << "\n";
std::cout << "Probability of extinction within " << TIME_LIMIT << " years (simulation): " << extinction_prob_sim << "\n";
if (extinct_within_limit > 0) {
    std::cout << "Mean extinction time (only extinct runs): " << mean_te << "\n"
       << "Variance of extinction time (only extinct runs): " << var_te << "\n";
} else {
    std::cout << "No extinctions within time limit in simulation.\n";
}

double sigma2 = SIGMA * SIGMA;
double w_val = w(MU, XD, sigma2, TIME_LIMIT);
double z_val = z(MU, XD, sigma2, TIME_LIMIT);
double extinction_prob_theory = pr(MU, XD, sigma2, TIME_LIMIT);

std::cout << "\n=== Analytical result ===\n";
std::cout << "w = " << w_val << "\n";
std::cout << "z = " << z_val << "\n";
std::cout << "Probability of extinction within " << TIME_LIMIT << " years (analysis): " << extinction_prob_theory << "\n";

return 0;
}
