#include <iostream>
#include <cmath>
#include <vector>
#include <iomanip>

// u1' and u2'
void F(double t, const std::vector<double>& u, std::vector<double>& dudt) {
    dudt[0] = 9.0 * u[0] + 24.0 * u[1] + 5.0 * cos(t) -  sin(t) / 3.0;
    dudt[1] = -24.0 * u[0] - 52.0 * u[1] - 9.0 * cos(t) + sin(t) / 3.0;
}

// exact value
std::vector<double> exactSolution(double t) {
    std::vector<double> u(2);
    u[0] = 2.0 * exp(-3.0 * t) - exp(-39.0 * t) + cos(t) / 3.0;
    u[1] = -exp(-3.0 * t) + 2.0 * exp(-39.0 * t) - cos(t) / 3.0;
    return u;
}

// RK4 method
void rungeKutta4(double t0, double tf, const std::vector<double>& u0, double h, 
                 std::vector<double>& t_values, 
                 std::vector<std::vector<double>>& u_values) {
    
    std::vector<double> u = u0;
    double t = t0;
    
    int n = static_cast<int>((tf - t0) / h) + 1;
    
    t_values.reserve(n);
    u_values.reserve(n);
    
    std::vector<double> k1(u.size()), k2(u.size()), k3(u.size()), k4(u.size());
    std::vector<double> temp(u.size());
    
    t_values.push_back(t);
    u_values.push_back(u);

    while (t + h <= tf + 1e-10) {
        // k1 = F(t, u)
        F(t, u, k1);
        k1[0] *= h;
        k1[1] *= h;
        // printf("k1: %.6f, %.6f\n", k1[0], k1[1]);
        
        // k2 = F(t + h/2, u + h*k1/2)
        for (size_t j = 0; j < u.size(); j++) {
            temp[j] = u[j] + 0.5 * k1[j];
        }
        F(t + 0.5 * h, temp, k2);
        k2[0] *= h;
        k2[1] *= h;
        // printf("k2: %.6f, %.6f\n", k2[0], k2[1]);
        
        // k3 = F(t + h/2, u + h*k2/2)
        for (size_t j = 0; j < u.size(); j++) {
            temp[j] = u[j] + 0.5 * k2[j];
        }
        F(t + 0.5 * h, temp, k3);
        k3[0] *= h;
        k3[1] *= h;
        // printf("k3: %.6f, %.6f\n", k3[0], k3[1]);
        
        // k4 = F(t + h, u + h*k3)
        for (size_t j = 0; j < u.size(); j++) {
            temp[j] = u[j] + k3[j];
        }
        F(t + h, temp, k4);
        k4[0] *= h;
        k4[1] *= h;
        // printf("k4: %.6f, %.6f\n", k4[0], k4[1]);
        
        for (size_t j = 0; j < u.size(); j++) {
            u[j] = u[j] + (1 / 6.0) * (k1[j] + 2.0 * k2[j] + 2.0 * k3[j] + k4[j]);
        }
        
        t += h;
        
        t_values.push_back(t);
        u_values.push_back(u);
    }
}

int main() {
    double t0 = 0.0;
    double tf = 1.0;
    std::vector<double> u0 = {4.0/3.0, 2.0/3.0};
    
    std::vector<double> step_sizes = {0.1, 0.05};

    
    for (double h : step_sizes) {
        std::vector<double> t_values;
        std::vector<std::vector<double>> u_values;
        
        rungeKutta4(t0, tf, u0, h, t_values, u_values);
        
        std::cout << "step h = " << h << " result:\n";
        std::cout << "--------------------------------------------\n";
        std::cout << std::setw(8) << "t" << std::setw(16) << "u1 (RK4)" << std::setw(16) << "u2 (RK4)";
        std::cout << std::setw(16) << "u1 " << std::setw(16) << "u2 ";
        std::cout << std::setw(16) << "error u1" << std::setw(16) << "error u2" << std::endl;
        
        for (size_t i = 0; i < t_values.size(); i++) {
            double t = t_values[i];
            std::vector<double> exact = exactSolution(t);
            
            std::cout << std::fixed << std::setprecision(6);
            std::cout << std::setw(8) << t;
            std::cout << std::setw(16) << u_values[i][0] << std::setw(16) << u_values[i][1];
            std::cout << std::setw(16) << exact[0] << std::setw(16) << exact[1];
            std::cout << std::setw(16) << exact[0] - u_values[i][0];
            std::cout << std::setw(16) << exact[1] - u_values[i][1] << std::endl;
        }
        std::cout << std::endl;
    }
    
    return 0;
}