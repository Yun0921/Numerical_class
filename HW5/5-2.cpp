#include <iostream>
#include <cmath>
#include <vector>
#include <iomanip>

// 定義微分方程系統
void F(double t, const std::vector<double>& u, std::vector<double>& dudt) {
    dudt[0] = 9.0 * u[0] + 24.0 * u[1] + 5.0 * cos(t) - (1.0/3.0) * sin(t);
    dudt[1] = -24.0 * u[0] - 52.0 * u[1] - 9.0 * cos(t) + (1.0/3.0) * sin(t);
}

// 精確解
std::vector<double> exactSolution(double t) {
    std::vector<double> u(2);
    u[0] = 2.0 * exp(-3.0 * t) - exp(-39.0 * t) + (1.0/3.0) * cos(t);
    u[1] = -exp(-3.0 * t) + 2.0 * exp(-39.0 * t) - (1.0/3.0) * cos(t);
    return u;
}

// 四階 Runge-Kutta 方法
void rungeKutta4(double t0, double tf, const std::vector<double>& u0, double h, 
                 std::vector<double>& t_values, 
                 std::vector<std::vector<double>>& u_values) {
    
    std::vector<double> u = u0;  // 當前解
    double t = t0;               // 當前時間
    
    // 計算步數
    int n = static_cast<int>((tf - t0) / h) + 1;
    
    // 為向量預留空間，避免重新分配
    t_values.reserve(n);
    u_values.reserve(n);
    
    // RK4 計算所需的臨時向量
    std::vector<double> k1(u.size()), k2(u.size()), k3(u.size()), k4(u.size());
    std::vector<double> temp(u.size());
    
    // 存儲初始值
    t_values.push_back(t);
    u_values.push_back(u);
    
    // RK4 迭代
    while (t + h <= tf + 1e-10) {
        // k1 = F(t, u)
        F(t, u, k1);
        
        // k2 = F(t + h/2, u + h*k1/2)
        for (size_t j = 0; j < u.size(); j++) {
            temp[j] = u[j] + 0.5 * h * k1[j];
        }
        F(t + 0.5 * h, temp, k2);
        
        // k3 = F(t + h/2, u + h*k2/2)
        for (size_t j = 0; j < u.size(); j++) {
            temp[j] = u[j] + 0.5 * h * k2[j];
        }
        F(t + 0.5 * h, temp, k3);
        
        // k4 = F(t + h, u + h*k3)
        for (size_t j = 0; j < u.size(); j++) {
            temp[j] = u[j] + h * k3[j];
        }
        F(t + h, temp, k4);
        
        // 使用 RK4 公式更新 u
        for (size_t j = 0; j < u.size(); j++) {
            u[j] = u[j] + (h / 6.0) * (k1[j] + 2.0 * k2[j] + 2.0 * k3[j] + k4[j]);
        }
        
        t += h;
        
        // 存儲值
        t_values.push_back(t);
        u_values.push_back(u);
    }
}

int main() {
    double t0 = 0.0;  // 初始時間
    double tf = 1.0;  // 終止時間
    std::vector<double> u0 = {4.0/3.0, 2.0/3.0};  // 初始條件
    
    // 嘗試不同的步長
    std::vector<double> step_sizes = {0.1, 0.05};
    
    for (double h : step_sizes) {
        std::vector<double> t_values;
        std::vector<std::vector<double>> u_values;
        
        rungeKutta4(t0, tf, u0, h, t_values, u_values);
        
        std::cout << "step h = " << h << " result:\n";
        std::cout << "--------------------------------------------\n";
        std::cout << std::setw(8) << "t" << std::setw(15) << "u1 (RK4)" << std::setw(15) << "u2 (RK4)";
        std::cout << std::setw(15) << "u1 " << std::setw(15) << "u2 ";
        std::cout << std::setw(15) << "error u1" << std::setw(15) << "error u2" << std::endl;
        
        for (size_t i = 0; i < t_values.size(); i++) {
            double t = t_values[i];
            std::vector<double> exact = exactSolution(t);
            
            std::cout << std::fixed << std::setprecision(6);
            std::cout << std::setw(8) << t;
            std::cout << std::setw(15) << u_values[i][0] << std::setw(15) << u_values[i][1];
            std::cout << std::setw(15) << exact[0] << std::setw(15) << exact[1];
            std::cout << std::setw(15) << std::abs(u_values[i][0] - exact[0]);
            std::cout << std::setw(15) << std::abs(u_values[i][1] - exact[1]) << std::endl;
        }
        std::cout << std::endl;
    }
    
    return 0;
}
