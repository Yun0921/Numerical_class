#include <iostream>
#include <vector>
#include <cmath>
using namespace std;

const double PI = acos(-1);
const int nx = 10; // x: 0 ~ π, step = 0.1π, total 11 points
const int ny = 5;  // y: 0 ~ π/2, step = 0.1π, total 6 points
const double h = PI / nx;      // h = k = 0.1π
const double k = PI / (2 * ny);
const double alpha = (h * h) / (k * k);
const double tol = 1e-6;
const int MAX_ITER = 10000;

double f(double x, double y) {
    return x * y;
}

int main() {
    vector<vector<double>> u(nx + 1, vector<double>(ny + 1, 0.0));

    // 設定邊界條件
    for (int j = 0; j <= ny; ++j) {
        double y = j * k;
        u[0][j] = cos(y);         // u(0, y) = cos(y)
        u[nx][j] = -cos(y);       // u(π, y) = -cos(y)
    }

    for (int i = 0; i <= nx; ++i) {
        double x = i * h;
        u[i][0] = cos(x);         // u(x, 0) = cos(x)
        u[i][ny] = 0.0;           // u(x, π/2) = 0
    }

    for (int iter = 0; iter < MAX_ITER; ++iter) {
        double max_diff = 0.0;
        for (int i = 1; i < nx; ++i) {
            for (int j = 1; j < ny; ++j) {
                double x = i * h; // x = x0 + ih
                double y = j * k; // y = y0 + jk
                double rhs = h * h * f(x, y);
                double u_new = (u[i+1][j] + u[i-1][j] + alpha * (u[i][j+1] + u[i][j-1]) - rhs) / (2 * (1 + alpha));
                max_diff = max(max_diff, fabs(u_new - u[i][j]));
                u[i][j] = u_new;
            }
        }
        if (max_diff < tol) break;
    }

    // 輸出數值解
    
    // 先印出 i 的索引
    printf("  y\\x  ");
    for (int i = 0; i <= nx; ++i) {
        printf("%-8d ", i);
    }
    cout << "\n";
    
    for (int j = ny; j >= 0; --j) {
        printf("%5d ", j);  // 印出 j 的索引
        for (int i = 0; i <= nx; ++i) {
            printf("%8.4f ", u[i][j]);
        }
        cout << "\n";
    }

    return 0;
}
