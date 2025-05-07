#include <iostream>
#include <cmath>
using namespace std;

int main() {
    const int n = 6;
    double A[n][n] = {
        {4, -1,  0, -1,  0,  0},
        {-1, 4, -1,  0, -1,  0},
        {0, -1, 4,  0,  0, -1},
        {-1, 0,  0, 4, -1,  0},
        {0, -1, 0, -1, 4, -1},
        {0,  0, -1, 0, -1, 4}
    };
    double b[n] = {0, -1, 9, 4, 8, 6};
    double x[n] = {0, 0, 0, 0, 0, 0};
    double omega = 1.25; // 鬆弛參數，可調整
    double tol = 1e-6;
    int max_iter = 10000, iter = 0;

    while (iter++ < max_iter) {
        double err = 0;
        for (int i = 0; i < n; i++) {
            double sum = b[i];
            for (int j = 0; j < n; j++)
                if (j != i) sum -= A[i][j] * x[j];
            double x_old = x[i];
            x[i] = (1 - omega) * x[i] + omega * (sum / A[i][i]);
            err += fabs(x[i] - x_old);
        }
        if (err < tol) break;
    }
    cout << "(c) SOR Method Solution:\n";
    for (int i = 0; i < n; i++)
        cout << "x" << i+1 << " = " << x[i] << endl;
    cout << "iterations: " << iter << endl;
    return 0;
}
