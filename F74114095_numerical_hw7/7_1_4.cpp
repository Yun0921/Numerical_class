#include <iostream>
#include <vector>
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
    vector<double> x(n, 0);      // 當前解
    vector<double> x_old(n);     // 前一次的解
    vector<double> v(n);         // 殘差/方向向量
    vector<double> Av(n);        // A*v

    double tol = 1e-6;
    int max_iter = 10000;
    int iter = 0;

    // 初始化殘差 v = b - A*x
    for (int i = 0; i < n; i++) {
        v[i] = b[i];
        for (int j = 0; j < n; j++)
            v[i] -= A[i][j] * x[j];
    }

    for (int k = 0; k < max_iter; k++) {
        x_old = x; // 保存前一次的解

        // 計算 Av = A*v
        for (int i = 0; i < n; i++) {
            Av[i] = 0;
            for (int j = 0; j < n; j++)
                Av[i] += A[i][j] * v[j];
        }

        // 計算步長 t = (v·v)/(v·Av)
        double vTv = 0, vTAv = 0;
        for (int i = 0; i < n; i++) {
            vTv += v[i] * v[i];
            vTAv += v[i] * Av[i];
        }
        double t = vTv / vTAv;

        // 更新解 x = x + t*v
        for (int i = 0; i < n; i++) {
            x_old[i] = x[i]; // 保存前一次的解
            x[i] += t * v[i];
        }

        // 計算誤差 err = ||x - x_old||
        double err = 0;
        for (int i = 0; i < n; i++) {
            err += fabs(x[i] - x_old[i]);
        }
        if (err < tol) {
            iter = k + 1; // 計算實際迭代次數
            break;
        }

        // 更新殘差 v = b - A*x
        for (int i = 0; i < n; i++) {
            v[i] = b[i];
            for (int j = 0; j < n; j++)
                v[i] -= A[i][j] * x[j];
        }
    }

    cout << "(d) Conjugate Gradient Method Solution:\n";
    for (int i = 0; i < n; i++)
        cout << "x" << i+1 << " = " << x[i] << endl;
    //cout << "iterations: " << iter << endl;
    return 0;
}
