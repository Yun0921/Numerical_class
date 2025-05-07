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
    vector<double> x(n, 0);
    vector<double> v(n);
    vector<double> Av(n);

    double tol = 1e-6;
    int max_iter = 10000;
    int iter = 0;

    for (int i = 0; i < n; i++) {
        v[i] = b[i];
        for (int j = 0; j < n; j++)
            v[i] -= A[i][j] * x[j];
    }

    for (int k = 0; k < max_iter; k++) {
        for (int i = 0; i < n; i++) {
            Av[i] = 0;
            for (int j = 0; j < n; j++)
                Av[i] += A[i][j] * v[j];
        }
        double vTv = 0, vTAv = 0;
        for (int i = 0; i < n; i++) {
            vTv += v[i] * v[i];
            vTAv += v[i] * Av[i];
        }
        double t = vTv / vTAv;
        for (int i = 0; i < n; i++)
            x[i] += t * v[i];
        vector<double> v_new(n);
        for (int i = 0; i < n; i++) {
            v_new[i] = b[i];
            for (int j = 0; j < n; j++)
                v_new[i] -= A[i][j] * x[j];
        }
        // 收斂判斷
        double err = 0;
        for (int i = 0; i < n; i++)
            err += v_new[i] * v_new[i];
        err = sqrt(err);
        if (err < tol) {
            iter = k + 1;
            break;
        }
        v = v_new;
    }
    
    cout << "(d) Conjugate Gradient Method Solution :\n";
    for (int i = 0; i < n; i++)
        cout << "x" << i+1 << " = " << x[i] << endl;
    cout << "iterations: " << iter << endl;
    return 0;
}
