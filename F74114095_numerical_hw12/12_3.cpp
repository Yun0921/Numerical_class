#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
using namespace std;

const double PI = acos(-1);
const int nr = 5;     // r: 0.5 ~ 1, step = 0.1 → 6 points
const int nt = 6;     // theta: 0 ~ pi/3, step = pi/18 → 7 points
const double dr = 0.1;
const double dtheta = PI / 18.0;
const double alpha = (dr / dtheta) * (dr / dtheta);

// 高斯消去法解線性方程組
vector<double> gaussianElimination(vector<vector<double>> A, vector<double> b) {
    int n = A.size();
    for (int i = 0; i < n; i++) {
        int maxRow = i;
        for (int k = i + 1; k < n; k++) {
            if (abs(A[k][i]) > abs(A[maxRow][i])) {
                maxRow = k;
            }
        }
        swap(A[i], A[maxRow]);
        swap(b[i], b[maxRow]);

        for (int k = i + 1; k < n; k++) {
            double factor = A[k][i] / A[i][i];
            b[k] -= factor * b[i];
            for (int j = i; j < n; j++) {
                A[k][j] -= factor * A[i][j];
            }
        }
    }
    vector<double> x(n);
    for (int i = n - 1; i >= 0; i--) {
        double sum = 0.0;
        for (int j = i + 1; j < n; j++) {
            sum += A[i][j] * x[j];
        }
        x[i] = (b[i] - sum) / A[i][i];
    }
    return x;
}

int main() {
    const int unknowns = (nr - 1) * (nt - 1);
    vector<vector<double>> A(unknowns, vector<double>(unknowns, 0.0));
    vector<double> F(unknowns, 0.0);
    vector<vector<double>> T(nr + 1, vector<double>(nt + 1, 0.0));

    // 設定邊界條件
    for (int i = 0; i <= nr; ++i) {
        T[i][0] = 0.0;      // theta = 0
        T[i][nt] = 0.0;     // theta = pi/3
    }
    for (int j = 0; j <= nt; ++j) {
        T[0][j] = 50.0;     // r = 0.5
        T[nr][j] = 100.0;   // r = 1
    }

    // 建立係數矩陣 A 和向量 F
    int row = 0;
    for (int j = 1; j < nt; ++j) {
        for (int i = 1; i < nr; ++i) {
            double ri = 0.5 + i * dr;
            double coeff_l = (ri * ri - 0.5 * dr * ri);
            double coeff_m = -2.0 * (alpha + ri * ri);
            double coeff_r = (ri * ri + 0.5 * dr * ri);

            A[row][row] = coeff_m;

            if (i > 1)
                A[row][row - 1] = coeff_l;

            if (i < nr - 1)
                A[row][row + 1] = coeff_r;;

            if (j > 1)
                A[row][row - (nr - 1)] = alpha;

            if (j < nt - 1)
                A[row][row + (nr - 1)] = alpha;

            if (i == 1)
                F[row] -= coeff_l * T[0][j];
            if (i == nr-1)
                F[row] -= coeff_r * T[nr][j];
            if (j == 1)
                F[row] -= alpha * T[i][0];
            if (j == nt-1)
                F[row] -= alpha * T[i][nt];

            row++;
        }
    }

    // 解線性方程組
    vector<double> sol = gaussianElimination(A, F);

    // 回填數值解
    row = 0;
    for (int j = 1; j < nt; ++j) {
        for (int i = 1; i < nr; ++i) {
            T[i][j] = sol[row++];
        }
    }

    cout << setw(6) << "t/r";
    for (int i = 0; i <= nr; ++i) {
        double r = 0.5 + i * dr;
        cout << setw(12) << r;
    }
    cout << endl;

    // 輸出 T(r,θ)
    cout << fixed << setprecision(4);
    for (int j = nt; j >= 0; --j) {
        double theta = j * dtheta;
        cout << setw(6) << theta;
        for (int i = 0; i <= nr; ++i) {
            cout << setw(12) << T[i][j];
        }
        cout << endl;
    }

    return 0;
}