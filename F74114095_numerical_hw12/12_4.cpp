#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
using namespace std;

const double PI = acos(-1);
const double alpha = 1.0;
const double dx = 0.1;
const double dt = 0.1;
const int Nx = 11;  // 包含邊界點
const int Nt = 11;
const int N = Nx - 2;  // 內部節點數
const double lambda2 = pow(alpha * dt / dx, 2);

// 初始位置
double f(double x) {
    return cos(2 * PI * x);
}

// 初始速度
double g(double x) {
    return 2 * PI * sin(2 * PI * x);
}

// 左右邊界
double left_bc(double t) { return 1.0; }
double right_bc(double t) { return 2.0; }

// 矩陣與向量相乘
vector<double> mat_vec_mul(const vector<vector<double>>& A, const vector<double>& v) {
    int n = v.size();
    vector<double> result(n, 0.0);
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            result[i] += A[i][j] * v[j];
    return result;
}

int main() {
    // 完整解 p[t][x]
    vector<vector<double>> p(Nt, vector<double>(Nx, 0.0));

    // 初始條件 t = 0
    for (int i = 0; i < Nx; ++i) {
        double x = i * dx;
        p[0][i] = f(x);
    }
    p[0][0] = left_bc(0);
    p[0][Nx - 1] = right_bc(0);

    // 第一步：用初始速度計算 u^1
    for (int i = 1; i < Nx - 1; ++i) {
        double x = i * dx;
        p[1][i] = f(x) + dt * g(x)
                + 0.5 * lambda2 * (f(x - dx) - 2 * f(x) + f(x + dx));
    }
    p[1][0] = left_bc(dt);
    p[1][Nx - 1] = right_bc(dt);

    // 建立矩陣 A (size N x N)
    vector<vector<double>> A(N, vector<double>(N, 0.0));
    for (int i = 0; i < N; ++i) {
        A[i][i] = 2 * (1 - lambda2);
        if (i > 0) A[i][i - 1] = lambda2;
        if (i < N - 1) A[i][i + 1] = lambda2;
    }

    // 時間迴圈：從 j=1 開始計算 u^{j+1}
    for (int j = 1; j < Nt - 1; ++j) {
        // 取出內部點 (不含邊界) 的 u^j 和 u^{j-1}
        vector<double> u_now(N), u_prev(N);
        for (int i = 0; i < N; ++i) {
            u_now[i] = p[j][i + 1];
            u_prev[i] = p[j - 1][i + 1];
        }

        // 做 A * u^j - u^{j-1}
        vector<double> u_next = mat_vec_mul(A, u_now);
        for (int i = 0; i < N; ++i)
            u_next[i] -= u_prev[i];

        // 邊界修正
        double t_next = (j + 1) * dt;
        u_next[0] += lambda2 * left_bc(j * dt);
        u_next[N - 1] += lambda2 * right_bc(j * dt);

        // 寫回完整解 p[j+1][i]
        for (int i = 0; i < N; ++i)
            p[j + 1][i + 1] = u_next[i];
        p[j + 1][0] = left_bc(t_next);
        p[j + 1][Nx - 1] = right_bc(t_next);
    }

    // 輸出表頭
    cout << setw(6) << "t\\x";
    for (int i = 0; i < Nx; ++i)
        cout << setw(9) << fixed << setprecision(2) << i * dx;
    cout << endl;

    // 輸出結果表格
    for (int j = 0; j < Nt; ++j) {
        cout << setw(6) << fixed << setprecision(2) << j * dt;
        for (int i = 0; i < Nx; ++i)
            cout << setw(9) << p[j][i];
        cout << endl;
    }

    return 0;
}
