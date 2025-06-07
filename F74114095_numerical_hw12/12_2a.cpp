#include <iostream>
#include <vector>
#include <iomanip>
using namespace std;

const int N = 6;
const double dr = 0.1;
const double dt = 0.5;
const double K = 0.1;
const double alpha2 = 4 * K;
const double lambda = alpha2 * dt / (2 * dr * dr); // = 10
const int time_steps = 21;

// 初始條件 T(r,0) = 200(r - 0.5)
vector<double> initial_T() {
    vector<double> T(N);
    for (int i = 0; i < N; ++i) {
        double r = 0.5 + i * dr;
        T[i] = 200 * (r - 0.5);
    }
    return T;
}

// T(1,t) = 100 + 40t
double right_bc(double t) {
    return 100 + 40 * t;
}

// g_i = 4K/r_i * (T_{i+1} - T_{i-1}) / 2dr
vector<double> compute_g(const vector<double>& T) {
    vector<double> g(N, 0.0);
    for (int i = 1; i < N - 1; ++i) {
        double r = 0.5 + i * dr;
        g[i] = alpha2 * (1.0 / r) * (T[i + 1] - T[i - 1]) / (2 * dr);
    }
    return g;
}

// 建立 A 矩陣 (N x N) 明確存起來
vector<vector<double>> build_A_matrix() {
    vector<vector<double>> A(N, vector<double>(N, 0.0));
    for (int i = 1; i < N - 1; ++i) {
        A[i][i] = 1 - 2 * lambda;
        A[i][i - 1] = lambda;
        A[i][i + 1] = lambda;
    }
    // 處理 Neumann 邊界 (左側): T0 = T1 + 6h*T0 ⇒ 模擬虛擬點邏輯
    A[0][0] = 1;  // 用補點處理
    // Dirichlet 邊界 (右側): T_N-1 直接套值，不經由 A
    A[N - 1][N - 1] = 0;
    return A;
}

// 矩陣與向量相乘 A * T
vector<double> mat_vec_mul(const vector<vector<double>>& A, const vector<double>& T_prev) {
    vector<double> result(N, 0.0);
    for (int i = 0; i < N; ++i)
        for (int j = 0; j < N; ++j)
            result[i] += A[i][j] * T_prev[j];
    return result;
}

int main() {
    vector<vector<double>> T(time_steps, vector<double>(N));
    T[0] = initial_T();

    vector<vector<double>> A = build_A_matrix();

    for (int j = 1; j < time_steps; ++j) {
        double t_prev = (j - 1) * dt;
        double t_now = j * dt;

        const vector<double>& T_prev = T[j - 1];
        vector<double>& T_now = T[j];

        vector<double> g = compute_g(T_prev);
        vector<double> AT = mat_vec_mul(A, T_prev);

        // 虛擬點補值：T_{-1} = T1 + 6h*T0
        double T_virtual = T_prev[1] + 6 * dr * T_prev[0];
        AT[0] = T_virtual;  // 左邊界 Neumann → T_0 = T_{-1}

        // 加上 dt * g 項
        for (int i = 0; i < N; ++i)
            AT[i] += dt * g[i];

        // 加上右邊界 Dirichlet 項
        AT[N - 1] = right_bc(t_now);

        T_now = AT;  // 儲存更新後的 T
    }

    // 輸出 r 標頭
    cout << setw(6) << "t";
    for (int i = 0; i < N; ++i) {
        double r = 0.5 + i * dr;
        cout << setw(9) << "r=" << fixed << setprecision(2) << r;
    }
    cout << endl;

    // 輸出 T(r,t)
    cout << fixed << setprecision(2);
    for (int j = 0; j < time_steps; ++j) {
        double t = j * dt;
        cout << setw(6) << t;
        for (int i = 0; i < N; ++i)
            cout << setw(9) << T[j][i];
        cout << endl;
    }

    return 0;
}
