#include <iostream>
#include <vector>
#include <iomanip>
using namespace std;

const int N = 6;           // 範圍從 0.5 到 1.0, Δr = 0.1 共 6 點
const double dr = 0.1;
const double dt = 0.5;
const double K = 0.1;
const double alpha2 = 4 * K;
const double lambda = alpha2 * dt / (dr * dr);
const int time_steps = 21; // 從 t = 0 到 t = 10，共 21 個時間點

// 初始條件 f(r) = 200(r - 0.5)
vector<double> initial_T() {
    vector<double> T(N);
    for (int i = 0; i < N; ++i) {
        double r = 0.5 + i * dr;
        T[i] = 200 * (r - 0.5);
    }
    return T;
}

// 邊界條件 T(1,t) = 100 + 40t
double right_bc(double t) {
    return 100 + 40 * t;
}

// g_i = 4K * (1/r_i) * (T_{i+1} - T_{i-1}) / (2*dr)
vector<double> compute_g(const vector<double>& T) {
    vector<double> g(N, 0.0);
    for (int i = 1; i < N - 1; ++i) {
        double r = 0.5 + i * dr;
        g[i] = alpha2 * (1.0 / r) * (T[i + 1] - T[i - 1]) / (2 * dr);
    }
    return g;
}

// Thomas 解三對角系統
vector<double> solve_tridiagonal(const vector<double>& a, const vector<double>& b,
                                 const vector<double>& c, const vector<double>& d) {
    int n = b.size();
    vector<double> cp(n), dp(n), x(n);
    cp[0] = c[0] / b[0];
    dp[0] = d[0] / b[0];
    for (int i = 1; i < n; ++i) {
        double denom = b[i] - a[i] * cp[i - 1];
        cp[i] = c[i] / denom;
        dp[i] = (d[i] - a[i] * dp[i - 1]) / denom;
    }
    x[n - 1] = dp[n - 1];
    for (int i = n - 2; i >= 0; --i) {
        x[i] = dp[i] - cp[i] * x[i + 1];
    }
    return x;
}

int main() {
    vector<vector<double>> T(time_steps, vector<double>(N)); // T[j][i]
    T[0] = initial_T(); // 初始條件

    // 預先建立 A 矩陣的三對角項
    vector<double> a(N - 2, -lambda);
    vector<double> b(N - 2, 1 + 2 * lambda);
    vector<double> c(N - 2, -lambda);

    for (int j = 1; j < time_steps; ++j) {
        double t_prev = (j - 1) * dt;
        double t_now = j * dt;

        vector<double>& T_prev = T[j - 1];
        vector<double>& T_now = T[j];

        // 計算 g 項
        vector<double> g = compute_g(T_prev);

        // 建立 RHS
        vector<double> d(N - 2);
        for (int i = 1; i < N - 1; ++i) {
            d[i - 1] = T_prev[i] + dt * g[i];
        }

        // Neumann 邊界修正（左邊）
        d[0] += lambda * (T_prev[1] + 6 * dr * T_prev[0]);

        // Dirichlet 邊界修正（右邊）
        d[N - 3] += lambda * right_bc(t_now);

        // 解線性系統
        vector<double> T_inner = solve_tridiagonal(a, b, c, d);

        // 邊界點填入
        T_now[0] = T_prev[1] + 6 * dr * T_prev[0]; // Neumann 邊界
        for (int i = 1; i < N - 1; ++i)
            T_now[i] = T_inner[i - 1];
        T_now[N - 1] = right_bc(t_now); // Dirichlet 邊界
    }

    // 輸出 r 標頭
    cout << setw(6) << "t";
    for (int i = 0; i < N; ++i) {
        double r = 0.5 + i * dr;
        cout << setw(8) << "r=" << fixed << setprecision(2) << r;
    }
    cout << endl;

    // 輸出 T(r,t)
    cout << fixed << setprecision(2);
    for (int j = 0; j < time_steps; ++j) {
        double t = j * dt;
        cout << setw(6) << t;
        for (int i = 0; i < N; ++i) {
            cout << setw(12) << T[j][i];
        }
        cout << endl;
    }
    return 0;
}
