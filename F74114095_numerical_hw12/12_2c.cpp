#include <iostream>
#include <vector>
#include <iomanip>
using namespace std;

const int N = 6;
const double dr = 0.1;
const double dt = 0.5;
const double K = 0.1;
const double alpha2 = 4 * K;
const double lambda = alpha2 * dt / (dr * dr);
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

// 邊界條件 T(1,t) = 100 + 40t
double right_bc(double t) {
    return 100 + 40 * t;
}

// 建立 A_L（左矩陣）
void build_AL(vector<double>& a, vector<double>& b, vector<double>& c) {
    for (int i = 0; i < N - 2; ++i) {
        a[i] = -0.5 * lambda;
        b[i] = 1 + lambda;
        c[i] = -0.5 * lambda;
    }
}

// 建立 A_R（右矩陣）
vector<double> apply_AR(const vector<double>& T_prev, double t_now) {
    vector<double> b(N - 2, 0.0);
    for (int i = 1; i < N - 1; ++i) {
        double left = (i == 1) ? T_prev[1] + 6 * dr * T_prev[0] : T_prev[i - 1];
        double right = (i == N - 2) ? right_bc(t_now) : T_prev[i + 1];
        b[i - 1] = (1 - lambda) * T_prev[i] + 0.5 * lambda * (left + right);
    }
    return b;
}

// Thomas 解法
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
    vector<vector<double>> T(time_steps, vector<double>(N));
    T[0] = initial_T();

    vector<double> a(N - 2), b_diag(N - 2), c(N - 2);
    build_AL(a, b_diag, c);

    for (int j = 1; j < time_steps; ++j) {
        double t_now = j * dt;

        const vector<double>& T_prev = T[j - 1];
        vector<double>& T_now = T[j];

        vector<double> rhs = apply_AR(T_prev, t_now);

        // 右邊邊界修正
        rhs[N - 3] += 0.5 * lambda * right_bc(t_now);

        // 左邊邊界修正
        rhs[0] += 0.5 * lambda * (T_prev[1] + 6 * dr * T_prev[0]);

        // 解三對角系統
        vector<double> T_inner = solve_tridiagonal(a, b_diag, c, rhs);

        // 填入解
        T_now[0] = T_prev[1] + 6 * dr * T_prev[0]; // Neumann
        for (int i = 1; i < N - 1; ++i)
            T_now[i] = T_inner[i - 1];
        T_now[N - 1] = right_bc(t_now); // Dirichlet
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
