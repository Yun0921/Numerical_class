#include <iostream>
#include <vector>
#include <cmath>
#include <functional>
#include <iomanip>

using namespace std;

const double PI = acos(-1);
const double h = 0.1;
const int N = 10;
const double Y_0 = 1.0; // y(0) = 1
const double Y_1 = 2.0; // y(1) = 2

inline double p(double x) { return -(x + 1); }
inline double q(double x) { return 2.0; }
inline double r(double x) { return (1 - x * x) * exp(-x); }

// === Shooting Method ===
double f1(double x, double y1, double y2) { return y2; }
double f2(double x, double y1, double y2) {
    return p(x) * y2 + q(x) * y1 + r(x);
}

void rungeKutta(double y1_0, double y2_0, vector<double>& y) {
    double x = 0.0, y1 = y1_0, y2 = y2_0;
    y[0] = y1;
    for (int i = 1; i <= N; ++i) {
        double k1 = h * f1(x, y1, y2);
        double l1 = h * f2(x, y1, y2);
        double k2 = h * f1(x + h / 2, y1 + k1 / 2, y2 + l1 / 2);
        double l2 = h * f2(x + h / 2, y1 + k1 / 2, y2 + l1 / 2);
        double k3 = h * f1(x + h / 2, y1 + k2 / 2, y2 + l2 / 2);
        double l3 = h * f2(x + h / 2, y1 + k2 / 2, y2 + l2 / 2);
        double k4 = h * f1(x + h, y1 + k3, y2 + l3);
        double l4 = h * f2(x + h, y1 + k3, y2 + l3);
        y1 += (k1 + 2 * k2 + 2 * k3 + k4) / 6;
        y2 += (l1 + 2 * l2 + 2 * l3 + l4) / 6;
        x += h;
        y[i] = y1;
    }
}

vector<double> shootingMethod() {
    vector<double> y1_vec(N + 1), y2_vec(N + 1), y(N + 1);

    rungeKutta(Y_0, 0.0, y1_vec);
    rungeKutta(0.0, 1.0, y2_vec);

    double c = (Y_1 - y1_vec[N]) / y2_vec[N];

    for (int i = 0; i <= N; ++i) {
        y[i] = y1_vec[i] + c * y2_vec[i];
    }
    return y;
}

// === Finite-Difference Method ===
vector<double> thomasSolve(vector<double>& a, vector<double>& b, vector<double>& c, vector<double>& d) {
    int n = b.size();
    vector<double> cp(n), dp(n), x(n);
    cp[0] = c[0] / b[0]; dp[0] = d[0] / b[0];
    for (int i = 1; i < n; ++i) {
        double denom = b[i] - a[i] * cp[i - 1];
        cp[i] = c[i] / denom;
        dp[i] = (d[i] - a[i] * dp[i - 1]) / denom;
    }
    x[n - 1] = dp[n - 1];
    for (int i = n - 2; i >= 0; --i)
        x[i] = dp[i] - cp[i] * x[i + 1];
    return x;
}

vector<double> finiteDifferenceMethod() {
    vector<double> a(N - 1), b(N - 1), c(N - 1), d(N - 1);

    for (int i = 1; i < N; ++i) {
        double x = i * h;
        double pi = p(x), qi = q(x), ri = r(x);

        a[i - 1] = -1.0 - 0.5 * h * pi; // A_{i,i-1}
        b[i - 1] = 2.0 + h * h * qi; // A_{i,i}
        c[i - 1] = -1.0 + 0.5 * h * pi; // A_{i,i+1}
        d[i - 1] = -h * h * ri;
    }

    // printf the matrix A
    // for(int i = 0; i < N - 1; ++i) {
    //     cout << "Row " << i + 1 << ": ";
    //     for (int j = 0; j < N - 1; ++j) {
    //         if (j == i - 1) cout << a[i] << " ";
    //         else if (j == i) cout << b[i] << " ";
    //         else if (j == i + 1) cout << c[i] << " ";
    //         else cout << "0.0 ";
    //     }
    //     cout << "| " << d[i] << endl;
    // }

    d[0] += (1.0 + 0.5 * h * p(h)) * Y_0; 
    d[N - 2] += (1.0 - 0.5 * h * p(1.0 - h)) * Y_1;

    vector<double> y_inner = thomasSolve(a, b, c, d);
    vector<double> y = {Y_0};
    y.insert(y.end(), y_inner.begin(), y_inner.end());
    y.push_back(Y_1);
    return y;
}

// === Variation Method ===
double y1_interp(double x) { return (1 - x) * Y_0 + x * Y_1; } // 線形插值 y(0)=1, y(1)=2
double phi(int i, double x) { return sin(i * PI * x); } // 基底
double dphi(int i, double x) { return i * PI * cos(i * PI * x); } // 基底一階導數

double integrate(function<double(double)> func, double a = 0, double b = 1) {
    // Gaussian quadrature (5-point)
    vector<double> x = {-0.9061798459, -0.5384693101, 0, 0.5384693101, 0.9061798459};
    vector<double> w = {0.2369268850, 0.4786286705, 0.5688888889, 0.4786286705, 0.2369268850};
    double sum = 0;
    for (int i = 0; i < 5; ++i) {
        double xi = 0.5 * (b - a) * x[i] + 0.5 * (b + a);
        sum += w[i] * func(xi);
    }
    return 0.5 * (b - a) * sum;
}

vector<double> variationApproach() {
    const int M = 5;
    vector<vector<double>> A(M, vector<double>(M));
    vector<double> b(M), c(M);
    for (int i = 0; i < M; ++i) {
        for (int j = 0; j < M; ++j) {
            A[i][j] = integrate([=](double x) {
                return p(x) * dphi(i + 1, x) * dphi(j + 1, x) + q(x) * phi(i + 1, x) * phi(j + 1, x);
            });
        }
        
        double y1p = (Y_1 - Y_0) / (1-0); // y1'(x) = (y(1) - y(0)) / (1 - 0) = 1
        double dp = -1; // (-(1+x))'
        b[i] = integrate([=](double x) {
            return (r(x) + dp*y1p - q(x) * y1_interp(x)) * phi(i + 1, x);
        });
    }

    // Gaussian elimination
    for (int i = 0; i < M; ++i) {
        for (int k = i + 1; k < M; ++k) {
            double t = A[k][i] / A[i][i];
            for (int j = i; j < M; ++j) A[k][j] -= t * A[i][j];
            b[k] -= t * b[i];
        }
    }
    for (int i = M - 1; i >= 0; --i) {
        c[i] = b[i];
        for (int j = i + 1; j < M; ++j)
            c[i] -= A[i][j] * c[j];
        c[i] /= A[i][i];
    }

    vector<double> y(N + 1);
    for (int i = 0; i <= N; ++i) {
        double x = i * h, y2 = 0;
        for (int j = 0; j < M; ++j) y2 += c[j] * phi(j + 1, x);
        y[i] = y1_interp(x) + y2;
    }
    return y;
}

int main() {
    vector<double> y_a = shootingMethod();
    vector<double> y_b = finiteDifferenceMethod();
    vector<double> y_c = variationApproach();

    cout << fixed << setprecision(6);
    cout << left << setw(10) << "x" 
         << setw(15) << "Shooting" 
         << setw(15) << "FiniteDiff" 
         << setw(15) << "Variation" << endl;

    for (int i = 0; i <= N; ++i) {
        double x = i * h;
        cout << left << setw(10) << x 
             << setw(15) << y_a[i] 
             << setw(15) << y_b[i] 
             << setw(15) << y_c[i] << endl;
    }
    return 0;
}