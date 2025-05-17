#include <iostream>
#include <vector>
#include <cmath>
#include <functional>
using namespace std;

const double PI = acos(-1);

// f(x)
double f(double x) {
    return x * x * sin(x);
}

// mapping x to z
double x_to_z(double x, double c, double d) {
    return PI * (2 * (x - c) / (d - c) - 1);
}

// s_n(x) = a0/2 + Σ (a_k * cos(kz) + b_k * sin(kz)) + a_n * cos(nz)
double S_m(double z, const vector<double>& a, const vector<double>& b, int m) {
    double sum = a[0] / 2.0;
    for (int k = 1; k < m; ++k) {
        sum += a[k] * cos(k * z) + b[k] * sin(k * z);
    }
    sum += a[m] * cos(m * z);
    return sum;
}

// trapezoidal integration
double trapezoidal_integration(const function<double(double)>& func, double a, double b, int steps) {
    double h = (b - a) / steps;
    double s = 0.5 * (func(a) + func(b));
    for (int i = 1; i < steps; ++i) {
        s += func(a + i * h);
    }
    return s * h;
}

int main() {
    int m = 16;
    int N = 2 * m;
    int n = 4;
    double c = 0.0, d = 1.0;

    double dx = (d - c) / (N - 1);

    // split the interval [c, d] into N points
    vector<double> x(N), y(N), z(N);
    for (int i = 0; i < N; ++i) {
        x[i] = c + i * dx;
        y[i] = f(x[i]);
        z[i] = x_to_z(x[i], c, d);
    }

    // compute a0, a_k, b_k
    double a0 = 0;
    for (int i = 0; i < N; ++i) a0 += y[i];
    a0 /= m;

    vector<double> a(m, 0), b(m, 0);
    a[0] = a0;
    for (int k = 1; k <= n; ++k) {
        double sum_a = 0, sum_b = 0;
        for (int i = 0; i < N; ++i) {
            sum_a += y[i] * cos(k * z[i]);
            sum_b += y[i] * sin(k * z[i]);
        }
        a[k] = sum_a / m;
        b[k] = sum_b / m;
    }

    // b. compute ∫_0^1 S_4(x) dx
    double integral_S4 = trapezoidal_integration([&](double x_val) {
        double z_val = x_to_z(x_val, c, d);
        return S_m(z_val, a, b, n);
    }, c, d, 10000);

    // c. compute ∫_0^1 f(x) dx
    double integral_f = trapezoidal_integration(f, c, d, 10000);

    // d. compute the error E(S_4)
    double error = trapezoidal_integration([&](double x_val) {
        double diff = f(x_val) - S_m(x_to_z(x_val, c, d), a, b, n);
        return diff * diff;
    }, c, d, 10000);

    //cout trigonometric polynomial S_m(x) = a0/2 + Σ (a_k * cos(kz) + b_k * sin(kz)) + a_n * cos(nz)
    cout << "Discrete least squares trigonometric polynomial\n S_4(x) = ";
    cout << a[0] / 2 << " + ";
    for (int k = 1; k < n; ++k) {
        cout << a[k] << " * cos(" << k << "z) + " << b[k] << " * sin(" << k << "z) + ";
    }
    cout << a[4] << " * cos(4z)\n";

    cout << "\n(b) Integral of S_4(x) over [0,1]: " << integral_S4 << "\n";
    cout << "(c) Integral of f(x) = x^2 sin x over [0,1]: " << integral_f << "\n";
    cout << "(d) Error E(S_4): " << error << "\n";

    return 0;
}
