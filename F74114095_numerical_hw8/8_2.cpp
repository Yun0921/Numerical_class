#include <iostream>
#include <vector>
#include <cmath>
#include <functional>
using namespace std;

// f(x)
double f(double x) {
    return 0.5 * cos(x) + 0.25 * sin(2 * x);
}


// P0 = 1
// P1 = x
// Pn = x*x - 1/3.0
double legendre(int n, double x) {
    if (n == 0) return 1.0;
    if (n == 1) return x;
    if (n == 2) return x*x - 1/3.0;
}

// trapezoidal rule for numerical integration
double trapezoidal(const function<double(double)>& func, double a, double b, int n) {
    double h = (b - a) / n;
    double sum = 0.5 * (func(a) + func(b));
    for (int i = 1; i < n; ++i) {
        sum += func(a + i*h);
    }
    return sum * h;
}

// P2(x) = c0*P0(x) + c1*P1(x) + c2*P2(x)
double p(double x, const vector<double>& c) {
    double val = 0;
    for (int k = 0; k < (int)c.size(); ++k) {
        val += c[k] * legendre(k, x);
    }
    return val;
}

// compute the sum of squared errors (SSE)
double compute_sse(const vector<double>& c, int n_points = 10000) {
    double sse = 0;
    double h = 2.0 / n_points;
    for (int i = 0; i <= n_points; ++i) {
        double x = -1.0 + i * h;
        double diff = f(x) - p(x, c);
        sse += diff * diff;
    }
    return sse * h;
}

int main() {
    int degree = 2;
    int n_points = 10000;

    vector<double> c(degree + 1);

    for (int k = 0; k <= degree; ++k) {
        // integral of f(x) * P_k(x) from -1 to 1
        double numerator = trapezoidal([k](double x) {
            return f(x) * legendre(k, x);
        }, -1, 1, n_points);

        // integral of P_k(x) * P_k(x) from -1 to 1
        double denominator = trapezoidal([k](double x) {
            double val = legendre(k, x);
            return val * val;
        }, -1, 1, n_points);

        c[k] = numerator / denominator;
    }

    cout << "legendre polynomial approximation:\n";
    for (int k = 0; k <= degree; ++k) {
        cout << fabs(c[k]) << " * ";
        switch (k)
        {
        case 0:
            cout << "1 ";
            if(c[k+1] < 0) {
                cout << "- ";
            }
            else {
                cout << "+ ";
            }
            break;
        case 1:
            cout << "x ";
            if(c[k+1] < 0) {
                cout << "- ";
            }
            else {
                cout << "+ ";
            }
            break;
        case 2:
            cout << "(x^2 - 1/3) ";
            break;
        default:
            break;
        }
    }
    cout << "\n";

    // compute the sum of squared errors (SSE)
    double sse = compute_sse(c, n_points);
    cout << "Sum of Squared Errors (SSE):" << sse;

    return 0;
}
