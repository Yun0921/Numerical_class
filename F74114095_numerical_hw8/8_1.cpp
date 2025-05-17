#include <iostream>
#include <vector>
#include <cmath>
#include <numeric>
using namespace std;

// degree 2 polynomial fitting
void poly2_fit(const vector<double>& x, const vector<double>& y, vector<double>& coeffs) {
    int n = x.size();
    double Sx = 0, Sx2 = 0, Sx3 = 0, Sx4 = 0;
    double Sy = 0, Sxy = 0, Sx2y = 0;

    for (int i = 0; i < n; ++i) {
        double xi = x[i], yi = y[i];
        double xi2 = xi*xi;
        double xi3 = xi2*xi;
        double xi4 = xi3*xi;
        Sx += xi;
        Sx2 += xi2;
        Sx3 += xi3;
        Sx4 += xi4;
        Sy += yi;
        Sxy += xi * yi;
        Sx2y += xi2 * yi;
    }

    // | n    Sx   Sx2  |   |a0|   |Sy  |
    // | Sx   Sx2  Sx3  | * |a1| = |Sxy |
    // | Sx2  Sx3  Sx4  |   |a2|   |Sx2y|

    // use Cramer's rule to solve the system of equations
    double D = n*(Sx2*Sx4 - Sx3*Sx3) - Sx*(Sx*Sx4 - Sx3*Sx2) + Sx2*(Sx*Sx3 - Sx2*Sx2);
    double Da0 = Sy*(Sx2*Sx4 - Sx3*Sx3) - Sx*(Sxy*Sx4 - Sx3*Sx2y) + Sx2*(Sxy*Sx3 - Sx2*Sx2y);
    double Da1 = n*(Sxy*Sx4 - Sx3*Sx2y) - Sy*(Sx*Sx4 - Sx3*Sx2) + Sx2*(Sx*Sx2y - Sxy*Sx2);
    double Da2 = n*(Sx2*Sx2y - Sxy*Sx3) - Sx*(Sx*Sx2y - Sxy*Sx2) + Sy*(Sx*Sx3 - Sx2*Sx2);

    coeffs[0] = Da0 / D; // a0
    coeffs[1] = Da1 / D; // a1
    coeffs[2] = Da2 / D; // a2
}

// compute sum of squared errors (SSE)
double sse(const vector<double>& y_true, const vector<double>& y_pred) {
    double sum = 0;
    int n = y_true.size();
    for (int i = 0; i < n; ++i) {
        double diff = y_true[i] - y_pred[i];
        sum += diff * diff;
    }
    return sum;
}

// exponential fitting y = b * exp(a*x), transform to ln(y) = a*x + ln(b) linear fitting
void exp_fit(const vector<double>& x, const vector<double>& y, double& a, double& b) {
    int n = x.size();
    vector<double> ln_y(n);
    for (int i = 0; i < n; ++i) {
        ln_y[i] = log(y[i]);
    }

    double Sx = accumulate(x.begin(), x.end(), 0.0);
    double Sy = accumulate(ln_y.begin(), ln_y.end(), 0.0);
    double Sx2 = 0, Sxy = 0;
    for (int i = 0; i < n; ++i) {
        Sx2 += x[i]*x[i];
        Sxy += x[i]*ln_y[i];
    }

    // | Sx    n  |   | a  |   |Sy  |
    // | Sx2   Sx | * |ln_b| = |Sxy |

    // use Cramer's rule to solve the system of equations
    a = (n * Sxy - Sx * Sy) / (n * Sx2 - Sx * Sx);
    double ln_b = (Sy - a * Sx) / n;
    b = exp(ln_b);
}

// exponential fitting y = b * x^n, transform to ln(y) = n*ln(x) + ln(b) linear fitting
void power_fit(const vector<double>& x, const vector<double>& y, double& n, double& b) {
    int size = x.size();
    vector<double> ln_x(size), ln_y(size);
    for (int i = 0; i < size; ++i) {
        ln_x[i] = log(x[i]);
        ln_y[i] = log(y[i]);
    }

    double Sx = accumulate(ln_x.begin(), ln_x.end(), 0.0);
    double Sy = accumulate(ln_y.begin(), ln_y.end(), 0.0);
    double Sx2 = 0, Sxy = 0;
    for (int i = 0; i < size; ++i) {
        Sx2 += ln_x[i]*ln_x[i];
        Sxy += ln_x[i]*ln_y[i];
    }

    // | Sx    size|   | n  |   |Sy  |
    // | Sx2   Sx  | * |ln_b| = |Sxy |

    // use Cramer's rule to solve the system of equations
    n = (size * Sxy - Sx * Sy) / (size * Sx2 - Sx * Sx);
    double ln_b = (Sy - n * Sx) / size;
    b = exp(ln_b);
}

int main() {
    vector<double> x = {4.0, 4.2, 4.5, 4.7, 5.1, 5.5, 5.9, 6.3};
    vector<double> y = {102.6, 113.2, 130.1, 142.1, 167.5, 195.1, 224.9, 256.8};

    // a.
    vector<double> coeffs(3, 0);
    poly2_fit(x, y, coeffs);
    vector<double> y_poly2(x.size());
    for (size_t i = 0; i < x.size(); ++i) {
        y_poly2[i] = coeffs[2]*x[i]*x[i] + coeffs[1]*x[i] + coeffs[0];
    }
    double error_poly2 = sse(y, y_poly2);
    cout << "a. degree 2: a0 = " << coeffs[0] << ", a1 = " << coeffs[1] << ", a2 = " << coeffs[2] << endl;
    cout << "a. error: " << error_poly2 << endl;

    // b.
    double a_exp, b_exp;
    exp_fit(x, y, a_exp, b_exp);
    vector<double> y_exp(x.size());
    for (size_t i = 0; i < x.size(); ++i) {
        y_exp[i] = b_exp * exp(a_exp * x[i]);
    }
    double error_exp = sse(y, y_exp);
    cout << "b. a = " << a_exp << ", b = " << b_exp << endl;
    cout << "b. error: " << error_exp << endl;

    // c.
    double n_power, b_power;
    power_fit(x, y, n_power, b_power);
    vector<double> y_power(x.size());
    for (size_t i = 0; i < x.size(); ++i) {
        y_power[i] = b_power * pow(x[i], n_power);
    }
    double error_power = sse(y, y_power);
    cout << "c. n = " << n_power << ", b = " << b_power << endl;
    cout << "c. error: " << error_power << endl;

    return 0;
}
