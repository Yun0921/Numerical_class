#include<iostream>
#include<cmath>
#include<vector>
using namespace std;

double x[] = {0.3, 0.4, 0.5, 0.6};
double y[] = {-0.440818, -0.27032, -0.106531, 0.051188};
int n = 4;

double Lagrange(double x0) {
    double result = 0;
    for (int i = 0; i < n; i++) {
        double temp = 1;
        for (int j = 0; j < n; j++) {
            if (j != i) {
                temp *= (x0 - x[j]) / (x[i] - x[j]);
            }
        }
        result += temp * y[i];
    }
    return result;
}

pair<double, int> secant(double x0, double x1, double tolerance, int max_iterations) {
    double x_prev = x0;
    double x_curr = x1;
    double x_next;

    for (int i = 0; i < max_iterations; ++i) {
        x_next = x_curr - Lagrange(x_curr) * (x_curr - x_prev) / (Lagrange(x_curr) - Lagrange(x_prev));

        if (abs(x_next - x_curr) < tolerance) {
            return make_pair(x_next, i + 1);
        }
        x_prev = x_curr;
        x_curr = x_next;
        //printf("Iteration %d: x = %.10f, Lagrange(x) = %.10f\n", i + 1, x_curr, Lagrange(x_curr));
    }
    return make_pair(NAN, max_iterations); // Did not converge
}

int main() {
    double x0_secant = 0;
    double x1_secant = 1;
    double tolerance = 1e-6;
    int max_iterations = 1000;
    vector<double> roots;
    
    cout << "Inverse Interpolation using Secant" << endl;
    cout << "Secant Method (Initial guesses " << x0_secant << ", " << x1_secant << "):" << endl;
    pair<double, int> result_secant = secant(x0_secant, x1_secant, tolerance, max_iterations);
    if (!isnan(result_secant.first)) {
       cout << "Converged to root: " << result_secant.first << " in " << result_secant.second << " iterations" << endl;
        roots.push_back(result_secant.first);

    } else {
      cout << "Did not converge within " << max_iterations << " iterations." << endl;
    }
    //printf("%lf\n", Lagrange(0.6));

    
    return 0;
}