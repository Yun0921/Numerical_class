#include <iostream>
#include <cmath>
using namespace std;

double Lagrange(double x[], double y[], int n, double x0) {
    double result = 0;
    for (int i = 0; i <= n; i++) {
        double temp = 1;
        for (int j = 0; j <= n; j++) {
            if (j != i) {
                temp *= (x0 - x[j]) / (x[i] - x[j]);
            }
        }
        result += temp * y[i];
    }
    return result;
}

double ErrorBound(double x[], int n, double x0, double (*f_derivative)(double, int)) {
    //set C  as the last one x[i]
    double M = 0;
    for (int i = 0; i <= n; i++) {  // Iterate over all nodes
        M = max(M, fabs(f_derivative(x[i], n+1))); // Compute (n+1)th derivative
    }
    double product = 1;
    for (int i = 0; i <= n; i++) {
        product *= (x0 - x[i]);
    }
    return abs(M * fabs(product)) / tgamma(n + 2); // n+1! = tgamma(n+2)
}

double cos_derivative(double x, int order) {
    // Derivatives of cos(x) for higher orders
    if (order == 1) {
        return -sin(x);  // 1st derivative
    } else if (order == 2) {
        return -cos(x);  // 2nd derivative
    } else if (order == 3) {
        return sin(x);   // 3rd derivative
    } else if (order == 4) {
        return cos(x);   // 4th derivative
    }
    return 0;
}

int main() {
    double x[] = {0.698, 0.733, 0.768, 0.803};
    double y[] = {0.7661, 0.7432, 0.7193, 0.6946};
    double x0 = 0.75;
    double Y = cos(0.75);  // True value of cos(0.75)
    
    for (int d = 1; d <= 3; d++) {
        double y0 = Lagrange(x, y, d, x0);
        printf("Degree: %d, Lagrange Approximation: %lf\n", d, y0);
        
        double error_bound = ErrorBound(x, d, x0, cos_derivative);
        printf("Error Bound: %.10f\n", error_bound);
        printf("\n");
    }

    return 0;
}
