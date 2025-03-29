#include<iostream>
#include<cmath>
using namespace std;

double Lagrange(double x[], double y[], int n, double x0) {
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

double ErrorBound(double x[], int n, double x0, double (*f_derivative)(double)) {
    double M = 0;
    for (int i = 0; i < n; i++) {
        M = max(M, fabs(f_derivative(x[i])));
    }
    double product = 1;
    for (int i = 0; i < n; i++) {
        product *= (x0 - x[i]);
    }
    return (M * fabs(product)) / tgamma(n + 1); // n! = tgamma(n+1)
}

double cos_derivative(double x) {
    return -sin(x);  // Derivative of cos(x)
}

int main() {
    double x[] = {0.698, 0.733, 0.768, 0.803};
    double y[] = {0.7661, 0.7432, 0.7193, 0.6946};
    double x0 = 0.75;
    double Y = cos(0.75);
    
    for(int d = 1; d <= 4; d++) {
        double y0 = Lagrange(x, y, d, x0);
        printf("Degree: %d, cos(0.75) = %lf\n", d, y0);
        
        double error_bound = ErrorBound(x, d, x0, cos_derivative);
        printf("Error Bound: %.10f\n", error_bound);
        printf("\n");
    }

    return 0;
}
