#include<iostream>
#include<cmath>

using namespace std;

double f(double x){
    return pow(x, 2)*log(x);
}

double gaussianQuadrature3(double a, double b){
    double x1 = -sqrt(3.0/5);
    double x2 = 0;
    double x3 = sqrt(3.0/5);
    double w1 = 5.0/9;
    double w2 = 8.0/9;
    double w3 = 5.0/9;
    double integral = 0;
    integral = w1*f((b-a)/2*x1 + (a+b)/2) + w2*f((b-a)/2*x2 + (a+b)/2) + w3*f((b-a)/2*x3 + (a+b)/2);
    integral *= (b-a)/2;
    return integral;
}

double gaussianQuadrature4(double a, double b){
    double x1 = -sqrt(3.0/7 + 2.0/7*sqrt(6.0/5));
    double x2 = -sqrt(3.0/7 - 2.0/7*sqrt(6.0/5));
    double x3 = sqrt(3.0/7 - 2.0/7*sqrt(6.0/5));
    double x4 = sqrt(3.0/7 + 2.0/7*sqrt(6.0/5));
    //printf("%f %f %f %f\n", x1, x2, x3, x4);
    double w1 = (18 - sqrt(30))/36;
    double w2 = (18 + sqrt(30))/36;
    double integral = 0;
    integral = w1*f((b-a)/2*x1 + (a+b)/2) + w2*f((b-a)/2*x2 + (a+b)/2) + w2*f((b-a)/2*x3 + (a+b)/2) + w1*f((b-a)/2*x4 + (a+b)/2);
    integral *= (b-a)/2;
    return integral;
}

int main(){
    double a = 1, b = 1.5;
    int n1 = 3;
    int n2 = 4;

    double exact_integral = 0.1922593577;
    
    double result1 = gaussianQuadrature3(a, b);
    double result2 = gaussianQuadrature4(a, b);
    printf("the integral of f(x) = x^2*log(x) using Gaussian Quadrature with n = 3: %.10f\n", result1);
    printf("Absolute error of n = 3: %.10f\n", abs(exact_integral - result1));
    printf("the integral of f(x) = x^2*log(x) using Gaussian Quadrature with n = 4: %.10f\n", result2);
    printf("Absolute error of n = 4: %.10f\n", abs(exact_integral - result2));
    //cout << "The integral of f(x) = x^2*log(x) using Gaussian Quadrature with n = 3: " << result1 << endl;
    //cout << "The integral of f(x) = x^2*log(x) using Gaussian Quadrature with n = 4: " << result2 << endl;



    return 0;
}