#include <iostream>
#include <cmath>

using namespace std;

// Define the function for the first part of the problem
double f(double x) {
    return sin(x) * pow(x, -0.25);
}

// Define the function for the second part of the problem
double g(double t) {
    return sin(1/t) * pow(t, 4);
}

// Simpson's Rule implementation
double composite_simpson_rule(double a, double b, double h, char func_type = 'f') {
    double integral = 0;
    if(func_type == 'f') {
        // integral += f(a+1e-6);
        // integral += f(b);
        // for(double x = a+h; x < b; x += 2*h){
        //     integral += 4*f(x);
        // }
        // for(double x = a+2*h; x < b-h; x += 2*h){
        //     integral += 2*f(x);
        // }
        integral = 7*f(a+1e-6) + 32*f(a+h) + 12*f(a+2*h) + 32*f(a+3*h) + 7*f(b);
    } else {
        // integral += g(a+1e-6);
        // integral += g(b);
        // for(double x = a+h; x < b; x += 2*h){
        //     integral += 4*g(x);
        // }
        // for(double x = a+2*h; x < b-h; x += 2*h){
        //     integral += 2*g(x);
        // }
        integral = 7*g(a+1e-6) + 32*g(a+h) + 12*g(a+2*h) + 32*g(a+3*h) + 7*g(b);
    }
    
    integral = integral * (2*h / 45);
    return integral;
}


int main() {
    double a = 0, b = 1;
    double h1 = (b-a)/4;
    double integral = 0;
    integral = composite_simpson_rule(a, b, h1, 'f');
    cout << "The integral of f(x) = sin(x)*x^(-0.25) using Simpson's rule: " << integral << endl;
    //--------------------------------- 
    double c = 0, d = 1;
    double h2 = (d-c)/4;
    double integral2 = 0;
    integral2 = composite_simpson_rule(c, d, h2, 'g');
    cout << "The integral of g(t) = sin(1/t)*t^4 using Simpson's rule: " << integral2 << endl;
    

    return 0;
}
