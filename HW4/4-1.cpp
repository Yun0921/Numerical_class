#include<iostream>
#include<cmath>

using namespace std;

double f(double x){
    return exp(x)*sin(4*x);
}
double composite_trapezoidal_rule(double a, double b, double h){
    double integral = 0;
    int cnt = 0;
    integral += f(a);
    integral += f(b);
    for(double x = a+h; x < b; x += h){
        integral += 2*f(x);
    }
    integral *= h / 2;
    return integral;
}

double composite_simpson_rule(double a, double b, double h){
    double integral = 0;
    int cnt = 0;
    for(double x = a; x < b; x += h){
        integral += (f(x) + 4*f(x+h/2) + f(x+h));
        
    }
    integral = integral * h / 6;
    return integral;
}

double composite_midpoint_rule(double a, double b, double h){
    double integral = 0;
    for(double x = a+h; x < b; x += 2*h){
        integral += f(x);
    }
    integral = integral * 2*h;
    return integral;
}

int main(){
    int a = 0, b = 2;
    double h = 0.1;
    double intergral = 0;
    //cout << "check" << endl;
    intergral = composite_trapezoidal_rule(a, b, h);
    cout << "The integral of f(x) = exp(x)*sin(4*x) using Trapezopidal_rule: " << intergral << endl;
    intergral = composite_simpson_rule(a, b, h);
    cout << "The integral of f(x) = exp(x)*sin(4*x) using Simpson_rule: " << intergral << endl;
    intergral = composite_midpoint_rule(a, b, h);
    cout << "The integral of f(x) = exp(x)*sin(4*x) using Midpoint_rule: " << intergral << endl;

    return 0;
}