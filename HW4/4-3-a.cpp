#include<iostream>
#include<cmath>

using namespace std;

double f(double x, double y){
    return 2*y*sin(x) + cos(x)*cos(x);
}

double k(double x){
    return (cos(x)-sin(x))/4;
}

double F(double x){
    double y0 = sin(x);
    double y1 = sin(x)+k(x);
    double y2 = sin(x)+2*k(x);
    double y3 = sin(x)+3*k(x);
    double y4 = sin(x)+4*k(x);
    double Fx = k(x)*(7*f(x, y0) + 32*f(x, y1) + 12*f(x, y2) + 32*f(x, y3) + 7*f(x, y4));

    return Fx;
}


int main(){
    double x0 = 0;
    double x1 = M_PI/16;
    double x2 = M_PI/8;
    double x3 = 3*M_PI/16;
    double x4 = M_PI/4;
    double h = (x4-x0)/4;
    double integral = 0;
    integral = pow(2.0/45, 2) * M_PI/16 * (7*F(x0) + 32*F(x1) + 12*F(x2) + 32*F(x3) + 7*F(x4));
    cout << "The integral of f(x) = 2y*sin(x) + cos(x)*cos(x) using Simpson's rule: " << integral << endl;

    double exact_integral = 0.511845;
    cout << "Absolute error: " << abs(exact_integral - integral) << endl;
    

    return 0;
}