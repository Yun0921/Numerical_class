#include<iostream>
#include<cmath>

using namespace std;

double f(double x, double y){
    return 2*y*sin(x) + cos(x)*cos(x);
}

double C1(double x){
    return (sin(x) + cos(x))/2;
}
double C2(double x){
    return (cos(x) - sin(x))/2;
}

double gaussianQuadrature3(double a, double b){
    double x[3] = {-sqrt(3.0/5), 0, sqrt(3.0/5)};
    double w[3] = {5.0/9, 8.0/9, 5.0/9};

    double integral = 0;
    for(int i = 0; i < 3; i++){
        double inner_integral = 0;
        int tmp = (b-a)/2.0*x[i] + (a+b)/2.0;
        for(int j = 0; j < 3; j++){ 
            inner_integral += w[j]*f(tmp, C1(tmp)+C2(tmp)*x[j]);
        }
        integral += w[i]*C2(tmp)*inner_integral;
        //integral += w[i]*inner_integral;
    }
    return integral*(b-a)/2;
}



int main(){
    double a = 0;
    double b = M_PI/4;
    double integral = gaussianQuadrature3(a, b);
    cout << "The integral of f(x) = 2y*sin(x) + cos(x)*cos(x) using Gaussian Quadrature: " << integral << endl;



    return 0;
}