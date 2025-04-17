#include<iostream>
#include<iomanip>
#include<cmath>

using namespace std;

//f(t, y)
double f(double t, double y) {
    double yt = y / t;
    return 1.0 + yt + pow(yt, 2);
}

//df_dt
double df_dt(double t, double y) {
    double yt = y / t;
    return -y / (t * t) + 2 * yt * (-y / (t * t));
}

//df_dy
double df_dy(double t, double y) {
    double yt = y / t;
    return 1.0 / t + 2.0 * y / pow(t, 2);
}

//exact value
double exact(double t) {
    return t * tan(log(t));
}

//Euler method
void euler(double t0, double y0, double h, int step) {
    double t = t0;
    double y = y0;

    cout << "t\ty\tvalue\tabsolute error" << endl;
    cout << fixed << setprecision(4);
    for (int i = 0; i <= step; i++) {
        double exact_value = exact(t);
        double absolute_error = exact_value - y;
        printf("%.1f\t%.4f\t%.4f\t%.6f\n", t, y, exact_value, absolute_error);
        y += h * f(t, y);
        t += h;
    }
}

//y''
double y2(double t, double y) {
    return df_dt(t, y)  + df_dy(t, y) * f(t, y);
}

//Taylor method order 2
void taylor2(double t0, double y0, double h, int step) {
    double t = t0;
    double y = y0;

    cout << "t\ty\tvalue\tabsolute error" << endl;
    cout << fixed << setprecision(4);
    for (int i = 0; i <= step; i++) {
        double exact_value = exact(t);
        double absolute_error = exact_value - y;
        printf("%.1f\t%.4f\t%.4f\t%.6f\n", t, y, exact_value, absolute_error);
        double y_2 = y2(t, y);
        y += h * (f(t, y) + (h / 2) * y_2);
        t += h;
    }
}


int main() {
    double t0 = 1.0, y0 = 0.0, h = 0.1;
    int step = 10;
    printf("Euler method\n");
    euler(t0, y0, h, step);
    printf("\n");
    printf("Taylor method order 2\n");
    taylor2(t0, y0, h, step);
    return 0;
}