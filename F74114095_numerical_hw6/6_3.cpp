#include <iostream>
#include <vector>
using namespace std;

int main() {
    const int n = 4;
    double a[n] = {0, -1, -1, -1};
    double b[n] = {3, 3, 3, 3};
    double c[n] = {-1, -1, -1, 0};
    double d[n] = {2, 3, 4, 1};

    double l[n], u[n];
    l[0] = b[0];
    u[0] = c[0] / l[0];
    for (int i = 1; i < n; ++i) {
        l[i] = b[i] - a[i] * u[i-1];
        if (i < n-1)
            u[i] = c[i] / l[i];
    }

    //Ly = d
    double y[n];
    y[0] = d[0] / l[0];
    for (int i = 1; i < n; ++i) {
        y[i] = (d[i] - a[i] * y[i-1]) / l[i];
    }

    //Ux = y
    double x[n];
    x[n-1] = y[n-1];
    for (int i = n-2; i >= 0; --i) {
        x[i] = y[i] - u[i] * x[i+1];
    }

    cout << "result:" << endl;
    for (int i = 0; i < n; ++i) {
        cout << "x" << i+1 << " = " << x[i] << endl;
    }
    return 0;
}
