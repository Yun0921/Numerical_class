#include <iostream>
#include <cmath>
using namespace std;

const int N = 4;

void printMatrix(double a[N][N+1]) {
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j <= N; ++j) {
            if(fabs(a[i][j]) < 1e-6) {
                printf("%10.4f ", 0.0);
            }
            else{
                printf("%10.4f ", a[i][j]);
            }
        }
        printf("\n");
    }
    printf("------------------------------------------------------\n");
}

int main() {
    // Augmented matrix [A|b]
    double a[N][N+1] = {
        {1.19, 2.11, -100.0, 1.0, 1.12},
        {14.2, -0.112, 12.2, -1.0, 3.44},
        {0.0, 100.0, -99.9, 1.0, 2.15},
        {15.3, 0.110, -13.1, -1.0, 4.16}
    };

    // Partial Pivoting and Forward Elimination
    for (int k = 0; k < N-1; ++k) {
        // Find the row with the largest pivot
        int maxRow = k;
        double maxVal = fabs(a[k][k]);
        for (int i = k+1; i < N; ++i) {
            if (fabs(a[i][k]) > maxVal) {
                maxVal = fabs(a[i][k]);
                maxRow = i;
            }
        }
        // Swap rows if needed
        if (maxRow != k) {
            for (int j = 0; j <= N; ++j) {
                swap(a[k][j], a[maxRow][j]);
            }
        }

        // Elimination
        for (int i = k+1; i < N; ++i) {
            double factor = a[i][k] / a[k][k];
            for (int j = k; j <= N; ++j) {
                a[i][j] -= factor * a[k][j];
            }
        }

        printMatrix(a); // Print the matrix after each elimination step
    }

    // Back Substitution
    double x[N];
    for (int i = N-1; i >= 0; --i) {
        x[i] = a[i][N];
        for (int j = i+1; j < N; ++j) {
            x[i] -= a[i][j] * x[j];
        }
        x[i] /= a[i][i];
    }

    // Output the result
    cout << "result:" << endl;
    for (int i = 0; i < N; ++i) {
        cout << "x" << i+1 << " = " << x[i] << endl;
    }

    return 0;
}
