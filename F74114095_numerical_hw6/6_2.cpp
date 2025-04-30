#include <iostream>
#include <iomanip>
#include <cmath>
using namespace std;

const int N = 4;

void printMatrix(double a[N][N]) {
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j)
            cout << setw(10) << fixed << setprecision(4) << a[i][j] << " ";
        cout << endl;
    }
}

int main() {
    double A[N][N] = {
        {4, 1, -1, 0},
        {1, 3, -1, 0},
        {-1, -1, 6, 2},
        {0, 0, 2, 5}
    };

    double I[N][N] = {0};
    for (int i = 0; i < N; ++i)
        I[i][i] = 1;

    // Gauss-Jordan Elimination
    for (int i = 0; i < N; ++i) {
        double maxA = fabs(A[i][i]);
        int maxRow = i;
        for (int k = i + 1; k < N; ++k) {
            if (fabs(A[k][i]) > maxA) {
                maxA = fabs(A[k][i]);
                maxRow = k;
            }
        }

        if (maxRow != i) {
            for (int j = 0; j < N; ++j) {
                swap(A[i][j], A[maxRow][j]);
                swap(I[i][j], I[maxRow][j]);
            }
        }


        double diag = A[i][i];
        if (fabs(diag) < 1e-12) {
            cout << "matrix can't be inversed" << endl;
            return 0;
        }
        for (int j = 0; j < N; ++j) {
            A[i][j] /= diag;
            I[i][j] /= diag;
        }

        for (int k = 0; k < N; ++k) {
            if (k == i) continue;
            double factor = A[k][i];
            for (int j = 0; j < N; ++j) {
                A[k][j] -= factor * A[i][j];
                I[k][j] -= factor * I[i][j];
            }
        }
    }

    cout << "A's inverse matrix:" << endl;
    printMatrix(I);

    return 0;
}
