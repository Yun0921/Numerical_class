#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
using namespace std;

const double PI = acos(-1);
const int nx = 10; // x: 0 ~ π, step = 0.1π, total 11 points
const int ny = 5;  // y: 0 ~ π/2, step = 0.1π, total 6 points
const double h = PI / nx;      // h = k = 0.1π
const double k = PI / (2 * ny);
const double alpha = (h * h) / (k * k);
const double tol = 1e-6;
const int MAX_ITER = 10000;

// 高斯消去法解線性方程組
vector<double> gaussianElimination(vector<vector<double>> A, vector<double> b) {
    int n = A.size();
    
    // 前向消去
    for (int i = 0; i < n; i++) {
        // 找主元
        int maxRow = i;
        for (int k = i + 1; k < n; k++) {
            if (abs(A[k][i]) > abs(A[maxRow][i])) {
                maxRow = k;
            }
        }
        
        // 交換行
        if (maxRow != i) {
            swap(A[i], A[maxRow]);
            swap(b[i], b[maxRow]);
        }
        
        // 消去
        for (int k = i + 1; k < n; k++) {
            double factor = A[k][i] / A[i][i];
            b[k] -= factor * b[i];
            for (int j = i; j < n; j++) {
                A[k][j] -= factor * A[i][j];
            }
        }
    }
    
    // 回代
    vector<double> x(n);
    for (int i = n - 1; i >= 0; i--) {
        double sum = 0.0;
        for (int j = i + 1; j < n; j++) {
            sum += A[i][j] * x[j];
        }
        x[i] = (b[i] - sum) / A[i][i];
    }
    
    return x;
}

double f(double x, double y) {
    return x * y;
}

int main() {
    const int interior_points = (nx-1) * (ny-1);  // 內部點的數量
    vector<vector<double>> A(interior_points, vector<double>(interior_points, 0.0));
    vector<double> F(interior_points, 0.0);
    vector<vector<double>> u(nx + 1, vector<double>(ny + 1, 0.0));

    // 設定邊界條件
    for (int j = 0; j <= ny; ++j) {
        double y = j * k;
        u[0][j] = cos(y);         // u(0, y) = cos(y)
        u[nx][j] = -cos(y);       // u(π, y) = -cos(y)
    }

    for (int i = 0; i <= nx; ++i) {
        double x = i * h;
        u[i][0] = cos(x);         // u(x, 0) = cos(x)
        u[i][ny] = 0.0;           // u(x, π/2) = 0
    }

    // 建立係數矩陣 A 和向量 F
    int row = 0;
    for (int j = 1; j < ny; ++j) {
        for (int i = 1; i < nx; ++i) {
            double x = i * h;
            double y = j * k;
            
            // 當前點的係數
            A[row][row] = -2.0 * (1.0 + alpha);
            
            // 計算相鄰點的係數
            if (i > 1)  // 左點
                A[row][row-1] = 1.0;
            if (i < nx-1)  // 右點
                A[row][row+1] = 1.0;
            if (j > 1)  // 下點
                A[row][row-(nx-1)] = alpha;
            if (j < ny-1)  // 上點
                A[row][row+(nx-1)] = alpha;
            
            // 計算右側向量 F
            F[row] = h * h * f(x, y);
            
            // 加入邊界條件的影響
            if (i == 1)
                F[row] -= u[0][j];
            if (i == nx-1)
                F[row] -= u[nx][j];
            if (j == 1)
                F[row] -= alpha * u[i][0];
            if (j == ny-1)
                F[row] -= alpha * u[i][ny];
                
            row++;
        }
    }

    // 解線性方程組
    vector<double> solution = gaussianElimination(A, F);

    // 將解放回 u 矩陣
    row = 0;
    for (int j = 1; j < ny; ++j) {
        for (int i = 1; i < nx; ++i) {
            u[i][j] = solution[row++];
        }
    }

    cout << setw(6) << "y/x";
    for (int i = 0; i <= nx; ++i) {
        double r = 0 + i * h;
        cout << setw(12) << r;
    }
    cout << endl;

    cout << fixed << setprecision(4);
    for (int j = ny; j >= 0; --j) {
        double theta = j * k;
        cout << setw(6) << theta;
        for (int i = 0; i <= nx; ++i) {
            cout << setw(12) << u[i][j];
        }
        cout << endl;
    }

    return 0;
}
