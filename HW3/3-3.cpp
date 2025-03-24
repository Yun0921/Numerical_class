#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

double table[10][10];
vector<double> x_values = {0, 3, 5, 8, 13};
//vector<double> f_values = {0, 225, 383, 623 ,993}; 
vector<double> f_values = {0, 200, 375, 620, 990}; 
vector<double> df_values = {75, 77, 80, 74, 72};


void Divided_difference_table(){
    for(int i = 0; i < 10; i++){
        for(int j = 0; j < 10; j++){
            table[i][j] = 0;
        }
    }
    for(int i = 0; i < 10; i++){
        table[i][0] = f_values[i/2];
        if(i % 2 == 1){
            table[i][1] = df_values[i/2];
        }
        else{
            table[i][1] = (table[i][0] - table[i-1][0]) / (x_values[i/2] - x_values[i/2-1]);
        }
    }
    for(int col = 2; col < 10; col++){
        for(int row = col; row < 10; row++){
            table[row][col] = (table[row][col-1] - table[row-1][col-1]) / (x_values[row/2] - x_values[(row-col)/2]);
        }
    }
}

double Hermite(double x0){
    double result = 0;
    for(int i = 0; i < 10; i++){
        double temp = 1;
        for(int j = 0; j < i; j++){
            temp *= (x0 - x_values[j/2]);
        }
        result += temp * table[i][i];
    }
    return result;
}


int main() {
    
    Divided_difference_table();
    Hermite(10);
    cout << "The position at t = 10: " << Hermite(10) << " feet" << endl;

    double d = 0.0001;
    double df = (Hermite(10+d) - Hermite(10-d)) / (2*d);
    cout << "The velocity at t = 10: " << df << " feet/sec" << endl;

    cout << "Polynomial coefficients: ";
    for(int i = 0; i < 10; i++){
        cout << table[i][i] << ","<< " ";
    }

    
    //prnit the table
    // for(int i = 0; i < 10; i++){
    //     for(int j = 0; j < 10; j++){
    //         cout << table[i][j] << " ";
    //     }
    //     cout << endl;
    // }
    

    return 0;
}
