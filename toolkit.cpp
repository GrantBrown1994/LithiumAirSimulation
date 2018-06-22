/*
Grant Brown
Content:
1.  solveGauss  Uses the Gaussian Elimination technique to solve matrices.
                https://martin-thoma.com/solving-linear-equations-with-gaussian-elimination/
2.  exportdata  Data export through fstream.
3.  printAb     Prints matrix with rhs denoted with '|'. LHS of Ax = b
4.  printA      Prints a matrix.
5.  printx      Prints a vector.
*/
#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <algorithm>
#include <stdlib.h>
#include <stdio.h>

using namespace std;

vector<double> solveGauss(vector< vector<double> > AugMat)
{
    int n = AugMat.size();

    for (int i=0; i<n; i++) {
        // Search for maximum in this column
        double maxEl = abs(AugMat[i][i]);
        int maxRow = i;
        for (int k=i+1; k<n; k++) {
            if (abs(AugMat[k][i]) > maxEl) {
                maxEl = abs(AugMat[k][i]);
                maxRow = k;
            }
        }

        // Swap maximum row with current row (column by column)
        for (int k=i; k<n+1;k++) {
            double tmp = AugMat[maxRow][k];
            AugMat[maxRow][k] = AugMat[i][k];
            AugMat[i][k] = tmp;
        }

        // Make all rows below this one 0 in current column
        for (int k=i+1; k<n; k++) {
            double c = -AugMat[k][i]/AugMat[i][i];
            for (int j=i; j<n+1; j++) {
                if (i==j) {
                    AugMat[k][j] = 0;
                } else {
                    AugMat[k][j] += c * AugMat[i][j];
                }
            }
        }
    }

    // Solve equation Ax=b for an upper triangular matrix A
    vector<double> v(n);
    //cout << n-1 << endl;
    for (int i=n-1; i>=0; i--) {
        v[i] = AugMat[i][n]/AugMat[i][i];
        for (int k=i-1;k>=0; k--) {
            AugMat[k][n] -= AugMat[k][i] * v[i];
        }
    }
    return v;
}
void exportdata(const vector<double> v, string filename)
{
    ofstream output;
    output.clear();
	output.open(filename);
	
	int i = 0, n = v.size();
	for(auto itr: v)
	{
	    output << i << ", " << itr;
	    if(i < n-1)
	        output << endl;
	    i++;
	}
	output.close();
// 	int n = v.size();
// 	for(int i=0; i<n; i++)
// 	{
//         output << i << ", "<< v[i];
//         if(i < n - 1)
//         	output << "\n";
//     }
//     output.close();
}
void printAug(const vector< vector<double> > & AugMat)
{
    int n = AugMat.size();
    for (int i=0; i<n; ++i) {
        for (int j=0; j<n+1; ++j) {
            cout << AugMat[i][j] << " ";
            if (j == n-1) {
                cout << "| ";
            }
        }
        cout << "\n";
    }
    cout << endl;
}
void printVec(const vector <double> v, int cols = 0)
{
    int n = v.size();
    cout << "Results: \n";
    if (cols > 0)
    {
        for(int i = 0; i < n; ++i)
        {
            cout << fixed << setprecision(6) << right << setw(10);
            cout << v[i] << " ";
            if((i+1)%cols == 0)
                cout << "\n";
        }
    }
    else
    {
        for(auto itr: v)
            cout << itr << endl;
    }
    
}
double norm(vector <double> v)
{
	double r = 0;
	//int n = v.size();
	for(auto itr: v)
		r += fabs(itr);
	return r;
}
