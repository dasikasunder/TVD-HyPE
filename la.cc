/*
 * la.cc
 *      Author: sunder
 */

#include "la.hh"

//----------------------------------------------------------------------------
// Find the minimum value in a vector V of size N 
//----------------------------------------------------------------------------

double MinVal(const Vector& V, const int& N) {

    double min = V[0]; 
    
    for (int i = 1; i < N; ++i) {
        if (V[i] < min)
            min = V[i]; 
    }
    
    return min; 
}

//----------------------------------------------------------------------------
// Find the maximum value in a vector V of size N 
//---------------------------------------------------------------------------

double MaxVal(const Vector& V, const int& N) {

    double max = V[0]; 
    
    for (int i = 1; i < N; ++i) {
        if (V[i] > max)
            max = V[i]; 
    }
    
    return max; 
}

//----------------------------------------------------------------------------
// Given two vectors V1 and V2 of size N, find the L2 residue between them  
//---------------------------------------------------------------------------

double residue(const Vector& V1, const Vector& V2, const int& N) {
    
    double sum = 0.0;
    double diff; 
    
    for (int i = 0; i < N; ++i) {
        diff = V1[i] - V2[i];
        sum += diff*diff; 
    }
    
    return std::sqrt(sum); 
}

//----------------------------------------------------------------------------
// Matrix-Vector Multiplication
// Does the operation: y := A*x
// A -> (m x n) matrix
// x -> (n) vector
// y -> (m) vector
//----------------------------------------------------------------------------

void MatVecMult(const Matrix& A, const Vector& x, Vector& y, int m, int n) {
    for (int i = 0; i < m; i++ ) {
        y[i]= 0.0;
        for (int j = 0; j < n; j++ )
            y[i]+= A[i][j]*x[j];
    }
}

//----------------------------------------------------------------------------
// Matrix-Matrix Multiplication
// Does the operation: C := A*B
// A -> (n x m) matrix
// B -> (m x p) matrix
// C -> (n x p) matrix
//----------------------------------------------------------------------------

void MatMatMult(const Matrix& A, const Matrix& B, Matrix& C, int n, int m, int p) {
    
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < p; ++j) {
            C[i][j] = 0.0; 
            for (int k = 0; k < m; ++k) {
                C[i][j] += A[i][k]*B[k][j]; 
            }
        }
    }
}