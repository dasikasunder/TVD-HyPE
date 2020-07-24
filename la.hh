/*
 * la.hh
 *      Author: sunder
 */

#ifndef LA_HH_
#define LA_HH_

// Functions related to linear algebra 

# include "headers.hh"

double MinVal(const Vector&, const int&);
double MaxVal(const Vector&, const int&);
double residue(const Vector&, const Vector&, const int&);
void MatVecMult(const Matrix&, const Vector&, Vector&, int, int);
void MatMatMult(const Matrix&, const Matrix&, Matrix&, int, int, int);

//----------------------------------------------------------------------------
// LU decomposition class
//----------------------------------------------------------------------------

class LUdcmp {
    int n;                        // Size of the matrix
    Matrix LU;                    // Stores the decomposition
    std::vector<int> indx;        // Stores the permutation
    double d;                     // Used in determinant calculation
    public:
LUdcmp(const Matrix&, int);                         // Constructor
    void solve(const Vector&, Vector&) const;       // Solve for single RHS
    void solve(const Matrix&, int, Matrix&) const;  // Solve for multiple RHS
    void inverse(Matrix&) const;                    // Invert the matrix A
    double det() const;
};

//----------------------------------------------------------------------------
// QR decomposition class
//----------------------------------------------------------------------------

class QRdcmp {

    Matrix QR_;
    int m, n;
    Vector Rdiag;

public:
    QRdcmp(); 
    QRdcmp(Matrix&, int, int);
    void initialize(Matrix&, int, int);
    bool is_full_rank() const;
    void get_R(Matrix&) const;
    void get_Q(Matrix&) const;
    void solve(const Vector&, Vector&) const;
    void solve(const Matrix&, int, Matrix&) const;
};


#endif /* LA_HH_ */