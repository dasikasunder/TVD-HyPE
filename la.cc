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


//----------------------------------------------------------------------------
// Constructor for LUdcmp class. The input is a square (n x n) matrix
//----------------------------------------------------------------------------

LUdcmp::LUdcmp(const Matrix& A, int n_) :
	n(n_),
	LU(A),
	indx(n_),
	d(1.0)
{
    const double TINY = 1.0e-40;
    int i, imax, j, k;
    double big, temp;
    double* vv = new double[n];

    imax = 0; 

    // Loop over rows to get the implicit scaling information

    for (i=0; i<n; i++) {

        big = 0.0;

        for (j=0; j<n; j++) {

            if ((temp=std::abs(LU[i][j])) > big) {
                big=temp;
            }
        }

        if (big == 0.0) {
            std::cerr << "Singular matrix in LUdcmp" << std::endl;
            std::exit(EXIT_FAILURE); // No nonzero largest element.
        }

        vv[i] = 1.0/big; // Save the scaling
    }

    for (k = 0; k < n; k++) {

        big = 0.0;  // Initialize for the search for largest pivot element

        for (i = k; i < n; i++) {

            temp=vv[i]*std::abs(LU[i][k]);

            if (temp > big) {

                big = temp;
                imax = i;
            }

        }

        if (k != imax) {

            for (j = 0; j < n; j++) {

                temp = LU[imax][j];
                LU[imax][j] = LU[k][j];
                LU[k][j] = temp;


            }

            d = -d;
            vv[imax] = vv[k];
        }

        indx[k] = imax;

        if (LU[k][k] == 0.0)
            LU[k][k] = TINY;

        for (i=k+1;i<n;i++)  {

            temp = LU[i][k]/= LU[k][k]; // Divide by pivot element

            for (j=k+1;j<n;j++) {

                LU[i][j] -= temp*LU[k][j];
            }
        }
    }

    delete[] vv;
}

//----------------------------------------------------------------------------
// Solve for a single rhs b, a vector of size n and store the result in x
//----------------------------------------------------------------------------

void LUdcmp::solve(const Vector& b, Vector& x) const {

	int i, ii=0, ip, j;

    double sum;

    for (i = 0; i < n; i++)
    	x[i] = b[i];

    for (i = 0; i < n; i++) {

        ip=indx[i];
        sum=x[ip];
        x[ip]=x[i];

        if (ii != 0) {
            for (j = ii-1; j<i; j++) {
                sum -= LU[i][j]*x[j];
            }
        }

        else if (sum != 0.0) {
            ii=i+1;
        }

        x[i] = sum;
    }

    for (int k = n-1; k >= 0; k--) { // Back substitution

        sum=x[k];
        for (j=k+1;j<n;j++) {
            sum -= LU[k][j]*x[j];
        }
        x[k]=sum/LU[k][k];
    }
}

//----------------------------------------------------------------------------
// Solve for multiple rhs B (n x m), and store the result in X (n x m)
//----------------------------------------------------------------------------

void LUdcmp::solve(const Matrix& B, int m, Matrix& X) const {

    Vector xx(extents[n]);
    int i, j;
    
    for (j = 0; j < m; j++) {
        for (i = 0; i < n; i++) {
            xx[i] = B[i][j];
        }

        solve(xx, xx);

        for (i = 0;i < n; i++) {
            X[i][j] = xx[i];
        }
    }
}

//----------------------------------------------------------------------------
// Find the inverse of a matrix A
//----------------------------------------------------------------------------

void LUdcmp::inverse(Matrix& A_inv) const {

    int i,j;

    for (i=0;i<n;i++) {
    	for (j=0;j<n;j++)
            A_inv[i][j] = 0.;

    	A_inv[i][i] = 1.;
    }

    solve(A_inv, n, A_inv);
}

//----------------------------------------------------------------------------
// Find the determinant of the matrix
//----------------------------------------------------------------------------

double LUdcmp::det() const {

    //WARNING: For large systems, the determinant of the matrix may overflow or underflow

    double dd = d;

    for (int i=0; i < n; i++)
        dd *= LU[i][i];

    return dd;
}

//----------------------------------------------------------------------------
// Default constructor
//----------------------------------------------------------------------------


QRdcmp::QRdcmp() :
    m(0),
	n(0)
{
    
    
}

//----------------------------------------------------------------------------
// Constructor for QRdcmp class. The input is a rectangular (m x n, m >= n)
// matrix
//----------------------------------------------------------------------------

QRdcmp::QRdcmp(Matrix& A, int m_, int n_) :
    QR_(A),
	m(m_),
	n(n_),
    Rdiag(extents[n_])
{

    double s, nrm;
	int i = 0, j = 0, k = 0;

    // Main loop.
	
    for (k = 0; k < n; k++) {

		// Compute 2-norm of k-th column without under/overflow.

    	nrm = 0;
    	for (i = k; i < m; i++) {
    		nrm = std::hypot(nrm,QR_[i][k]);
    	}

    	if (nrm != 0.0) {
          // Form k-th Householder vector.
          if (QR_[k][k] < 0) {
             nrm = -nrm;
          }
          for (i = k; i < m; i++) {
             QR_[i][k] /= nrm;
          }
          QR_[k][k] += 1.0;

          // Apply transformation to remaining columns.
          for (j = k+1; j < n; j++) {
             s = 0.0;
             for (i = k; i < m; i++) {
                s += QR_[i][k]*QR_[i][j];
             }
             s = -s/QR_[k][k];
             for (i = k; i < m; i++) {
                QR_[i][j] += s*QR_[i][k];
             }
          }
       }

    	Rdiag[k] = -nrm;
    }
}

//----------------------------------------------------------------------------
// Reinitialize the class to a different value
//----------------------------------------------------------------------------

void QRdcmp::initialize(Matrix& A, int m_, int n_)
{   
    
    m = m_; n = n_; 
    
    QR_.resize(extents[m][n]);
	Rdiag.resize(extents[n]);
    
    double s, nrm;
	int i = 0, j = 0, k = 0;

    for (i = 0; i < m; ++i)
        for (j = 0; j < n; ++j)
            QR_[i][j] = A[i][j]; 
    
    // Main loop.
	
    for (k = 0; k < n; k++) {

		// Compute 2-norm of k-th column without under/overflow.

    	nrm = 0;
    	for (i = k; i < m; i++) {
    		nrm = std::hypot(nrm,QR_[i][k]);
    	}

    	if (nrm != 0.0) {
          // Form k-th Householder vector.
          if (QR_[k][k] < 0) {
             nrm = -nrm;
          }
          for (i = k; i < m; i++) {
             QR_[i][k] /= nrm;
          }
          QR_[k][k] += 1.0;

          // Apply transformation to remaining columns.
          for (j = k+1; j < n; j++) {
             s = 0.0;
             for (i = k; i < m; i++) {
                s += QR_[i][k]*QR_[i][j];
             }
             s = -s/QR_[k][k];
             for (i = k; i < m; i++) {
                QR_[i][j] += s*QR_[i][k];
             }
          }
       }

    	Rdiag[k] = -nrm;
    }
    
}

//----------------------------------------------------------------------------
// Check wether the given system is a full rank system
//----------------------------------------------------------------------------

bool QRdcmp::is_full_rank() const {

	for (int j = 0; j < n; j++) {
       if (Rdiag[j] == 0.0)
          return false;
    }

    return true;
}

//----------------------------------------------------------------------------
// Get the R (n x n) matrix
//----------------------------------------------------------------------------

void QRdcmp::get_R(Matrix& R) const {

	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			if (i < j) {
				R[i][j] = QR_[i][j];
			} else if (i == j) {
				R[i][j] = Rdiag[i];
			} else {
				R[i][j] = 0.0;
			}
       }
    }
}

//----------------------------------------------------------------------------
// Get the Q (m x n) matrix
//----------------------------------------------------------------------------

void QRdcmp::get_Q(Matrix& Q) const {

	  int i = 0, j = 0, k = 0;
	  double s;

	  for (k = n-1; k >= 0; k--) {

		  for (i = 0; i < m; i++)
			  Q[i][k] = 0.0;

    	Q[k][k] = 1.0;

    	for (j = k; j < n; j++) {
    		if (QR_[k][k] != 0) {

    			s = 0.0;

    			for (i = k; i < m; i++)
    				s += QR_[i][k]*Q[i][j];

    			s = -s/QR_[k][k];

    			for (i = k; i < m; i++)
    				Q[i][j] += s*QR_[i][k];

    		}
    	}
    }
}

//----------------------------------------------------------------------------
// Solve least squares system A*x = b. b is m length vector and x is n length
// vector
//----------------------------------------------------------------------------

void QRdcmp::solve(const Vector& b, Vector& x_) const {

	Vector x(b);
	double s;

    // Compute Y = transpose(Q)*b

	for (int k = 0; k < n; k++) {
		s = 0.0;
		for (int i = k; i < m; i++)
			s += QR_[i][k]*x[i];

		s = -s/QR_[k][k];

		for (int i = k; i < m; i++)
			x[i] += s*QR_[i][k];
    }

    // Solve R*X = Y;

    for (int k = n-1; k >= 0; k--) {

    	x[k] /= Rdiag[k];

       for (int i = 0; i < k; i++)
             x[i] -= x[k]*QR_[i][k];

    }

    // return n x nx portion of X

    for (int i=0; i<n; i++)
    	x_[i] = x[i];
}

//----------------------------------------------------------------------------
// Solve least squares system A*X = B. B is (m x k) array and
// x is (n x k) array
//----------------------------------------------------------------------------

void QRdcmp::solve(const Matrix& B, int k_, Matrix& X_) const {

	Matrix X(B);

    int nx = k_;
    double s;

	int i = 0, j = 0, k = 0;

    // Compute Y = transpose(Q)*B

	for (k = 0; k < n; k++) {
		for (j = 0; j < nx; j++) {

			s = 0.0;

			for (i = k; i < m; i++)
				s += QR_[i][k]*X[i][j];


			s = -s/QR_[k][k];

			for (i = k; i < m; i++)
				X[i][j] += s*QR_[i][k];

		}
    }

	// Solve R*X = Y;

	for (k = n-1; k >= 0; k--) {

		for (j = 0; j < nx; j++)
			X[k][j] /= Rdiag[k];


		for (i = 0; i < k; i++)
			for (j = 0; j < nx; j++)
				X[i][j] -= X[k][j]*QR_[i][k];
	}


	//  return n x nx portion of X 

	  for (i=0; i<n; i++)
	  	for (j=0; j<nx; j++)
			X_[i][j] = X[i][j];
}