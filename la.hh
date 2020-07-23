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


#endif /* LA_HH_ */