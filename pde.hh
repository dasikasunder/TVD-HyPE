/*
 * pde.hh
 *      Author: sunder
 */

#ifndef PDE_HH_
#define PDE_HH_

# include "headers.hh"
#include "la.hh"

const int nVar = 4;       // Number of variables in the PDE system 
const int nLin = 2;       // Number of linearly intermediate fields 
const double GAMMA = 1.4; // Ratio of specific heats of the gas 

void PDECons2Prim(const Vector&, Vector&);
void PDEPrim2Cons(const Vector&, Vector&);
double PDEConsFlux(const Vector&, const double&, const double&, const double&, const double&, Vector&);
void PDEEigenvalues(const Vector&, const double&, const double&, Vector&);
void PDEIntermediateFields(const Vector&, const double&, const double&, Matrix&, Matrix&, Matrix&);
double PDERusanovFlux(const Vector&, const Vector&,
					  const double&, const double&, const double&, const double&,
					  Vector&);
double PDEHLLEMFlux(const Vector&, const Vector&,
					  const double&, const double&, const double&, const double&,
					  Vector&);

#endif /* PDE_HH_ */