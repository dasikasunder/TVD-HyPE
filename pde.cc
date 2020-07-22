/*
 * pde.cc
 *      Author: sunder
 */

// Definitions of PDE related functions 

#include "pde.hh"

//----------------------------------------------------------------------------
// Convert conserved variables to primitive variables 
//----------------------------------------------------------------------------

void PDECons2Prim(const Vector& Q, Vector& V) {

	V[0] = Q[0];
	V[1] = Q[1]/Q[0];
	V[2] = Q[2]/Q[0];
	V[3] = (GAMMA -1.0)*( Q[3] - 0.5*(Q[1]*Q[1] + Q[2]*Q[2])/Q[0] );
}

//----------------------------------------------------------------------------
// Convert primitive variables to conserved variables 
//----------------------------------------------------------------------------

void PDEPrim2Cons(const Vector& V, Vector& Q) {
	
	double e = (V[3])/(GAMMA - 1.0);
	double k = 0.5*V[0]*(V[1]*V[1] + V[2]*V[2]);

	Q[0] = V[0];
	Q[1] = V[0]*V[1];
	Q[2] = V[0]*V[2];
	Q[3] = k + e;
}

//----------------------------------------------------------------------------
// Find the conservative flux in the normal direction F.n 
//----------------------------------------------------------------------------

double PDEConsFlux(const Vector& Q, const double& nx, const double& ny,
									const double& x, const double& y,
									Vector& F) {

	double p = (GAMMA -1.0)*( Q[3] - 0.5*(Q[1]*Q[1] + Q[2]*Q[2])/Q[0] );

	if (p < 0.0) {
		std::cerr << "Neagative pressure, p = " << p << std::endl;
		std::cerr << "At x = " << x << ", y = " << y << std::endl;  
		std::exit(EXIT_FAILURE);
	}

	if (Q[0] < 0.0) {
		std::cerr << "Neagative density, rho = " <<  Q[0] << std::endl;
		std::cerr << "At x = " << x << ", y = " << y << std::endl; 
		std::exit(EXIT_FAILURE);
	}

	F[0] = nx*Q[1] + ny*Q[2];
	F[1] = nx*(Q[1]*Q[1]/Q[0] + p) + ny*(Q[1]*Q[2]/Q[0]);
	F[2] = nx*(Q[1]*Q[2]/Q[0]) + ny*(Q[2]*Q[2]/Q[0] + p);
	F[3] = nx*(Q[1]*(Q[3] + p)/Q[0]) + ny*(Q[2]*(Q[3] + p)/Q[0]);
	
	// Also obtain the maximum eigen value 
	
	double s_max = std::abs(Q[1]*nx + Q[2]*ny)/Q[0] + std::sqrt(GAMMA*p/Q[0]);
	
	return s_max;
}

//----------------------------------------------------------------------------
// Find the Rusanov (LLF) flux  
//----------------------------------------------------------------------------

double PDERusanovFlux(const Vector& QL, const Vector& QR,
					  const double& nx, const double& ny,
					  const double& x, const double& y,
					  Vector& Flux) {

    Vector FL(extents[nVar]); 
    Vector FR(extents[nVar]); 

    double smaxL = PDEConsFlux(QL, nx, ny, x, y, FL); 
    double smaxR = PDEConsFlux(QR, nx, ny, x, y, FR);

    double smax = std::max(smaxL, smaxR);

    for (int c = 0; c < nVar; ++c)
        Flux[c] = 0.5*(FR[c] + FL[c] - smax*(QR[c] - QL[c]));

    return smax; 
}