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
// For a given state Q, find all the eigenvalues L in the normal direction 
//----------------------------------------------------------------------------

void PDEEigenvalues(const Vector& Q, const double& nx, const double& ny, Vector& L) {
    
    double rho = Q[0];
    double u = Q[1]/Q[0]; 
    double v = Q[2]/Q[0]; 
    double p = (GAMMA - 1.0)*(Q[3] - 0.5*(Q[1]*Q[1] + Q[2]*Q[2])/Q[0] );
    double un = u*nx + v*ny;
    double a = std::sqrt(GAMMA*p/rho);

    L[0] = un - a;
    L[1] = un + a;
    L[2] = un;
    L[3] = un; 
}

//----------------------------------------------------------------------------
// For a given state Q, find the intermediate states 
//----------------------------------------------------------------------------

void PDEIntermediateFields(const Vector& Q, const double& nx, const double& ny, Matrix& Lambda, Matrix& RS, Matrix& LS) {
    
    int i, j; 
    Matrix DQ_DV(extents[nVar][nVar]); Matrix DV_DQ(extents[nVar][nVar]);      // Jacobian matrices
    Matrix RS_p(extents[nVar][nLin]), LS_p(extents[nLin][nVar]);               // Right and left eigenvectors in primitive variables ; 
    
    // First extract states 
    
    double rho = Q[0];
    double u = Q[1]/Q[0]; 
    double v = Q[2]/Q[0]; 
    double p = (GAMMA - 1.0)*(Q[3] - 0.5*(Q[1]*Q[1] + Q[2]*Q[2])/Q[0] );
    double un = u*nx + v*ny;
    double tempa = GAMMA-1.0; 
    
    // First fill the diagonal intermediate characteristic fields matrix Lambda
    
    for (i = 0; i < nLin; ++i) 
        for (j = 0; j < nLin; ++j)
            Lambda[i][j] = 0.0; 
    
    Lambda[0][0] = un; 
    Lambda[1][1] = un; 
    
    //-------------------------------------------------------------------------------------------------------
    
    // Fill the Jacobian matrices 
    
    for (i = 0; i < nVar; ++i) { 
        for (j = 0; j < nVar; ++j) {
            DQ_DV[i][j] = 0.0; DV_DQ[i][j] = 0.0; 
        }
    }
    
    DQ_DV[0][0] = 1.0; // Row 0
    
    DQ_DV[1][0] = u; 
    DQ_DV[1][1] = rho; // Row 1
    
    DQ_DV[2][0] = v; 
    DQ_DV[2][2] = rho; // Row 2
    
    DQ_DV[3][0] = 0.5*(u*u + v*v); 
    DQ_DV[3][1] = rho*u;
    DQ_DV[3][2] = rho*v;
    DQ_DV[3][3] = 1.0/tempa; // Row 3

    //-------------------------------------------------------------------------------------------------------

    DV_DQ[0][0] = 1.0; // Row 0
    
    DV_DQ[1][0] = -u/rho; 
    DV_DQ[1][1] = 1.0/rho; // Row 1
    
    DV_DQ[2][0] = -v/rho; 
    DV_DQ[2][2] = 1.0/rho; // Row 2
    
    DV_DQ[3][0] = 0.5*tempa*(u*u + v*v); 
    DV_DQ[3][1] = -u*tempa;
    DV_DQ[3][2] = -v*tempa;
    DV_DQ[3][3] = tempa; // Row 3

    //-------------------------------------------------------------------------------------------------------
    
    // Form the right eigenvectors in primitive variables 
    
    
    for (i = 0; i < nVar; ++i) 
        for (j = 0; j < nLin; ++j) 
            RS_p[i][j] = 0.0;  
    
    RS_p[0][0] = 1.0; 
    
    RS_p[1][1] = -ny; 
    RS_p[2][1] =  nx; 

    //-------------------------------------------------------------------------------------------------------
    
    // Form the left eigenvectors in primitive variables
    
    for (i = 0; i < nLin; ++i) 
        for (j = 0; j < nVar; ++j) 
            LS_p[i][j] = 0.0;
    
    LS_p[0][0] = 1.0; 
    LS_p[0][3] = -rho/(GAMMA*p);
    
    LS_p[1][1] = -ny; 
    LS_p[1][2] =  nx; 

    //-------------------------------------------------------------------------------------------------------
    
    // Now obtain the eigenvector matrices in conserved variables 
    
    MatMatMult(LS_p, DV_DQ, LS, nLin, nVar, nVar);
    MatMatMult(DQ_DV, RS_p, RS, nVar, nVar, nLin);
}

//----------------------------------------------------------------------------
// HLLEM Flux in normal direction F.n 
//----------------------------------------------------------------------------

double PDEHLLEMFlux(const Vector& QL, const Vector& QR, 
                    const double& nx, const double& ny, 
                    const double& x, const double& y, 
                    Vector& F) {
        
    // Declare variables 
    
    int i, j;
    Vector QM(extents[nVar]), FL(extents[nVar]), FR(extents[nVar]), FM(extents[nVar]); 
    Vector LL(extents[nVar]), LR(extents[nVar]), LM(extents[nVar]); 
    Vector QHLL(extents[nVar]);
    Vector AD(extents[nVar]), Q_jump(extents[nVar]), Work1(extents[nVar]), Work2(extents[nVar]), Work3(extents[nLin]), Work4(extents[nLin]); 
    Matrix RS(extents[nVar][nLin]), LS(extents[nLin][nVar]), Lambda(extents[nLin][nLin]), Id(extents[nLin][nLin]);
    Matrix Delta(extents[nLin][nLin]), Lap(extents[nLin][nLin]), Lam(extents[nLin][nLin]);
    double sL, sR, smaxL, smaxR;
    const double flattener = 1.0; 
    
    // Compute the intermediate state
    
    for (i = 0; i < nVar; ++i)
        QM[i] = 0.5*(QL[i] + QR[i]);
    
    // Compute the fluxes FL and FR 
    
    smaxL = PDEConsFlux(QL, nx, ny, x, y, FL);
    smaxR = PDEConsFlux(QR, nx, ny, x, y, FR);
    
    // Compute the eigenvalues of QM, QL and QR 

    PDEEigenvalues(QM, nx, ny, LM);
    PDEEigenvalues(QL, nx, ny, LL);
    PDEEigenvalues(QR, nx, ny, LR);
    
    // Compute the left and right  signal speed 
    
    sL = std::min(0.0, std::min(MinVal(LM, nVar), MinVal(LL, nVar)) ); 
    sR = std::max(0.0, std::max(MaxVal(LM, nVar), MaxVal(LR, nVar)) );
    
    // Find the HLL state  
        
    for (i = 0; i < nVar; ++i)
        QHLL[i] = ( QR[i]*sR - QL[i]*sL - (FR[i] - FL[i]) ) / (sR-sL); 
    
    // Find the HLL Flux 

    for (i = 0; i < nVar; ++i)
        F[i] = (sR*FL[i] - sL*FR[i])/(sR - sL) + ((sL*sR)/(sR - sL))*(QR[i] - QL[i]);
    
    //-------------------------------------------------------------------------------------------------------
    // Add the anti-diffusive term of HLLEM-RS
    
    // Form identity matrix 
    
    for (i = 0; i < nLin; ++i) {
        for (j = 0; j < nLin; ++j) {
            if (i == j)
                Id[i][j] = 1.0; 
            else
                Id[i][j] = 0.0;
        }
    }
     
    // Compute intermediate eigenvalues and eigenvectors 

    PDEIntermediateFields(QM, nx, ny, Lambda, RS, LS);
    
    // Compute Lambda+ and Lambda- 
    
    for (i = 0; i < nLin; ++i) {
        for (j = 0; j < nLin; ++j) {
            Lap[i][j] = 0.5*(Lambda[i][j] + std::abs(Lambda[i][j]));
            Lam[i][j] = 0.5*(Lambda[i][j] - std::abs(Lambda[i][j]));
        }
    }
    
    // Compute Delta*
    
    for (i = 0; i < nLin; ++i) {
        for (j = 0; j < nLin; ++j) {
            Delta[i][j] = Id[i][j] - Lam[i][j]/(sL - 1.0e-14) - Lap[i][j]/(sR + 1.0e-14);
        }
    }
    
    for (i = 0; i < nVar; ++i)
        Q_jump[i] = QR[i] - QL[i]; 
    
    MatVecMult(LS, Q_jump, Work3, nLin, nVar);
    MatVecMult(Delta, Work3, Work4, nLin, nLin);
    MatVecMult(RS, Work4, AD, nVar, nLin);
    
    for (i = 0; i < nVar; ++i)
        AD[i] = (sR*sL)/(sR - sL)*AD[i]; 

    for (i = 0; i < nVar; ++i)
        F[i] = F[i] - flattener*AD[i];
    
    return std::max(smaxL, smaxR); 
}

//----------------------------------------------------------------------------
// Rusanov (LLF) flux in normal direction   
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