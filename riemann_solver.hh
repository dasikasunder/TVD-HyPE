/*
 * riemann_solver.hh
 *      Author: sunder
 */

#ifndef RIEMANN_SOLVER_HH_
#define RIEMANN_SOLVER_HH_

#include "headers.hh"

class Riemann_Solver {

	/* Properties of the gas */
	
	double gamma; 
	
public:
	
	/* Number of components in the system */
	
	static const int n_comp = 4;
	
	/* Constructors */
	
	Riemann_Solver();
	Riemann_Solver(const double&); 
	
	/* Convert between conserved and primitive variables */
	
	void conserved_to_primitive(const std::vector<double>&, std::vector<double>&);
	void primitive_to_conserved(const std::vector<double>&, std::vector<double>&);

	/* Find conservative and non-conservative fluxes */
	
	double conservative_flux(const std::vector<double>&, const double&, const double&, const double&,
							 const double&, std::vector<double>&);
	
	void non_conservative_flux(const std::vector<double>&, const std::vector<double>&, const std::vector<double>&,
							   std::vector<double>&);
	
	void matrix_B(const std::vector<double>&, const double&, const double&, double**); 

	
	void roe_matrix(const std::vector<double>&, const std::vector<double>&,
							const double&, const double&,
							double**);
	
	double llf_riemann_solver(const std::vector<double>&, 
								    const std::vector<double>&,
									const double&, 
									const double&,
									const double&,
									const double&,
								    std::vector<double>&); 
	 
	

};





#endif /* CLAW_HH_ */
