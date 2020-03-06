/*
 * hype.hh
 *      Author: sunder
 */ 

#ifndef HYPE_HH_
#define HYPE_HH_

#include "riemann_solver.hh"
#include "headers.hh"
#include "triangulation.hh"

struct AppCtx {
	std::string mesh_file_name; 
	double CFL;
	double initial_time = 0.0;
	double final_time;
	bool apply_limiter = true;
	int write_interval = 5;
	double gamma = 1.4; 
};

class Hype {
	
	AppCtx Params;          /* Parameters controlling the solution */
	Triangulation tria;     /* Object representing the mesh */
	
	std::vector<double> (*init_func)(double, double);  /* Initial condition function */
	Riemann_Solver riemann;
	
	double time; 
	double dt;
	double h_min; 
	
	double*** U;            /* Store the conservative variables along with its components in each cell */
	double**  W;            /* Store the primitive variables in each cell */
	double**  F;            /* Store the flux on each in each face */
	double**  RHS;          /* Store the RHS value for each cell */
	
	void initialize();     
	void assign_boundary_ids();
	void compute_rhs();
	void compute_primitive_variables(); 
	void solve();
	void plot_vtk(int = 0, const int = 4);
	
public:
	
	Hype(std::vector<double>(*)(double, double), AppCtx);
	~Hype(); 
	void run();
};



#endif /* HYPE_H_ */
