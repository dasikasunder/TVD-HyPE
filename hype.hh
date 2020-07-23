/*
 * hype.hh
 *      Author: sunder
 */ 

#ifndef HYPE_HH_
#define HYPE_HH_

#include "pde.hh"
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

class HyPE_2D {

	// Typedefs

    typedef boost::multi_array<Vector, 1> array_type;
    typedef array_type::index index;
	
	AppCtx Params;          /* Parameters controlling the solution */
	Triangulation tria;     /* Object representing the mesh */
	
	Vector (*init_func)(double, double);  /* Initial condition function */
	
	double time; 
	double dt;
	double h_min; 
	
	multi_array<double, 3> U;   /* Conservative variables along with its components in each cell */
	multi_array<double, 2> W;   /* Primitive variables in each cell */
	multi_array<double, 2> F;   /* Upwind flux on each in each face */
	multi_array<double, 2> RHS; /* Update coefficient for each cell */
	
	void initialize();     
	void assign_boundary_ids();
	void compute_rhs();
	void compute_primitive_variables(); 
	void solve();
	void plot_vtk(int = 0, const int = 4);
	
public:
	
	HyPE_2D(Vector(*)(double, double), AppCtx);
	void run();
};



#endif /* HYPE_H_ */
