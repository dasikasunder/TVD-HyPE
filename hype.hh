/*
 * hype.hh
 *      Author: sunder
 */ 

#ifndef HYPE_HH_
#define HYPE_HH_

#include "la.hh"
#include "pde.hh"
#include "headers.hh"
#include "triangulation.hh"

struct AppCtx {
	std::string mesh_file_name; 
	double CFL;
	double tend;
	int write_interval = 100;
};

class HyPE_2D {

	// Typedefs

    typedef boost::multi_array<Vector, 1> array_type;
    typedef array_type::index index;
	
	AppCtx Params;          // Parameters controlling the solution 
	Triangulation tria;     // Object representing the mesh 
	
	Vector (*init_func)(double, double);  // Initial condition function 
	
	double time;                  // Current time of the simulation 
	double dt;                    // Time step 
	int time_step;                // Step number 
	double h_min;                 // Size of the smallest element in mesh 
	
	multi_array<double, 3> uh;    // Conservative variables along with its components in each cell */
	multi_array<double, 2> wh;    // Primitive variables in each cell 
	multi_array<double, 2> Fu;    // Upwind flux on each in each face in normal direction 
	multi_array<double, 2> duh;   // Update coefficient for each cell 
	
	multi_array<QRdcmp, 1> qr;    // QR decomposition class for reconstruction

	// 

	void initialize();  
	void assign_boundary_ids();
	void reconstruct();
	void compute_rhs();
	void compute_primitive_variables(); 
	void plot_vtk(int = 0, const int = 4);
	
public:
	
	HyPE_2D(Vector(*)(double, double), AppCtx);
	void run();
};



#endif /* HYPE_H_ */
