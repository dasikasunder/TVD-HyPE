#include "hype.hh"

std::vector<double> init_cond(double, double);


int main() {
	
	try {
		
		AppCtx Params; 
		
		Params.mesh_file_name = "../mesh.su2"; 
		Params.CFL = 0.4;
		Params.final_time = 4.0;
		Params.write_interval = 20;
		Params.gamma = 1.4;
		Hype Problem(init_cond, Params); 
		
		Problem.run();
	}
	
	catch (ErrorBase& E) {
		E.print_log_data(std::cout);
		E.print_debug();
	}
}

std::vector<double> init_cond(double x, double y) {
	
	x = x - 0.0;
	y = y - 0.0;

	std::vector<double> W(4); 
	
	W[0] = 1.4;
    W[1] = 3.0;
    W[2] = 0.0; 
    W[3] = 1.0;
	
	return W; 
}
