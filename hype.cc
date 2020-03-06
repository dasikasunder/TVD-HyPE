#include "hype.hh"

//----------------------------------------------------------------------------
// Constructor - Takes input parameters and sets memory 
//----------------------------------------------------------------------------

Hype::Hype(std::vector<double> (*f)(double, double), AppCtx params):
	Params(params),
	tria(params.mesh_file_name),
	init_func(f),
	riemann(params.gamma),
	time(params.initial_time),
	dt(0.0)
{
	
	std::cout << "----------------------------------------------------------------------------" << std::endl;
	std::cout << "Allocating memory " << std::endl;
	std::cout << "No. of cells in the triangulation = " << tria.no_cells() << std::endl; 
	
	// Allocate memory for all the variables 
	
	Utilities::allocate_mem_3d(U,   tria.no_cells(), Riemann_Solver::n_comp, 4);
	Utilities::allocate_mem_2d(W,   tria.no_cells(), Riemann_Solver::n_comp);
	Utilities::allocate_mem_2d(RHS, tria.no_cells(), Riemann_Solver::n_comp);
	Utilities::allocate_mem_2d(F,   tria.no_faces(), Riemann_Solver::n_comp);
	
	h_min = tria.min_incircle_dia();

	std::cout << "Done !" << std::endl;
	std::cout << "----------------------------------------------------------------------------" << std::endl;
}

//----------------------------------------------------------------------------
// Destructor - Free all the memory  
//----------------------------------------------------------------------------

Hype::~Hype() {
	
	// Free all the memory 
	
	Utilities::free_mem_3d(U,   tria.no_cells(), Riemann_Solver::n_comp);
	Utilities::free_mem_2d(W,   tria.no_cells());
	Utilities::free_mem_2d(RHS, tria.no_cells());
	Utilities::free_mem_2d(F,   tria.no_faces());
	
}

//----------------------------------------------------------------------------
// Initialize the solution using the initial conditions 
//----------------------------------------------------------------------------

void Hype::initialize() {
	
	std::cout << "Initializing the solution" << std::endl; 

	std::vector<double> V(Riemann_Solver::n_comp);
	std::vector<double> Q(Riemann_Solver::n_comp); 
	
	
	for (int c = 0; c < tria.no_cells(); ++c) {
		
		V = init_func(tria.xc(c), tria.yc(c)); 
		
		riemann.primitive_to_conserved(V, Q); 
		
		for (int n_c = 0; n_c < Riemann_Solver::n_comp; ++n_c)
			U[c][n_c][0] = Q[n_c]; 
			
	}
	
	std::cout << "Done !" << std::endl;
	std::cout << "----------------------------------------------------------------------------" << std::endl; 
}

//----------------------------------------------------------------------------
// To each face in the mesh assign boundary ids 
//----------------------------------------------------------------------------

void Hype::assign_boundary_ids() {
	
	std::cout << "Assigning boundary id's" << std::endl; 
	
	// Usual convention: 0 -> Transmissive/Outlet; 1 -> Inlet/Dirichlet; 2 -> Reflective 

	// Loop over all the faces 
	
	for (int f = 0; f < tria.no_faces(); ++f) {
		
		if (tria.face_at_boundary(f)) {
			
			if (std::abs(tria.xf(f)) < 1.0e-12)              // Inlet of the boundary 
				tria.set_boundary_id(f, 1);
				
			else if (std::abs(tria.xf(f) - 3.0) < 1.0e-12)   // Outflow  
				tria.set_boundary_id(f, 0);
			
			else                                             // Reflective 
				tria.set_boundary_id(f, 2);
		}
		
	}
	
	std::cout << "Done !" << std::endl;
	std::cout << "----------------------------------------------------------------------------" << std::endl; 
}

//----------------------------------------------------------------------------
// In each cell, find the primitive variables 
//----------------------------------------------------------------------------

void Hype::compute_primitive_variables() {
	
	std::vector<double> V(Riemann_Solver::n_comp);
	std::vector<double> Q(Riemann_Solver::n_comp);  
	
	for (int c = 0; c < tria.no_cells(); ++c) {
		
		for (int n_c = 0; n_c < Riemann_Solver::n_comp; ++n_c)
			Q[n_c] = U[c][n_c][0]; 
		
		riemann.conserved_to_primitive(Q, V);
		
		for (int n_c = 0; n_c < Riemann_Solver::n_comp; ++n_c)
			W[c][n_c] = V[c]; 
	}
}

//----------------------------------------------------------------------------
// Find the RHS in each cell and also compute the time step size 
//----------------------------------------------------------------------------

void Hype::compute_rhs() {
	
	std::vector<double> QL(Riemann_Solver::n_comp); 
	std::vector<double> QR(Riemann_Solver::n_comp); 
	std::vector<double> Flux(Riemann_Solver::n_comp);
	
	std::vector<double> VR(Riemann_Solver::n_comp);
	
	double s, s_max = 0.0; 
	
	// Find the upwind flux on each face 
	
	int L_cell_index, R_cell_index;
	
	for (int f = 0; f < tria.no_faces(); ++f) {
		
		L_cell_index = tria.link_face_to_cell(f, 0);
		
		for (int n_c = 0; n_c < Riemann_Solver::n_comp; ++n_c)
			QL[n_c] = U[L_cell_index][n_c][0];
		
		// Apply boundary conditions for faces at boundaries  
		
		if (tria.face_at_boundary(f)) {
			
			if (tria.get_boundary_id(f) == 1) { // Inlet 
				
				VR[0] = 1.4; // Density  
				VR[1] = 3.0; // x-Velocity 
				VR[2] = 0.0; // y-Velocity 
				VR[3] = 1.0; // Pressure 
			
				riemann.primitive_to_conserved(VR, QR);
			}
			
			else if (tria.get_boundary_id(f) == 2) { // Reflecting wall 
				
				QR[0] = QL[0];
				QR[1] = QL[1] - 2.0*QL[1]*tria.nx(f)*tria.nx(f) - 2.0*QL[2]*tria.nx(f)*tria.ny(f); 
				QR[2] = QL[2] - 2.0*QL[1]*tria.nx(f)*tria.ny(f) - 2.0*QL[2]*tria.ny(f)*tria.ny(f);
				QR[3] = QL[3];
			}
			
			else { // Transmissive
			
				QR[0] = QL[0]; 
				QR[1] = QL[1];
				QR[2] = QL[2]; 
				QR[3] = QL[3];
			}
		}
		
		else {
			
			R_cell_index = tria.link_face_to_cell(f, 1);
			
			for (int n_c = 0; n_c < Riemann_Solver::n_comp; ++n_c)
				QR[n_c] = U[R_cell_index][n_c][0];
		}
		
		s = riemann.llf_riemann_solver(QL, QR, tria.nx(f), tria.ny(f), tria.xf(f), tria.yf(f), Flux);
		
		if (s > s_max)
			s_max = s; 
		
		for (int n_c = 0; n_c < Riemann_Solver::n_comp; ++n_c)
			F[f][n_c] = Flux[n_c];
	}
	
	// Now find the RHS in each cell 
	
	int global_f; double r1_v; 
	
	for (int c = 0; c < tria.no_cells(); ++c) {
		
		for (int n_c = 0; n_c < Riemann_Solver::n_comp; ++n_c)
			RHS[c][n_c] = 0.0; 
			
		r1_v = 1./tria.vol(c);
		
		for (int local_f = 0; local_f < GeometryInfo::faces_per_cell; ++local_f) {
			
			global_f = tria.link_cell_to_face(c, local_f);
			
			for (int n_c = 0; n_c < Riemann_Solver::n_comp; ++n_c)
				RHS[c][n_c] += -r1_v*tria.n_sign(c, local_f)*F[global_f][n_c]*tria.areaf(global_f);
		}
	}
	
	// Compute the time step value 
	
	dt = (Params.CFL*h_min)/s_max;
	
	//  Check size of dt to avoid exceeding output time

	if((time + dt)>Params.final_time)
		dt = Params.final_time - time;
}

//----------------------------------------------------------------------------
// Solve the system 
//----------------------------------------------------------------------------

void Hype::solve() {

	while (time < Params.final_time) {
		
		compute_rhs(); 
		
		printf ("time = %4.3e, dt = %4.3e, final time = %4.3e\n", time, dt, Params.final_time);
		
		for (int c = 0; c < tria.no_cells(); ++c) {
			
			for (int n_c = 0; n_c < Riemann_Solver::n_comp; ++n_c)
				U[c][n_c][0] += dt*RHS[c][n_c]; 
		}
		
		time += dt; 
	}
	
	printf ("time = %4.3e, dt = %4.3e, final time = %4.3e\n", time, dt, Params.final_time);
	
	std::cout << "Done !" << std::endl;
	std::cout << "----------------------------------------------------------------------------" << std::endl;
}

//----------------------------------------------------------------------------
// Plot solution in vtk format 
//----------------------------------------------------------------------------

void Hype::plot_vtk(int i, const int digits) {

	std::ofstream vtk;
	const std::string filename = "../plots/plot_" + Utilities::int_to_string (i, digits) + ".vtk";
	vtk.open (filename);
	vtk.flags( std::ios::dec | std::ios::scientific );
	vtk.precision(6);

	Assert(vtk.is_open(), ErrIO(filename));
	
	vtk << "# vtk DataFile Version 3.0" << "\n";
	vtk << "2D Euler" << "\n";
	vtk << "ASCII" << "\n";
	vtk << "\nDATASET UNSTRUCTURED_GRID" << "\n";
	vtk << "\nFIELD FieldData 1" << "\n";
	vtk << "TIME 1 1 double" << "\n";
	vtk << time << "\n";
	
	// Point data 
	
	vtk << "POINTS " << tria.no_vertices() << " double" << "\n"; 
	
	for (int v = 0; v < tria.no_vertices() ; ++v)
		vtk << tria.xv(v) << " " << tria.yv(v) << " " << 0.0 << "\n"; 
	
	// Cell connectivity 
	
	vtk << "CELLS " << tria.no_cells() << " " << tria.no_cells()*4  << "\n";
	
	for (int c = 0; c < tria.no_cells(); ++c) {
		
		vtk << 3 << " ";
		
		for (int local_v = 0; local_v < GeometryInfo::vertices_per_cell; ++local_v )
			vtk << tria.link_cell_to_vertex(c, local_v) << " ";

		vtk << "\n"; 
	}
	
	// Cell type - 5 for triangles 
	
	vtk << "CELL_TYPES " << tria.no_cells() << "\n";
	
	for (int c = 0; c < tria.no_cells(); ++c) 
		vtk << 5 << "\n";
	
	// --------------------------------------------------------------------------------------
	// Add scalars or vectors here. No need to change anything above this 
	//--------------------------------------------------------------------------------------
	
	compute_primitive_variables();
	
	double* RHO_vtrx = new double[tria.no_vertices()];
	double* Vx_vtrx  = new double[tria.no_vertices()];
	double* Vy_vtrx  = new double[tria.no_vertices()];
	double* P_vtrx   = new double[tria.no_vertices()];
	
	for (int v = 0; v < tria.no_vertices(); ++v) {
	
		RHO_vtrx[v] = 0.0; Vx_vtrx[v] = 0.0; 
		Vy_vtrx[v]  = 0.0; P_vtrx[v]  = 0.0; 
		
		for (int i = 0; i < tria.no_cells_sharing_vertex(v); ++i) {
			RHO_vtrx[v] += W[tria.cell_sharing_vertex(v,i)][0];
			Vx_vtrx[v]  += W[tria.cell_sharing_vertex(v,i)][1];
			Vy_vtrx[v]  += W[tria.cell_sharing_vertex(v,i)][2];
			P_vtrx[v]   += W[tria.cell_sharing_vertex(v,i)][3];
		}
		
		RHO_vtrx[v] /= static_cast<double>(tria.no_cells_sharing_vertex(v));
		Vx_vtrx[v]  /= static_cast<double>(tria.no_cells_sharing_vertex(v));
		Vy_vtrx[v]  /= static_cast<double>(tria.no_cells_sharing_vertex(v));
		P_vtrx[v]   /= static_cast<double>(tria.no_cells_sharing_vertex(v));
	}
	
	vtk << "POINT_DATA " << tria.no_vertices() << "\n";
	vtk << "SCALARS Density double" << "\n"; 
	vtk << "LOOKUP_TABLE default" << "\n"; 
	
	for (int v = 0; v < tria.no_vertices(); ++v)
		vtk << RHO_vtrx[v] << "\n";
	
	vtk << "\n"; 
	
	vtk << "SCALARS Pressure double" << "\n"; 
	vtk << "LOOKUP_TABLE default" << "\n"; 
	for (int v = 0; v < tria.no_vertices(); ++v)
		vtk << P_vtrx[v] << "\n"; 

	vtk << "\n"; 
	
	vtk << "VECTORS Velocity double" << "\n"; 
	for (int v = 0; v < tria.no_vertices(); ++v)
		vtk << Vx_vtrx[v] << " " << Vy_vtrx[v] << " " << 0.0 << "\n";

	vtk.close(); 
	
	delete[] RHO_vtrx;
	delete[] Vx_vtrx;
	delete[] Vy_vtrx;
	delete[] P_vtrx;
}

//----------------------------------------------------------------------------
// Put everything together and run the problem 
//----------------------------------------------------------------------------

void Hype::run() {

	initialize(); 
	assign_boundary_ids();
	solve();
	plot_vtk();
	
}
