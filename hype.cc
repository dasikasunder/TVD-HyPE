#include "hype.hh"

//----------------------------------------------------------------------------
// Constructor - Takes input parameters and sets memory 
//----------------------------------------------------------------------------

HyPE_2D::HyPE_2D(Vector (*f)(double, double), AppCtx params):
	Params(params),
	tria(params.mesh_file_name),
	init_func(f),
	time(params.initial_time),
	dt(0.0)
{
	
	std::cout << "----------------------------------------------------------------------------" << std::endl;
	std::cout << "Allocating memory " << std::endl;
	std::cout << "No. of cells in the triangulation = " << tria.no_cells() << std::endl; 
	
	// Allocate memory for all the variables 
	
	U.resize(extents[tria.no_cells()][nVar][4]); 
	W.resize(extents[tria.no_cells()][nVar]);
	RHS.resize(extents[tria.no_cells()][nVar]);
	F.resize(extents[tria.no_faces()][nVar]);
	
	h_min = tria.min_incircle_dia();

	std::cout << "Done !" << std::endl;
	std::cout << "----------------------------------------------------------------------------" << std::endl;
}

//----------------------------------------------------------------------------
// Initialize the solution using the initial conditions 
//----------------------------------------------------------------------------

void HyPE_2D::initialize() {
	
	std::cout << "Initializing the solution" << std::endl; 

	Vector V(extents[nVar]), Q(extents[nVar]);
	
	for (int iCell = 0; iCell < tria.no_cells(); ++iCell) {
		
		V = init_func(tria.xc(iCell), tria.yc(iCell)); 
		
		PDEPrim2Cons(V, Q); 
		
		for (int iVar = 0; iVar < nVar; ++iVar)
			U[iCell][iVar][0] = Q[iVar]; 
			
	}
	
	std::cout << "Done !" << std::endl;
	std::cout << "----------------------------------------------------------------------------" << std::endl; 
}

//----------------------------------------------------------------------------
// To each face in the mesh assign boundary ids 
//----------------------------------------------------------------------------

void HyPE_2D::assign_boundary_ids() {
	
	std::cout << "Assigning boundary id's" << std::endl; 
	
	// Usual convention: 0 -> Transmissive/Outlet; 1 -> Inlet/Dirichlet; 2 -> Reflective 

	// Loop over all the faces 
	
	for (int iFace = 0; iFace < tria.no_faces(); ++iFace) {
		
		if (tria.face_at_boundary(iFace)) {
			
			if (std::abs(tria.xf(iFace)) < 1.0e-12)              // Inlet of the boundary 
				tria.set_boundary_id(iFace, 1);
				
			else if (std::abs(tria.xf(iFace) - 3.0) < 1.0e-12)   // Outflow  
				tria.set_boundary_id(iFace, 0);
			
			else                                             // Reflective 
				tria.set_boundary_id(iFace, 2);
		}
	}
	
	std::cout << "Done !" << std::endl;
	std::cout << "----------------------------------------------------------------------------" << std::endl; 
}

//----------------------------------------------------------------------------
// In each cell, find the primitive variables 
//----------------------------------------------------------------------------

void HyPE_2D::compute_primitive_variables() {
	
	Vector V(extents[nVar]), Q(extents[nVar]);
		
	for (int iCell = 0; iCell < tria.no_cells(); ++iCell) {
		
		for (int iVar = 0; iVar < nVar; ++iVar)
			Q[iVar] = U[iCell][iVar][0]; 
		
		PDECons2Prim(Q, V);
		
		for (int iVar = 0; iVar < nVar; ++iVar)
			W[iCell][iVar] = V[iVar]; 
	}
}

//----------------------------------------------------------------------------
// Find the RHS in each cell and also compute the time step size 
//----------------------------------------------------------------------------

void HyPE_2D::compute_rhs() {
	
	Vector QL(extents[nVar]), QR(extents[nVar]), VR(extents[nVar]), Flux(extents[nVar]); 

	double s, s_max = 0.0; 
	
	// Find the upwind flux on each face 
	
	int L_cell_index, R_cell_index;
	
	for (int iFace = 0; iFace < tria.no_faces(); ++iFace) {
		
		L_cell_index = tria.link_face_to_cell(iFace, 0);
		
		for (int iVar = 0; iVar < nVar; ++iVar)
			QL[iVar] = U[L_cell_index][iVar][0];
		
		// Apply boundary conditions for faces at boundaries  
		
		if (tria.face_at_boundary(iFace)) {
			
			if (tria.get_boundary_id(iFace) == 1) { // Inlet 
				
				VR[0] = 1.4; // Density  
				VR[1] = 3.0; // x-Velocity 
				VR[2] = 0.0; // y-Velocity 
				VR[3] = 1.0; // Pressure 
			
				PDEPrim2Cons(VR, QR);
			}
			
			else if (tria.get_boundary_id(iFace) == 2) { // Reflecting wall 
				
				QR[0] = QL[0];
				QR[1] = QL[1] - 2.0*QL[1]*tria.nx(iFace)*tria.nx(iFace) - 2.0*QL[2]*tria.nx(iFace)*tria.ny(iFace); 
				QR[2] = QL[2] - 2.0*QL[1]*tria.nx(iFace)*tria.ny(iFace) - 2.0*QL[2]*tria.ny(iFace)*tria.ny(iFace);
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
			
			R_cell_index = tria.link_face_to_cell(iFace, 1);
			
			for (int iVar = 0; iVar < nVar; ++iVar)
				QR[iVar] = U[R_cell_index][iVar][0];
		}
		
		s = PDERusanovFlux(QL, QR, tria.nx(iFace), tria.ny(iFace), tria.xf(iFace), tria.yf(iFace), Flux);
		
		if (s > s_max)
			s_max = s; 
		
		for (int iVar = 0; iVar < nVar; ++iVar)
			F[iFace][iVar] = Flux[iVar];
	}
	
	// Now find the RHS in each cell 
	
	int global_f; double r1_v; 
	
	for (int iCell = 0; iCell < tria.no_cells(); ++iCell) {
		
		for (int iVar = 0; iVar < nVar; ++iVar)
			RHS[iCell][iVar] = 0.0; 
			
		r1_v = 1./tria.vol(iCell);
		
		for (int local_f = 0; local_f < GeometryInfo::faces_per_cell; ++local_f) {
			
			global_f = tria.link_cell_to_face(iCell, local_f);
			
			for (int iVar = 0; iVar < nVar; ++iVar)
				RHS[iCell][iVar] += -r1_v*tria.n_sign(iCell, local_f)*F[global_f][iVar]*tria.areaf(global_f);
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

void HyPE_2D::solve() {

	while (time < Params.final_time) {
		
		compute_rhs(); 
		
		printf ("time = %4.3e, dt = %4.3e, final time = %4.3e\n", time, dt, Params.final_time);
		
		for (int iCell = 0; iCell < tria.no_cells(); ++iCell)
			for (int iVar = 0; iVar < nVar; ++iVar)
				U[iCell][iVar][0] += dt*RHS[iCell][iVar]; 
	
		time += dt; 
	}
	
	printf ("time = %4.3e, dt = %4.3e, final time = %4.3e\n", time, dt, Params.final_time);
	
	std::cout << "Done !" << std::endl;
	std::cout << "----------------------------------------------------------------------------" << std::endl;
}

//----------------------------------------------------------------------------
// Plot solution in vtk format 
//----------------------------------------------------------------------------

void HyPE_2D::plot_vtk(int i, const int digits) {

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

void HyPE_2D::run() {

	initialize(); 
	assign_boundary_ids();
	solve();
	plot_vtk();
	
}
