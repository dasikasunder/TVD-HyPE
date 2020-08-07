#include "hype.hh"

//----------------------------------------------------------------------------
// Constructor - Takes input parameters and sets memory 
//----------------------------------------------------------------------------

HyPE_2D::HyPE_2D(Vector (*f)(double, double), AppCtx params):
	Params(params),
	tria(params.mesh_file_name),
	init_func(f),
	time(0.0),
	dt(0.0),
	time_step(0)
{
	
	std::cout << "Allocating memory. ";
	std::cout << "No. of cells in the triangulation = " << tria.no_cells() << ". "; 
	
	// Allocate memory for all the variables 
	
	uh.resize(extents[tria.no_cells()][nVar][4]); 
	wh.resize(extents[tria.no_cells()][nVar]);
	duh.resize(extents[tria.no_cells()][nVar]);
	Fu.resize(extents[tria.no_faces()][nVar]);
	qr.resize(extents[tria.no_cells()]);
	
	h_min = tria.min_incircle_dia();

	//-------------------------------------------------------------------------------------
	// Now start computing gradient matrix 
	//-------------------------------------------------------------------------------------

	const int n_nb_cells = 3; 
	int Nb, count; 
	Matrix A(extents[n_nb_cells][nDim]);
	double x0, y0, xr, yr, x[3], y[3], r1_vol_nb;

	for (int iCell = 0; iCell < tria.no_cells(); ++iCell) {
		
		x0 = tria.xc(iCell); y0 = tria.yc(iCell);

		for (int iNb = 0; iNb < n_nb_cells; ++iNb) {

			Nb = tria.cell_neighbour(iCell, iNb);
			count = 1; 
			
			// If the neighbouring cell is present 

			if (Nb != -1) {

				// Get the three vertices of the triangle 

				for (int v = 0; v < GeometryInfo::vertices_per_cell; ++v) {
					x[v] = tria.xv(tria.link_cell_to_vertex(Nb, v)); 
					y[v] = tria.yv(tria.link_cell_to_vertex(Nb, v)); 
				}
				
				r1_vol_nb = 1.0/tria.vol(Nb);
			}	

			// If the neighbouring cell is absent 
			// We will mirror the current cell and construct a ghost cell  
			
			else {

				for (int v = 0; v < GeometryInfo::vertices_per_cell; ++v) {
					
					if( tria.vrtx_at_boundary(tria.link_cell_to_vertex(iCell, v)) ) {
						x[count] =  tria.xv(tria.link_cell_to_vertex(iCell, v));
						y[count] =  tria.yv(tria.link_cell_to_vertex(iCell, v));
						count++; 
					}

					else {
						xr =  tria.xv(tria.link_cell_to_vertex(iCell, v));
						yr =  tria.yv(tria.link_cell_to_vertex(iCell, v));
					}
				}

				reflect_point(xr, yr, x[1], y[1], x[2], y[2], x[0], y[0]); 
				r1_vol_nb = 1.0/tria.vol(iCell); 
			}
			
			A[iNb][0] = r1_vol_nb*integrate_xmx0(x0, x[0], y[0], x[1], y[1], x[2], y[2]); 
			A[iNb][1] = r1_vol_nb*integrate_ymy0(y0, x[0], y[0], x[1], y[1], x[2], y[2]);
		}

		qr[iCell].initialize(A, n_nb_cells, nDim);
	}

	std::cout << "Done !" << std::endl;
}

//----------------------------------------------------------------------------
// Initialize the solution using the initial conditions 
//----------------------------------------------------------------------------

void HyPE_2D::initialize() {
	
	std::cout << "Initializing the solution..."; 

	Vector V(extents[nVar]), Q(extents[nVar]);
	
	for (int iCell = 0; iCell < tria.no_cells(); ++iCell) {
		
		V = init_func(tria.xc(iCell), tria.yc(iCell)); 
		
		PDEPrim2Cons(V, Q); 
		
		for (int iVar = 0; iVar < nVar; ++iVar)
			uh[iCell][iVar][0] = Q[iVar]; 
	}
	
	std::cout << "Done !" << std::endl;
}

//----------------------------------------------------------------------------
// To each face in the mesh assign boundary condition 
//----------------------------------------------------------------------------

void HyPE_2D::set_boundary_conds() {
	
	std::cout << "Setting boundary conditions..."; 
	
	// Loop over all the faces 
	
	for (int iFace = 0; iFace < tria.no_faces(); ++iFace) {
		
		if (tria.face_at_boundary(iFace)) {

			tria.set_boundary_cond(iFace, reflective);

			if (std::abs((tria.nx(iFace) + 1.0)) < 1.0e-12) {
				tria.set_boundary_cond(iFace, transmissive);
			} 
			
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
			Q[iVar] = uh[iCell][iVar][0]; 
		
		PDECons2Prim(Q, V);
		
		for (int iVar = 0; iVar < nVar; ++iVar)
			wh[iCell][iVar] = V[iVar]; 
	}
}

//----------------------------------------------------------------------------
// Reconstruct the solution in each cell using Barth-Jesperson limiter 
//----------------------------------------------------------------------------

void HyPE_2D::reconstruct() {

	const int n_nb_cells = 3; 
	Vector rhs(extents[n_nb_cells]);
	Vector nb_value(extents[n_nb_cells]);
	Vector sol(extents[n_nb_cells]);
	Vector Qf(extents[nVar]);
	Vector Ff(extents[nVar]);
	Vector V_bnd(extents[nVar]);
	Vector Q_bnd(extents[nVar]);
	double min_val, max_val, node_val, x_node, y_node, xf, yf, x0, y0, nx, ny, Af, reduce_modes;
	double temp, r1_v; 
	int Nb, iFace;

	//-----------------------------------------------------
	V_bnd[0] = 1.4;  // Density  
	V_bnd[1] = 3.0;  // x-Velocity 
	V_bnd[2] = 0.0;  // y-Velocity 
	V_bnd[3] = 1.0;  // Pressure 
	//-----------------------------------------------------

	PDEPrim2Cons(V_bnd, Q_bnd);

	for (int iCell = 0; iCell < tria.no_cells(); ++iCell) {

		x0 = tria.xc(iCell); y0 = tria.yc(iCell);
		r1_v = 1.0/tria.vol(iCell);
		
		reduce_modes = 1.0; 

		for (int iVar = 0; iVar < nVar; ++iVar) {

			min_val = uh[iCell][iVar][0]; 
			max_val = uh[iCell][iVar][0]; 
			uh[iCell][iVar][3] = 0.0; 

			
			for (int iNb = 0; iNb < n_nb_cells; ++iNb) {
				
				Nb = tria.cell_neighbour(iCell, iNb);
			
				// If the neighbouring cell is present 

				if (Nb != -1) {
					nb_value[iNb] = uh[Nb][iVar][0];
				}
				
				// If neighbour cell is absent fill appropriate ghost cell value

				else {
					
					iFace = tria.link_cell_to_face(iCell, iNb);
					
					nx = tria.n_sign(iCell, iNb)*tria.nx(iFace); 
					ny = tria.n_sign(iCell, iNb)*tria.ny(iFace);

					if (tria.get_boundary_cond(iFace) == transmissive) {
						nb_value[iNb] = uh[Nb][iVar][0];
					}

					if (tria.get_boundary_cond(iFace) == reflective) {
						
						nb_value[iNb] = uh[Nb][iVar][0];
						
						if (iVar == 1)
							nb_value[iNb] = uh[iCell][iVar][1] - 2.0*uh[iCell][iVar][1]*nx*nx - 2.0*uh[iCell][iVar][2]*nx*ny; 
						
						if (iVar == 2) 
							nb_value[iNb] = uh[iCell][iVar][2] - 2.0*uh[iCell][iVar][1]*nx*ny - 2.0*uh[iCell][iVar][1]*ny*ny;
						
					}

					if (tria.get_boundary_cond(iFace) == inflow) {
						nb_value[iNb] = Q_bnd[iVar];
					}

					nb_value[iNb] = uh[iCell][iVar][0]; // Transmissive boundary 
				}

				// Find the minimum and maximum in the neighbourhood 

				if (nb_value[iNb] > max_val)
					max_val = nb_value[iNb];
				
				if (nb_value[iNb] < min_val)
					min_val = nb_value[iNb];

				rhs[iNb] = nb_value[iNb] - uh[iCell][iVar][0];
			}
			
			// Find the non-limited solution 

			qr[iCell].solve(rhs, sol); 
			
			// Limit the solution using Barth-Jesperson limiter 

			for (int v = 0; v < GeometryInfo::vertices_per_cell; ++v) {

				x_node = tria.xv(tria.link_cell_to_vertex(iCell,v));
				y_node = tria.yv(tria.link_cell_to_vertex(iCell,v));

				node_val = uh[iCell][iVar][0] + sol[0]*(x_node-x0) + sol[1]*(y_node-y0);
				
				temp = node_val - uh[iCell][iVar][0];  

				if (std::abs(temp) < 1.0e-12) {
					reduce_modes = std::min(1.0, reduce_modes);
				}

				else if (temp > 0.0) {
					reduce_modes = std::min(reduce_modes, (max_val-uh[iCell][iVar][0])/(temp+1.0e-12));
				}

				else {
					reduce_modes = std::min(reduce_modes, (min_val-uh[iCell][iVar][0])/(temp+1.0e-12));
				}
			}

			// Make sure reduce modes are greater than zero 

			if (reduce_modes < 0.0)
				reduce_modes = 0.0; 

			// Find the limited solution 

			uh[iCell][iVar][1] = sol[0]; 
			uh[iCell][iVar][2] = sol[1];
		}

		for (int iVar = 0; iVar < nVar; ++iVar) {
			uh[iCell][iVar][1]  = reduce_modes*uh[iCell][iVar][1];
			uh[iCell][iVar][2]  = reduce_modes*uh[iCell][iVar][2];
		}

		// After finding the spatial variation, we now need to find the temporal variation

		for (int f = 0; f < GeometryInfo::faces_per_cell; ++f) {
			
			iFace = tria.link_cell_to_face(iCell, f); 
			xf = tria.xf(iFace); yf = tria.yf(iFace); 
			
			nx = tria.n_sign(iCell, f)*tria.nx(iFace); 
			ny = tria.n_sign(iCell, f)*tria.ny(iFace);

			Af = tria.areaf(iFace);

			for (int iVar = 0; iVar < nVar; ++iVar) {
				Qf[iVar] = uh[iCell][iVar][0] + uh[iCell][iVar][1]*(xf-x0) + uh[iCell][iVar][2]*(yf-y0);
			}

			PDEConsFlux(Qf, nx, ny, xf, yf, Ff); 

			for (int iVar = 0; iVar < nVar; ++iVar) {
				uh[iCell][iVar][3] += -r1_v*Af*Ff[iVar];
			}
		}
	}
} 

//----------------------------------------------------------------------------
// Find the RHS in each cell and also compute the time step size 
//----------------------------------------------------------------------------

void HyPE_2D::compute_rhs() {
	
	Vector QL(extents[nVar]), QR(extents[nVar]), VR(extents[nVar]), Flux(extents[nVar]); 

	double x0_L, y0_L, x0_R, y0_R,  xf, yf, s, s_max = 0.0; 
	
	// Find the upwind flux on each face 
	
	int L_cell_index, R_cell_index;
	
	for (int iFace = 0; iFace < tria.no_faces(); ++iFace) {
		
		L_cell_index = tria.link_face_to_cell(iFace, 0);
		
		xf = tria.xf(iFace);
		yf = tria.yf(iFace); 

		x0_L = tria.xc(L_cell_index); 
		y0_L = tria.yc(L_cell_index); 
		
		for (int iVar = 0; iVar < nVar; ++iVar) {
			QL[iVar] = uh[L_cell_index][iVar][0] + 
			           uh[L_cell_index][iVar][1]*(xf-x0_L) +
					   uh[L_cell_index][iVar][2]*(yf-y0_L) + 
					   uh[L_cell_index][iVar][3]*0.5*dt;
		}
		
		// Apply boundary conditions for faces at boundaries  
		
		if (tria.face_at_boundary(iFace)) {
			
			if (tria.get_boundary_cond(iFace) == inflow) { 
				
				VR[0] = 1.4;  // Density  
				VR[1] = 3.0;  // x-Velocity 
				VR[2] = 0.0;  // y-Velocity 
				VR[3] = 1.0;  // Pressure 
			
				PDEPrim2Cons(VR, QR);
			}
			
			else if (tria.get_boundary_cond(iFace) == reflective) {
				
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

			x0_R = tria.xc(R_cell_index);
			y0_R = tria.yc(R_cell_index); 
			
			for (int iVar = 0; iVar < nVar; ++iVar) {
				
				QR[iVar] = uh[R_cell_index][iVar][0] + 
			               uh[R_cell_index][iVar][1]*(xf-x0_R) +
					       uh[R_cell_index][iVar][2]*(yf-y0_R) + 
					       uh[R_cell_index][iVar][3]*0.5*dt;
			}
		}
		
		s = PDERusanovFlux(QL, QR, tria.nx(iFace), tria.ny(iFace), tria.xf(iFace), tria.yf(iFace), Flux);
		
		if (s > s_max)
			s_max = s; 
		
		for (int iVar = 0; iVar < nVar; ++iVar)
			Fu[iFace][iVar] = Flux[iVar];
	}
	
	// Now find the RHS in each cell 
	
	int global_f; double r1_v; 
	
	for (int iCell = 0; iCell < tria.no_cells(); ++iCell) {
		
		for (int iVar = 0; iVar < nVar; ++iVar)
			duh[iCell][iVar] = 0.0; 
			
		r1_v = 1./tria.vol(iCell);
		
		for (int local_f = 0; local_f < GeometryInfo::faces_per_cell; ++local_f) {
			
			global_f = tria.link_cell_to_face(iCell, local_f);
			
			for (int iVar = 0; iVar < nVar; ++iVar)
				duh[iCell][iVar] += -r1_v*tria.n_sign(iCell, local_f)*Fu[global_f][iVar]*tria.areaf(global_f);
		}
	}
	
	// Compute the time step value 
	
	dt = (Params.CFL*h_min)/s_max;
	
	//  Check size of dt to avoid exceeding output time

	if((time + dt)>Params.tend)
		dt = Params.tend - time;
}

//----------------------------------------------------------------------------
// Plot solution in vtk format 
//----------------------------------------------------------------------------

void HyPE_2D::plot_vtk(int i, const int digits) {

	std::ofstream vtk;
	const std::string filename = "../plots/sol-" + Utilities::int_to_string (i, digits) + ".vtk";
	vtk.open (filename);
	vtk.flags( std::ios::dec | std::ios::scientific );
	vtk.precision(6);

	if ( !( vtk.is_open() ) ) {
		std::cerr << "Error. Unable to open file: " << filename << std::endl; 
		std::exit(EXIT_FAILURE);
	}
	
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
			RHO_vtrx[v] += wh[tria.cell_sharing_vertex(v,i)][0];
			Vx_vtrx[v]  += wh[tria.cell_sharing_vertex(v,i)][1];
			Vy_vtrx[v]  += wh[tria.cell_sharing_vertex(v,i)][2];
			P_vtrx[v]   += wh[tria.cell_sharing_vertex(v,i)][3];
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
	set_boundary_conds();
	compute_rhs(); // Just to set an initial time step size  
	
	// Main loop in time 

	while (time < Params.tend) {
		
		reconstruct(); 
		compute_rhs();

		if (Params.write_interval != 0 ) 
			if (time_step % Params.write_interval == 0)
				plot_vtk(time_step, 5); 
		
		
		printf ("time = %4.3e, dt = %4.3e, final time = %4.3e\n", time, dt, Params.tend);
		
		for (int iCell = 0; iCell < tria.no_cells(); ++iCell)
			for (int iVar = 0; iVar < nVar; ++iVar)
				uh[iCell][iVar][0] += dt*duh[iCell][iVar]; 
	
		time += dt; 
		time_step++; 
	}
	
	printf ("time = %4.3e, dt = %4.3e, final time = %4.3e\n", time, dt, Params.tend);
	
	std::cout << "Done !" << std::endl;
	std::cout << "----------------------------------------------------------------------------" << std::endl;
	
	plot_vtk(time_step, 5); 
	
}
