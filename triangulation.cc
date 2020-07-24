#include "triangulation.hh"

//----------------------------------------------------------------------------
// Constructor for the triangulation class. Takes the name of the grid file 
// as input. Allocates appropriate memory and forms the grid connectivity
//----------------------------------------------------------------------------

Triangulation::Triangulation(const std::string& file_name) :
	dim(2),
	n_vrtx_per_cell(3),
	n_face_per_cell(3),
	n_vrtx_per_face(2)
{

	// First read the file 
	
	std::ifstream input_data(file_name);

	if ( !(input_data.is_open()) ) {
		std::cerr << "Error. Unable to open file: " << file_name << std::endl; 
		std::exit(EXIT_FAILURE);
	}
	
	std::string line;   // Read a line in the text file
	std::string buffer; // Read a buffer in the text file

	std::getline(input_data, line);
	std::getline(input_data, line);

	// Get the number of cells in the mesh and vertices that make up each cell 
	
	std::istringstream instr(line);
	instr >> buffer >> n_cell;
	
	Utilities::allocate_mem_2d(cell_to_vrtx, n_cell, n_vrtx_per_cell);
	
	for (int c = 0; c < n_cell; ++c) {
		std::getline(input_data, line);
		instr.clear();
		instr.str(line);
		instr >> buffer 
			  >> cell_to_vrtx[c][0]
		      >> cell_to_vrtx[c][1] 
		      >> cell_to_vrtx[c][2];

	
		if ( buffer != "5" ) {
			std::cerr << "Error. Only triangular cells supported. Some cells are not triangular in the mesh" << std::endl; 
			std::exit(EXIT_FAILURE);
		}
	}
	
	// Get the number of vertices in the triangulation
	
	std::getline(input_data, line);
	instr.clear(); instr.str(line);
	instr >> buffer >> n_vrtx;
	
	// Initialize the vertex coordinate matrix and read the grid coordinates
	
	Utilities::allocate_mem_2d(vrtx_coords, n_vrtx, dim);
	
	for (int v = 0; v < n_vrtx; ++v) {
		std::getline(input_data, line);
		instr.clear();
		instr.str(line);
		instr >> vrtx_coords[v][0]
		      >> vrtx_coords[v][1] ;
	}
	
	input_data.close();
	
	// Find the volume (area) of each cell, incircle diameter and cell center coordinates 
	
	Utilities::allocate_mem_1d(cell_vol, n_cell);
	Utilities::allocate_mem_1d(dia_incircle, n_cell);
	Utilities::allocate_mem_2d(cell_coords, n_cell, dim);
	const double r1_3 = 1./3.;
	
	for (int c = 0; c < n_cell; ++c) {
		
		cell_vol[c] = area_of_triangle(vrtx_coords[cell_to_vrtx[c][0]][0], vrtx_coords[cell_to_vrtx[c][0]][1],
									   vrtx_coords[cell_to_vrtx[c][1]][0], vrtx_coords[cell_to_vrtx[c][1]][1],
									   vrtx_coords[cell_to_vrtx[c][2]][0], vrtx_coords[cell_to_vrtx[c][2]][1]);
		
		dia_incircle[c] = incircle_dia(vrtx_coords[cell_to_vrtx[c][0]][0], vrtx_coords[cell_to_vrtx[c][0]][1],
									   vrtx_coords[cell_to_vrtx[c][1]][0], vrtx_coords[cell_to_vrtx[c][1]][1],
									   vrtx_coords[cell_to_vrtx[c][2]][0], vrtx_coords[cell_to_vrtx[c][2]][1]);
		
		cell_coords[c][0] = r1_3*(vrtx_coords[cell_to_vrtx[c][0]][0] + 
								  vrtx_coords[cell_to_vrtx[c][1]][0] + 
								  vrtx_coords[cell_to_vrtx[c][2]][0] ); 
			
		cell_coords[c][1] = r1_3*(vrtx_coords[cell_to_vrtx[c][0]][1] + 
								  vrtx_coords[cell_to_vrtx[c][1]][1] + 
								  vrtx_coords[cell_to_vrtx[c][2]][1] ); 
	}
	
	// Identify the faces in the mesh 
	
	std::pair<int, int>  face1;
	std::pair<int, int>  face2;
	std::pair<int, int>  face3;

	std::pair<int, int>  unordered_face1;
	std::pair<int, int>  unordered_face2;
	std::pair<int, int>  unordered_face3;

	std::set< std::pair<int, int> > faces;

	for (int c = 0; c < n_cell; ++c) {

		// get the three faces of the triangle

		face1.first  = cell_to_vrtx[c][0];
		face1.second = cell_to_vrtx[c][1];

		face2.first  = cell_to_vrtx[c][1];
		face2.second = cell_to_vrtx[c][2];

		face3.first  = cell_to_vrtx[c][2];
		face3.second = cell_to_vrtx[c][0];

		unordered_face1.first  = std::min(face1.first, face1.second);
		unordered_face1.second = std::max(face1.first, face1.second);

		unordered_face2.first  = std::min(face2.first, face2.second);
		unordered_face2.second = std::max(face2.first, face2.second);

		unordered_face3.first  = std::min(face3.first, face3.second);
		unordered_face3.second = std::max(face3.first, face3.second);

		faces.insert(unordered_face1);
		faces.insert(unordered_face2);
		faces.insert(unordered_face3);
	}
	
	n_face = faces.size();
	
	// Build face->vertex link
	
	Utilities::allocate_mem_2d(face_to_vrtx, n_face, n_vrtx_per_face);
	
	int f_index = 0;
	std::set< std::pair<int, int> >::iterator  it;
	
	std::map< std::pair<int, int>, int > face_index_map;
	
	for (it=faces.begin(); it!=faces.end(); ++it) {
		
		face1 = *it;
		
		face_to_vrtx[f_index][0] = face1.first;
		face_to_vrtx[f_index][1] = face1.second;
		
		face_index_map[face1] = f_index;
		
		f_index++;
	}
	
	// Build cell->face link
	
	Utilities::allocate_mem_2d(cell_to_face, n_cell, n_face_per_cell);
	
	for (int c = 0; c < n_cell; ++c) {
		
		face1.first  = cell_to_vrtx[c][0];
		face1.second = cell_to_vrtx[c][1];

		face2.first  = cell_to_vrtx[c][1];
		face2.second = cell_to_vrtx[c][2];

		face3.first  = cell_to_vrtx[c][2];
		face3.second = cell_to_vrtx[c][0];

		unordered_face1.first  = std::min(face1.first, face1.second);
		unordered_face1.second = std::max(face1.first, face1.second);

		unordered_face2.first  = std::min(face2.first, face2.second);
		unordered_face2.second = std::max(face2.first, face2.second);

		unordered_face3.first  = std::min(face3.first, face3.second);
		unordered_face3.second = std::max(face3.first, face3.second);
		
		cell_to_face[c][0] = face_index_map[unordered_face1];
		cell_to_face[c][1] = face_index_map[unordered_face2];
		cell_to_face[c][2] = face_index_map[unordered_face3];
	}
	
	// Build face->cell link 
	
	Utilities::allocate_mem_2d(face_to_cell, n_face, 2);

	bool* did_not_count_face; Utilities::allocate_mem_1d(did_not_count_face, n_face); 
	
	for (int f = 0; f < n_face; ++f) {
		did_not_count_face[f] = true; 
		for (int s = 0; s < 2; ++s) {
			face_to_cell[f][s] = -1; // Initialize all to -1
		}
	}
	
	for (int c = 0; c < n_cell; ++c) {
		for (int f = 0; f < n_face_per_cell; ++f) {
			
			f_index = cell_to_face[c][f]; 
			
			if (did_not_count_face[f_index])
				face_to_cell[f_index][0] = c;
			else 
				face_to_cell[f_index][1] = c;
			
			did_not_count_face[f_index] = false;
		}
	}
	
	Utilities::free_mem_1d(did_not_count_face); 
	
	// Tag boundary faces 
	
	Utilities::allocate_mem_1d(boundary_face, n_face); 
	Utilities::allocate_mem_1d(boundary_id, n_face); 
	
	for (int f = 0; f < n_face; ++f) {
		
		boundary_id[f] = transmissive; 
		
		if (face_to_cell[f][1] == -1)
			boundary_face[f] = true; 
		else
			boundary_face[f] = false;
	}
	
	// Tag boundary cells
	
	Utilities::allocate_mem_1d(boundary_cell, n_cell); 
	
	for (int c = 0; c < n_cell; ++c) {
		boundary_cell[c] = false; 
		for (int f =0; f < n_face_per_cell; ++f) {
			if (boundary_face[cell_to_face[c][f]])
				boundary_cell[c] = true;
		}
	}

	// Tag boundary vertices 

	Utilities::allocate_mem_1d(boundary_vrtx, n_vrtx);

	for (int v = 0; v < n_vrtx; ++v) {
		boundary_vrtx[v] = false; // First assume no vertex is on boundary. 
	}
	
	
	for (int f = 0; f < n_face; ++f) {
		
		// If the face is on boundary, both the vertices are on its boundary. 
		
		if (boundary_face[f]) {
			boundary_vrtx[face_to_vrtx[f][0]] = true;
			boundary_vrtx[face_to_vrtx[f][1]] = true; 
		}
	}

	// Build cell neighbour link 

	Utilities::allocate_mem_2d(cell_nb, n_cell, GeometryInfo::faces_per_cell);

	for (int c = 0; c < n_cell; ++c) {
		for (int n_f = 0; n_f < GeometryInfo::faces_per_cell; ++n_f) {
			
			f_index = cell_to_face[c][n_f];
			
			if (face_to_cell[f_index][0] != c)
				cell_nb[c][n_f] = face_to_cell[f_index][0];
			else 
				cell_nb[c][n_f] = face_to_cell[f_index][1];
		}
	}

	// Calculate the area of each face and surface normal center 
	
	Utilities::allocate_mem_1d(face_area, n_face);
	Utilities::allocate_mem_2d(face_coords, n_face, dim);
	Utilities::allocate_mem_2d(sn, n_face, dim);
	
	double xc1, xc2, yc1, yc2, dx, dy,dxc, dyc; 
	
	for (int f = 0; f < n_face; ++f) {
		
		face_area[f] = area_of_face(vrtx_coords[face_to_vrtx[f][0]][0], vrtx_coords[face_to_vrtx[f][0]][1],
							        vrtx_coords[face_to_vrtx[f][1]][0], vrtx_coords[face_to_vrtx[f][1]][1]); 
		
		face_coords[f][0] = 0.5*(vrtx_coords[face_to_vrtx[f][0]][0] + vrtx_coords[face_to_vrtx[f][1]][0]);
		face_coords[f][1] = 0.5*(vrtx_coords[face_to_vrtx[f][0]][1] + vrtx_coords[face_to_vrtx[f][1]][1]);
		
		if (boundary_face[f]) {
			xc1 = cell_coords[face_to_cell[f][0]][0]; yc1 = cell_coords[face_to_cell[f][0]][1];
			xc2 = face_coords[f][0]; yc2 = face_coords[f][1]; 
		}
		
		else {
			xc1 = cell_coords[face_to_cell[f][0]][0]; yc1 = cell_coords[face_to_cell[f][0]][1]; 
			xc2 = cell_coords[face_to_cell[f][1]][0]; yc2 = cell_coords[face_to_cell[f][1]][1];
		}
		
		dx = vrtx_coords[face_to_vrtx[f][1]][0] - vrtx_coords[face_to_vrtx[f][0]][0];
		dy = vrtx_coords[face_to_vrtx[f][1]][1] - vrtx_coords[face_to_vrtx[f][0]][1];
		
		dxc = xc2 - xc1; dyc = yc2 - yc1; 
		
		sn[f][0] = -dy/face_area[f]; sn[f][1] = dx/face_area[f]; 
		
		if ( (sn[f][0]*dxc + sn[f][1]*dyc) < 0.0) {
			sn[f][0] = -sn[f][0]; 
			sn[f][1] = -sn[f][1];
		}
	}
	
	// Assign a sign to the surface normal from the view point of the cell
	
	Utilities::allocate_mem_2d(sn_sign, n_cell, n_face_per_cell);
	
	int gf, lc; 
	
	for (int c = 0; c < n_cell; ++c) {
		for (int f =0; f < n_face_per_cell; ++f) {
			gf = cell_to_face[c][f]; 
			lc = face_to_cell[gf][0];
			
			if (lc == c)
				sn_sign[c][f] =  1.0;
			else
				sn_sign[c][f] = -1.0;
		}
	}
	
	// Vertex to cell connectivity 
	
	Utilities::allocate_mem_1d(vrtx_to_cell1, n_cell*3);
	Utilities::allocate_mem_1d(vrtx_to_cell2, n_vrtx + 1);
	
	int ipoin, ipoi1, istor; 
	
	// Initialize
	
	for (int v = 0; v < n_vrtx+1; ++v) {
		vrtx_to_cell2[v]  = 0; 
	}
	
	for (int c = 0; c < n_cell; ++c) { 		
		for (int v = 0; v < n_vrtx_per_cell; ++v) {
			ipoi1 = cell_to_vrtx[c][v] + 1;
			vrtx_to_cell2[ipoi1] = vrtx_to_cell2[ipoi1] + 1;
		}
	}
	
	// Storage/reshuffling pass 1 
	
	for (int v = 1; v < n_vrtx + 1; ++v) {
		// Loop over all vertices 
		// Update storage counter and store
		vrtx_to_cell2[v] = vrtx_to_cell2[v] + vrtx_to_cell2[v-1]; 
	}
	
	// Cell pass 2: Store the cells in vert_to_cell1

	for (int c = 0; c < n_cell; ++c) { 		
		for (int v = 0; v < n_vrtx_per_cell; ++v) {
			// Update storage counter, storing in vert_to_cell1
			ipoin = cell_to_vrtx[c][v]; 
			istor = vrtx_to_cell2[ipoin] + 1;   
			vrtx_to_cell2[ipoin] = istor; 
			vrtx_to_cell1[istor-1] = c; 
		}
	}
	
	// Storage/reshuffling pass 2:
	
	for (int v = n_vrtx; v >=1 ; --v) 
		vrtx_to_cell2[v] = vrtx_to_cell2[v-1]; 
	
	vrtx_to_cell2[0] = 0; 
	
	// Subtract 1 from all 
	
	for (int v = 0; v < n_vrtx+1; ++v) 
		vrtx_to_cell2[v]  -= 1; 
}

//----------------------------------------------------------------------------
// Destructor: free all the allocated memory of the class 
//----------------------------------------------------------------------------

Triangulation::~Triangulation() {
	
	
	Utilities::free_mem_2d(vrtx_coords, n_vrtx);
	Utilities::free_mem_2d(cell_coords, n_cell);
	Utilities::free_mem_2d(face_coords, n_face);
	Utilities::free_mem_2d(cell_to_vrtx, n_cell);
	Utilities::free_mem_1d(cell_vol);
	Utilities::free_mem_1d(dia_incircle);
	Utilities::free_mem_2d(face_to_vrtx, n_face);
	Utilities::free_mem_2d(cell_to_face, n_cell);
	Utilities::free_mem_2d(face_to_cell, n_face);
	Utilities::free_mem_2d(cell_nb, n_cell);
	Utilities::free_mem_1d(face_area);
	Utilities::free_mem_1d(boundary_face);
	Utilities::free_mem_1d(boundary_cell);
	Utilities::free_mem_1d(boundary_vrtx); 
	Utilities::free_mem_1d(boundary_id);
	Utilities::free_mem_2d(sn, n_face);
	Utilities::free_mem_2d(sn_sign, n_cell);
	Utilities::free_mem_1d(vrtx_to_cell1);
	Utilities::free_mem_1d(vrtx_to_cell2);
}

//----------------------------------------------------------------------------
// Set the boundary id for a face. If the face is not at boundary, throw an
// error 
//----------------------------------------------------------------------------

void Triangulation::set_boundary_cond(const int& f, const boundary_condition& b_cond) {
	
	if (!(boundary_face[f])) {
		std::cerr << "Error. You are trying to set boundary id for a face which is not at boundary" << std::endl; 
		std::exit(EXIT_FAILURE);
	}
	
	boundary_id[f] = b_cond; 
}

//----------------------------------------------------------------------------
// Find the minimum volume ratio of cells 
//----------------------------------------------------------------------------

double Triangulation::min_vol_ratio() const {

	double min = cell_vol[0]; double max = cell_vol[0];
	
	for (int c = 1; c < n_cell; ++c) {
		min = std::min(min, cell_vol[c]); 
		max = std::max(max, cell_vol[c]);
	}
	
	return min/max; 
}

//----------------------------------------------------------------------------
// Find the minimum incircle diameter of the triangulation  
//----------------------------------------------------------------------------

double Triangulation::min_incircle_dia() const {

	double min = dia_incircle[0];
	
	for (int c = 1; c < n_cell; ++c) 
		min = std::min(min, dia_incircle[c]); 
	
	return min; 
}


//----------------------------------------------------------------------------
// Given the three coordinates of a triangle, find its area
//----------------------------------------------------------------------------

double area_of_triangle(const double& x1, const double& y1,
						const double& x2, const double& y2,
						const double& x3, const double& y3) {
	
	// use shoe-lace formula
	
	return 0.5*std::abs(x1*y2 + x2*y3 + x3*y1 - x1*y3 - x3*y2 - x2*y1);
}

//----------------------------------------------------------------------------
// Find the diameter of the incircle of a triangle 
//----------------------------------------------------------------------------

double incircle_dia(const double& x1, const double& y1,
				    const double& x2, const double& y2,
				    const double& x3, const double& y3) {
	
	// use heron's formula
	
	double a = std::hypot(x2-x1, y2-y1);
	double b = std::hypot(x3-x2, y3-y2);
	double c = std::hypot(x1-x3, y1-y3);
	double s = 0.5*(a + b + c); 
	
	double r = std::sqrt(s*(s-a)*(s-b)*(s-c))/s; 
	
	return 2.0*r;
}

//----------------------------------------------------------------------------
// Given the two coordinates of a face, find its length
//----------------------------------------------------------------------------

double area_of_face(const double& x1, const double& y1,
					const double& x2, const double& y2) {
	
	return std::hypot(x2-x1, y2-y1); 
}

//----------------------------------------------------------------------------
// Integrate the function (x-x_0) on a given cell 
//----------------------------------------------------------------------------

double integrate_xmx0(const double& x0, 
	                  const double& x1, const double& y1,
				      const double& x2, const double& y2,
				      const double& x3, const double& y3) {
	
	double term = (1./2.)*x0*(-x1*y2 + x1*y3 + x2*y1 - x2*y3 - x3*y1 + x3*y2) + 
	              (1./6.)*(x1*(x1*y2 - x1*y3 - x2*y1 + x2*y2 + x3*y1 - x3*y3) + 
	                       x2*(-x2*y1 + x2*y3 - x3*y2 + x3*y3) + 
						   x3*(x3*y1 - x3*y2));
	return term;
}

//----------------------------------------------------------------------------
// Integrate the function (y-y_0) on a given cell 
//----------------------------------------------------------------------------

double integrate_ymy0(const double& y0, 
	                  const double& x1, const double& y1,
				      const double& x2, const double& y2,
				      const double& x3, const double& y3) {
	

	double term = (1./2.)*y0*(-x1*y2 + x1*y3 + x2*y1 - x2*y3 - x3*y1 + x3*y2) + 
	              (1./6.)*(y1*(x1*y2 - x1*y3 - x2*y1 - x2*y2 + x3*y1 + x3*y3) + 
				           y2*(x1*y2 + x2*y3 - x3*y2 - x3*y3) + 
						   y3*(x2*y3 - x1*y3));
	return term;
}

//----------------------------------------------------------------------------
// Given a point A(x1, y1) and a line defined by L((x2, y2), (x3, y3)), find 
// the mirror image of the point A wrt to the line L
//----------------------------------------------------------------------------

void reflect_point(const double& x1, const double& y1,
				   const double& x2, const double& y2,
				   const double& x3, const double& y3,
				   double& xr, double& yr) {
	double a, b, c, rhs;  

	// If the line is parallel to y-axis, directly use geometric result 

	if (std::abs(x3-x2) < 1.0e-12) {
		xr = x1 - 2.0*(x1-x3);
		yr = y1; 
	}

	// Else, use formula from coordinate geometry 

	else {
		a = (y3-y2)/(x3-x2);
		b = -1.0; 
		c = y2 - a*x2; 
		rhs = -2.0*(a*x1 + b*y1 + c)/(a*a + b*b);

		xr = x1 + a*rhs;
		yr = y1 + b*rhs;  
	}
}
