/*
 * triangulation.hh
 *      Author: sunder
 */  

#ifndef TRIANGULATION_HH_
#define TRIANGULATION_HH_

#include "headers.hh"

struct GeometryInfo {
	static const int vertices_per_cell = 3;
	static const int faces_per_cell = 3;
};

class Triangulation {
	int dim;              /* No. of dimensions - 2 in this case */
	int n_vrtx_per_cell;  /* No. of vertices per cell */
	int n_face_per_cell;  /* No. of faces in each cell */
	int n_vrtx_per_face;  /* No. of faces in each cell */
	int n_vrtx;           /* Total no. of vertices in the mesh */ 
	int n_cell;           /* Total no. of cells in the mesh */
	int n_face;           /* Total no. of faces in the mesh */
	double** vrtx_coords; /* Coordinates of the vertices in the triangulation */
	double** cell_coords; /* Coordinates of the cell centers in the triangulation */
	double** face_coords; /* Coordinates of the face centers in the triangulation */
	int** cell_to_vrtx;   /* Cell to vertex connectivity */
	int** face_to_vrtx;   /* Face to vertex connectivity */
	int** cell_to_face;   /* Cell to face connectivity */
	int** face_to_cell;   /* Face to cell connectivity */
	int* vrtx_to_cell1;   /* Vertex to cell connectivity */
	int* vrtx_to_cell2;   /* Vertex to cell connectivity */
	double* cell_vol;     /* Volume (in this case area) of each cell */ 
	double* face_area;    /* Area of each face */
	double* dia_incircle; /* Diameter of incircle of each cell */
	double** sn;          /* Surface normal of a face */
	double** sn_sign;     /* Sign of the surface normal from cell view point */
	bool* boundary_face;  /* List of all the faces at boundary */
	bool* boundary_cell;  /* List of all the cells at boundary */
	int* boundary_id;     /* Boundary id's for faces on boundary (default 0) */ 
public:
	
	/* Constructor and destructor */
	
	Triangulation (const std::string& file_name);  /* Constructor */
	~Triangulation();                              /* Destructor */
	
	/* Basic information about the grid */
	
	int no_dim() const {return dim;} 
	int no_vertices() const { return n_vrtx; }
	int no_cells() const { return n_cell; }
	int no_faces() const { return n_face; }
	
	/* Connectivity information */
	
	/* For cell c, get the global vertex index of local vertex n_v */ 
	int link_cell_to_vertex(const int& c, const int& n_v) const {return cell_to_vrtx[c][n_v];}
	/* For cell c, get the global face index of local face n_f */
	int link_cell_to_face(const int& c, const int& n_f) {return cell_to_face[c][n_f];}
	/* For face f, get the global face indices of the two cells straddling the face */
	int link_face_to_cell(const int& f, const int& n_c) {return face_to_cell[f][n_c];}
	/* For face f, get the global face indices of the two vertices that make the cell */
	int link_face_to_vertex(const int& f, const int& n_v) {return face_to_vrtx[f][n_v];}
	/* Number of cells that share a vertex */
	int no_cells_sharing_vertex(const int& v) const {return (vrtx_to_cell2[v+1] - (vrtx_to_cell2[v]+1)) + 1;}
	/* I'th cell sharing vertex v */
	int cell_sharing_vertex(const int& v, const int& i) const {return vrtx_to_cell1[vrtx_to_cell2[v] + 1 + i];}
	/* Check if the given cell is at boundary */
	bool cell_at_boundary(const int& c) const {return boundary_cell[c];}
	/* Check if the given face is at boundary */
	bool face_at_boundary(const int& f) const {return boundary_face[f];}
	/* Set the boundary id of face f */
	void set_boundary_id(const int&, const int&);
	/* Get the boundary id of face f */
	int get_boundary_id(const int& f) const {return boundary_id[f];}
	
	/* Metric data */
	
	double xv(const int& v) const { return vrtx_coords[v][0]; }               /* x-coordinate of vertex v */
	double yv(const int& v) const { return vrtx_coords[v][1]; }               /* y-coordinate of vertex v */
	double xc(const int& c) const { return cell_coords[c][0]; }               /* x-coordinate of cell center c */
	double yc(const int& c) const { return cell_coords[c][1]; }               /* y-coordinate of cell center c */
	double xf(const int& f) const { return face_coords[f][0]; }               /* x-coordinate of face center f */
	double yf(const int& f) const { return face_coords[f][1]; }               /* y-coordinate of face center f */
	double vol(const int& c) const {return cell_vol[c];}                      /* Volume of cell c */
	double h(const int& c) const {return dia_incircle[c];}                    /* Incircle diameter of cell c */
	double areaf(const int& f) const {return face_area[f];}                   /* Area of a face f */
	double nx(const int& f) const {return sn[f][0];}                          /* Normal in the x-direction for face f */
	double ny(const int& f) const {return sn[f][1];}                          /* Normal in the y-direction for face f */
	double n_sign(const int& c, const int& f) const {return sn_sign[c][f];}   /* Sign of normal from local view point */
	
	/* Mesh quality */
	
	double min_vol_ratio() const;
	double min_incircle_dia() const;
	
};

// Additional geometry related functions 

double area_of_triangle(const double&, const double&,
						const double&, const double&,
						const double&, const double&); 

double incircle_dia(const double&, const double&,
						const double&, const double&,
						const double&, const double&); 

double area_of_face(const double&, const double&,
				    const double&, const double&); 

#endif /* TRIANGULATION_HH_ */
