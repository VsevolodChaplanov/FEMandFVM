#ifndef __FINITE_ELEMENTS_MESH__
#define __FINITE_ELEMENTS_MESH__

#include <iostream>
#include <algorithm>
#include <string>
#include <vector>
#include <fstream>

class Mesh {
private:
  // Number of vertices of mesh
  size_t Nvert;
  // Number of elements of mesh
  size_t Nelem;
  // Vertices 
  // {x0, y0, z0, x1, y1, z1, ...}
  std::vector<double> vertices;
  // as in vtk format
  // {
  //  {number_of_vertices vert1 vert2 vert3 ...},
  //  ... 
  // }
  // number_of_vertices is not included
  std::vector<std::vector<size_t>> cells;
  // Types of finite elements
  std::vector<size_t> cell_types;
  // Return number of related nodes for
  // CELLS Nelem N_related_nodes
  size_t number_of_related_nodes() const;
public:
  Mesh(size_t, size_t, const std::vector<double>&,
      const std::vector<std::vector<size_t>>&, const std::vector<size_t>&);
  size_t get_vertices_number() const;
  size_t get_elements_number() const;
  const std::vector<double>& get_vertices() const;
  const std::vector<std::vector<size_t>>& get_cells() const;
  const std::vector<size_t>& get_cell_types() const;
  const double* get_vertex(size_t) const;
  void save_mesh_vtk(const std::string& filename) const;
  ~Mesh();
};

#endif
