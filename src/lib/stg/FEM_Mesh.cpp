#ifndef __FINITE_ELEMENTS_MESH__
#define __FINITE_ELEMENTS_MESH__

#include "FEM_Mesh.hpp"

FEM::Mesh::Mesh(size_t Nvert, size_t Nelem, 
  const std::vector<double>& vertices,
  const std::vector<std::vector<size_t>>& cells, 
  const std::vector<size_t>& cell_types
) : Nvert(Nvert), Nelem(Nelem),
    vertices(vertices),
    cells(cells),
    cell_types(cell_types)
{ }

size_t FEM::Mesh::get_vertices_number() const
{
  return Nvert;
}

size_t FEM::Mesh::get_elements_number() const
{
  return Nelem;
}

const std::vector<double>& FEM::Mesh::get_vertices() const
{
  return vertices;
}

const std::vector<std::vector<size_t>>& FEM::Mesh::get_cells() const
{
  return cells;
}

const std::vector<size_t>& FEM::Mesh::get_cell_types() const
{
  return cell_types;
}

const double* FEM::Mesh::get_vertex(size_t i) const
{
  return &vertices[3 * i];
}

void FEM::Mesh::save_mesh_vtk(const std::string& filename) const
{
  size_t rel_nodes = number_of_related_nodes();
  std::ofstream File;
  File.open(filename);	

  // --------- Header of vtk file ---------
  // required version of vtk datafile
  File << "# vtk DataFile Version 3.0" << std::endl; 
  // name of the dataset
  File << "Finite Elements Method" << std::endl;
  // format of dataset
  File << "ASCII" << std::endl;
  // --------- Definition of the FEM type ---------
  File << "DATASET UNSTRUCTURED_GRID" << std::endl;
  File << std::endl;
  // --------- Part to define finite elements mesh ---------
  File << "POINTS " << Nvert << " " << "double" << std::endl;

  for (size_t i = 0; i < Nvert; i++)
  {
    for (size_t j = 0; j < 3; j++)
    {
      File << vertices[3 * i + j] << " ";	
    }
    File << std::endl;	
  }
  
  // --------- Part to define cells --------- Nelem + boundary_elements.size()
  File << std::endl;
  File << "CELLS " << cells.size() << " " << rel_nodes << std::endl; 

  for (const std::vector<size_t>& cell : cells)
  {
    File << cell.size() << " ";
    for (size_t number : cell)
    {
      File << number << " ";
    }
    File << std::endl;
  }
  // -------- Looping through elements to define their type ---------
  File << std::endl;
  File << "CELL_TYPES " << cells.size() << std::endl;
  for (size_t cell_type : cell_types)
  {
    File << cell_type << std::endl;
  }
}

FEM::Mesh::~Mesh() { }

size_t FEM::Mesh::number_of_related_nodes() const
{
  size_t result = 0;
  for (const std::vector<size_t>& cell : cells)
  {
    result += cell.size();
  }
  return result + Nelem;
}

#endif