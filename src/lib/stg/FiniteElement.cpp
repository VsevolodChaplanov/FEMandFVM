#include "FiniteElement.hpp"

FEM::IFiniteElement::IFiniteElement(
  const std::vector<double> &vertices, 
  const std::vector<size_t> global_indices
) : global_indices(global_indices)
{ }

FEM::IFiniteElement::~IFiniteElement() { }

FEM::LinearLineElement::LinearLineElement(
  const std::vector<double> &vertices, 
  const std::vector<size_t> global_indices
) : IFiniteElement(vertices, global_indices),
    det_j(vertices[3] - vertices[0]),
    lumped_matrix({det_j / 2.0, det_j / 2.0}),
    mass_matrix({det_j / 3.0, det_j / 6.0, det_j / 6.0, det_j / 3.0}),
    stiffness_matrix({1.0 / det_j, - 1.0 / det_j, - 1.0 / det_j, 1.0 / det_j}),
    basis_functions( {
        [](const double* param_point)->double{ return 1 - param_point[0]; },
        [](const double* param_point)->double{ return param_point[0]; }
      } ),
    base_point(&vertices[0]),
    center_point(calc_center_point(vertices))
{ }

const std::array<double, 1>& FEM::LinearLineElement::calc_center_point(const std::vector<double>& vertices) const {
  std::array<double, 1> center = { 0. };
  center[0] = (vertices[3] - vertices[0]) / 2;
  return center;
}

double FEM::LinearLineElement::get_mass(size_t i, size_t j) const {
  return mass_matrix.at(i * 2 + j);
}

double FEM::LinearLineElement::get_stiffness(size_t i, size_t j) const {
  return stiffness_matrix.at(i * 2 + j);
}

double FEM::LinearLineElement::get_lumped(size_t i) const {
  return lumped_matrix.at(i);
}

void FEM::LinearLineElement::phys_to_param(const double* phys_in, double* param_out) const {
  param_out[0] = (phys_in[0] - base_point[0]) / det_j;
}

void FEM::LinearLineElement::param_to_phys(const double* param_in, double* phys_out) const {
  phys_out[0] = param_in[0] * det_j + base_point[0];
}

const double* FEM::LinearLineElement::get_center_coordinates() const {
  return &center_point[0];
}

double FEM::LinearLineElement::get_volume() const {
  return det_j;
}

size_t FEM::LinearLineElement::get_number_basis_func() const {
  return 2;
}

size_t FEM::LinearLineElement::get_element_type() const {
  return 1;
}

const std::vector<size_t>& FEM::LinearLineElement::get_global_indices() const {
  return global_indices;
}

FEM::LinearLineElement::~LinearLineElement() { 
  base_point = nullptr;
}

FEM::LinearTriangleElement::LinearTriangleElement(
  const std::vector<double> &vertices, 
  const std::vector<size_t> global_indices
) : IFiniteElement(vertices, global_indices),
    jacobi_matrix({ vertices[3] - vertices[0], vertices[6] - vertices[0],
			vertices[4] - vertices[1], vertices[7] - vertices[1]}),
    det_j(jacobi_matrix[0] * jacobi_matrix[3] - jacobi_matrix[1] * jacobi_matrix[2]),
    lumped_matrix({det_j/6, det_j/6, det_j/6}),
    mass_matrix({ det_j / 12, det_j / 24, det_j / 24,
			det_j / 24, det_j / 12, det_j / 24,
			det_j / 24, det_j / 24, det_j / 12 }),
    stiffness_matrix({
      1/2.0*((jacobi_matrix[0]-jacobi_matrix[1])*(jacobi_matrix[0]-jacobi_matrix[1])+(jacobi_matrix[2]-jacobi_matrix[3])*(jacobi_matrix[2]-jacobi_matrix[3]))/det_j, 
      1/2.0*((jacobi_matrix[0]-jacobi_matrix[1])*jacobi_matrix[1]+(jacobi_matrix[2]-jacobi_matrix[3])*jacobi_matrix[3])/det_j, 
      1/2.0*(jacobi_matrix[0]*(-jacobi_matrix[0]+jacobi_matrix[1])+jacobi_matrix[2]*(-jacobi_matrix[2]+jacobi_matrix[3]))/det_j,

			1/2.0*((jacobi_matrix[0]-jacobi_matrix[1])*jacobi_matrix[1]+(jacobi_matrix[2]-jacobi_matrix[3])*jacobi_matrix[3])/det_j, 
      1/2.0*(jacobi_matrix[1]*jacobi_matrix[1]+jacobi_matrix[3]*jacobi_matrix[3])/det_j, 
      1/2.0*(-jacobi_matrix[0]*jacobi_matrix[1]-jacobi_matrix[2]*jacobi_matrix[3])/det_j,

			1/2.0*(jacobi_matrix[0]*(-jacobi_matrix[0]+jacobi_matrix[1])+jacobi_matrix[2]*(-jacobi_matrix[2]+jacobi_matrix[3]))/det_j, 
      1/2.0*(-jacobi_matrix[0]*jacobi_matrix[1]-jacobi_matrix[2]*jacobi_matrix[3])/det_j, 
      1/2.0*(jacobi_matrix[0]*jacobi_matrix[0]+jacobi_matrix[2]*jacobi_matrix[2])/det_j}),
    basis_functions( {
        [](const double* param_point)->double{ return 1 - param_point[0] - param_point[1]; },
        [](const double* param_point)->double{ return param_point[0] ; },
        [](const double* param_point)->double{ return param_point[1] ; }
      } ),
    base_point(&vertices[0]),
    center_point(calc_center_point(vertices))
{ }

const std::array<double, 2>& FEM::LinearTriangleElement::calc_center_point(const std::vector<double>& vertices) const {
  double area = det_j / 2;
  // Тут утечка памяти и я не знаю как её решить чтобы поле класса осталось константным 
  std::array<double, 2> center_coordinates = { 0., 0. };
  for (size_t i = 0; i < global_indices.size(); i++) {
    size_t x_ind = i * 3;
    size_t y_ind = i * 3 + 1;
    size_t x_ind_next = x_ind + 3;
    size_t y_ind_next = y_ind + 3;
    double temp_diff = ( vertices[x_ind] * vertices[y_ind_next] - 
        vertices[x_ind_next] * vertices[y_ind]
      );
    center_coordinates[0] += (vertices[x_ind] + vertices[x_ind_next]) * temp_diff;
    center_coordinates[1] += (vertices[y_ind] + vertices[y_ind_next]) * temp_diff;
  }
  center_coordinates[0] /= 6 * area;
  center_coordinates[1] /= 6 * area;
  return center_coordinates;
}

double FEM::LinearTriangleElement::get_mass(size_t i, size_t j) const {
  return mass_matrix.at(i * 3 + j);
}

double FEM::LinearTriangleElement::get_stiffness(size_t i, size_t j) const {
  return stiffness_matrix.at(i * 3 + j);
}

double FEM::LinearTriangleElement::get_lumped(size_t i) const {
  return lumped_matrix.at(i);
}

void FEM::LinearTriangleElement::phys_to_param(const double* phys_in, double* param_out) const {
  param_out[0] = (phys_in[0] - base_point[0]) / det_j;
  param_out[2] = (phys_in[1] - base_point[1]) / det_j;
}

void FEM::LinearTriangleElement::param_to_phys(const double* param_in, double* phys_out) const {
  phys_out[0] = param_in[0] * det_j + base_point[0];
  phys_out[1] = param_in[1] * det_j + base_point[1];
}

const double* FEM::LinearTriangleElement::get_center_coordinates() const {
  return &center_point[0];
}

double FEM::LinearTriangleElement::get_volume() const {
  return det_j / 2;
}

size_t FEM::LinearTriangleElement::get_number_basis_func() const {
  return 3;
}

size_t FEM::LinearTriangleElement::get_element_type() const {
  return 2;
}

const std::vector<size_t>& FEM::LinearTriangleElement::get_global_indices() const {
  return global_indices;
}

FEM::LinearTriangleElement::~LinearTriangleElement() { 
  base_point = nullptr;
}

std::unique_ptr<FEM::IFiniteElement> FEM::FiniteElementFactory::create_element(
  const std::vector<double>& vertices,
  const std::vector<size_t>& global_indices,
  size_t element_type
) {
  if (element_type == 1) {
    return std::make_unique<FEM::LinearLineElement>(vertices, global_indices);
  }
  if (element_type == 2) {
    return std::make_unique<FEM::LinearTriangleElement>(vertices, global_indices);
  }
  throw std::runtime_error("Can't assemble element of this element type");
}