#include "FEM_FiniteElement.hpp"

IFiniteElement::IFiniteElement(
  const std::vector<double> &vertices, 
  const std::vector<size_t> global_indices
) : global_indices(global_indices)
{ }

IFiniteElement::~IFiniteElement() { }

LinearLineElement::LinearLineElement(
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

std::array<double, 1> LinearLineElement::calc_center_point(const std::vector<double>& vertices) const {
  std::array<double, 1> center = { 0. };
  center[0] = (vertices[3] - vertices[0]) / 2;
  return center;
}

double LinearLineElement::get_mass(size_t i, size_t j) const {
  return mass_matrix.at(i * 2 + j);
}

double LinearLineElement::get_stiffness(size_t i, size_t j) const {
  return stiffness_matrix.at(i * 2 + j);
}

double LinearLineElement::get_lumped(size_t i) const {
  return lumped_matrix.at(i);
}

void LinearLineElement::phys_to_param(const double* phys_in, double* param_out) const {
  param_out[0] = (phys_in[0] - base_point[0]) / det_j;
}

void LinearLineElement::param_to_phys(const double* param_in, double* phys_out) const {
  phys_out[0] = param_in[0] * det_j + base_point[0];
}

const double* LinearLineElement::get_center_coordinates() const {
  return &center_point[0];
}

double LinearLineElement::get_volume() const {
  return det_j;
}

size_t LinearLineElement::get_number_basis_func() const {
  return 2;
}

size_t LinearLineElement::get_element_type() const {
  return 1;
}

const std::vector<size_t>& LinearLineElement::get_global_indices() const {
  return global_indices;
}

LinearLineElement::~LinearLineElement() { 
  base_point = nullptr;
}

LinearTriangleElement::LinearTriangleElement(
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

std::array<double, 2> LinearTriangleElement::calc_center_point(const std::vector<double>& vertices) const {
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

double LinearTriangleElement::get_mass(size_t i, size_t j) const {
  return mass_matrix.at(i * 3 + j);
}

double LinearTriangleElement::get_stiffness(size_t i, size_t j) const {
  return stiffness_matrix.at(i * 3 + j);
}

double LinearTriangleElement::get_lumped(size_t i) const {
  return lumped_matrix.at(i);
}

void LinearTriangleElement::phys_to_param(const double* phys_in, double* param_out) const {
  param_out[0] = (phys_in[0] - base_point[0]) / det_j;
  param_out[2] = (phys_in[1] - base_point[1]) / det_j;
}

void LinearTriangleElement::param_to_phys(const double* param_in, double* phys_out) const {
  phys_out[0] = param_in[0] * det_j + base_point[0];
  phys_out[1] = param_in[1] * det_j + base_point[1];
}

const double* LinearTriangleElement::get_center_coordinates() const {
  return &center_point[0];
}

double LinearTriangleElement::get_volume() const {
  return det_j / 2;
}

size_t LinearTriangleElement::get_number_basis_func() const {
  return 3;
}

size_t LinearTriangleElement::get_element_type() const {
  return 2;
}

const std::vector<size_t>& LinearTriangleElement::get_global_indices() const {
  return global_indices;
}

LinearTriangleElement::~LinearTriangleElement() { 
  base_point = nullptr;
}

// ----------------------------- Rectangle ----------------------------- //
LinearRectangleElement::LinearRectangleElement(
  const std::vector<double> &vertices, 
  const std::vector<size_t> global_indices
) : IFiniteElement(vertices, global_indices),
    det_j((vertices[3] - vertices[0]) * (vertices[10] - vertices[1])),
    lumped_matrix( {(det_j)/4, (det_j)/4, (det_j)/4, (det_j)/4} ),
    mass_matrix( { det_j /36 * 4,	det_j /36 * 2,	det_j /36 * 1,	det_j /36 * 2,
      det_j /36 * 2,	det_j /36 * 4,	det_j /36 * 2,	det_j /36 * 1,
      det_j /36 * 1,	det_j /36 * 2,	det_j /36 * 4,	det_j /36 * 2,
      det_j /36 * 2,	det_j /36 * 1,	det_j /36 * 2,	det_j /36 * 4 } ),
    stiffness_matrix({
        ((vertices[3] - vertices[0]) * (vertices[3] - vertices[0]) + ((vertices[10] - vertices[1]) * (vertices[10] - vertices[1]))) / (det_j * 6) * 2,	
        (((vertices[3] - vertices[0]) * (vertices[3] - vertices[0])) + ((vertices[10] - vertices[1]) * (vertices[10] - vertices[1]))) / (det_j * 6) * (-2+(3 * ((vertices[3] - vertices[0]) * (vertices[3] - vertices[0])))/(((vertices[3] - vertices[0]) * (vertices[3] - vertices[0]))+((vertices[10] - vertices[1]) * (vertices[10] - vertices[1])))),	
        (((vertices[3] - vertices[0]) * (vertices[3] - vertices[0])) + ((vertices[10] - vertices[1]) * (vertices[10] - vertices[1]))) / (det_j * 6) * (-1),	
        (((vertices[3] - vertices[0]) * (vertices[3] - vertices[0])) + ((vertices[10] - vertices[1]) * (vertices[10] - vertices[1]))) / (det_j * 6) * (1-(3 * ((vertices[3] - vertices[0]) * (vertices[3] - vertices[0])))/(((vertices[3] - vertices[0]) * (vertices[3] - vertices[0]))+((vertices[10] - vertices[1]) * (vertices[10] - vertices[1])))),
        (((vertices[3] - vertices[0]) * (vertices[3] - vertices[0])) + ((vertices[10] - vertices[1]) * (vertices[10] - vertices[1]))) / (det_j * 6) * (-2+(3 * ((vertices[3] - vertices[0]) * (vertices[3] - vertices[0])))/(((vertices[3] - vertices[0]) * (vertices[3] - vertices[0]))+((vertices[10] - vertices[1]) * (vertices[10] - vertices[1])))),	
        (((vertices[3] - vertices[0]) * (vertices[3] - vertices[0])) + ((vertices[10] - vertices[1]) * (vertices[10] - vertices[1]))) / (det_j * 6) * 2,	
        (((vertices[3] - vertices[0]) * (vertices[3] - vertices[0])) + ((vertices[10] - vertices[1]) * (vertices[10] - vertices[1]))) / (det_j * 6) * (1-(3 * ((vertices[3] - vertices[0]) * (vertices[3] - vertices[0])))/(((vertices[3] - vertices[0]) * (vertices[3] - vertices[0]))+((vertices[10] - vertices[1]) * (vertices[10] - vertices[1])))),	
        (((vertices[3] - vertices[0]) * (vertices[3] - vertices[0])) + ((vertices[10] - vertices[1]) * (vertices[10] - vertices[1]))) / (det_j * 6) * (-1),
        (((vertices[3] - vertices[0]) * (vertices[3] - vertices[0])) + ((vertices[10] - vertices[1]) * (vertices[10] - vertices[1]))) / (det_j * 6) * (-1),	
        (((vertices[3] - vertices[0]) * (vertices[3] - vertices[0])) + ((vertices[10] - vertices[1]) * (vertices[10] - vertices[1]))) / (det_j * 6) * (1-(3 * ((vertices[3] - vertices[0]) * (vertices[3] - vertices[0])))/(((vertices[3] - vertices[0]) * (vertices[3] - vertices[0]))+((vertices[10] - vertices[1]) * (vertices[10] - vertices[1])))),	
        (((vertices[3] - vertices[0]) * (vertices[3] - vertices[0])) + ((vertices[10] - vertices[1]) * (vertices[10] - vertices[1]))) / (det_j * 6) * 2,	
        (((vertices[3] - vertices[0]) * (vertices[3] - vertices[0])) + ((vertices[10] - vertices[1]) * (vertices[10] - vertices[1]))) / (det_j * 6) * (-2+(3 * ((vertices[3] - vertices[0]) * (vertices[3] - vertices[0])))/(((vertices[3] - vertices[0]) * (vertices[3] - vertices[0]))+((vertices[10] - vertices[1]) * (vertices[10] - vertices[1])))),
        (((vertices[3] - vertices[0]) * (vertices[3] - vertices[0])) + ((vertices[10] - vertices[1]) * (vertices[10] - vertices[1]))) / (det_j * 6) * (1-(3 * ((vertices[3] - vertices[0]) * (vertices[3] - vertices[0])))/(((vertices[3] - vertices[0]) * (vertices[3] - vertices[0]))+((vertices[10] - vertices[1]) * (vertices[10] - vertices[1])))),	
        (((vertices[3] - vertices[0]) * (vertices[3] - vertices[0])) + ((vertices[10] - vertices[1]) * (vertices[10] - vertices[1]))) / (det_j * 6) * (-1),	
        (((vertices[3] - vertices[0]) * (vertices[3] - vertices[0])) + ((vertices[10] - vertices[1]) * (vertices[10] - vertices[1]))) / (det_j * 6) * (-2+(3 * ((vertices[3] - vertices[0]) * (vertices[3] - vertices[0])))/(((vertices[3] - vertices[0]) * (vertices[3] - vertices[0]))+((vertices[10] - vertices[1]) * (vertices[10] - vertices[1])))),	
        (((vertices[3] - vertices[0]) * (vertices[3] - vertices[0])) + ((vertices[10] - vertices[1]) * (vertices[10] - vertices[1]))) / (det_j * 6) * 2 }),
    basis_functions( {
        [](const double* param_point)->double{ return (1 - param_point[0]) * (1 - param_point[1]); },
        [](const double* param_point)->double{ return param_point[0] * (1 - param_point[1]); },
        [](const double* param_point)->double{ return param_point[0] * param_point[1] ; },
        [](const double* param_point)->double{ return (1 - param_point[0]) * param_point[1]; }
      } ),
    base_point(&vertices[0]),
    center_point(calc_center_point(vertices))
{ }

std::array<double, 2> LinearRectangleElement::calc_center_point(const std::vector<double>& vertices) const {
  double area = det_j;
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

double LinearRectangleElement::get_mass(size_t i, size_t j) const {
  return mass_matrix.at(i * 4 + j);
}

double LinearRectangleElement::get_stiffness(size_t i, size_t j) const {
  return stiffness_matrix.at(i * 4 + j);
}

double LinearRectangleElement::get_lumped(size_t i) const {
  return lumped_matrix.at(i);
}

void LinearRectangleElement::phys_to_param(const double* phys_in, double* param_out) const {
  param_out[0] = (phys_in[0] - base_point[0]) / det_j;
  param_out[2] = (phys_in[1] - base_point[1]) / det_j;
}

void LinearRectangleElement::param_to_phys(const double* param_in, double* phys_out) const {
  phys_out[0] = param_in[0] * det_j + base_point[0];
  phys_out[1] = param_in[1] * det_j + base_point[1];
}

const double* LinearRectangleElement::get_center_coordinates() const {
  return &center_point[0];
}

double LinearRectangleElement::get_volume() const {
  return det_j;
}

size_t LinearRectangleElement::get_number_basis_func() const {
  return 4;
}

size_t LinearRectangleElement::get_element_type() const {
  return 3;
}

const std::vector<size_t>& LinearRectangleElement::get_global_indices() const {
  return global_indices;
}

LinearRectangleElement::~LinearRectangleElement() { 
  base_point = nullptr;
}

// ----------------------------- Factory ----------------------------- //

std::shared_ptr<IFiniteElement> FiniteElementFactory::create_element(
  const std::vector<double>& vertices,
  const std::vector<size_t>& global_indices,
  size_t element_type
) {
  if (element_type == 1) {
    return std::make_shared<LinearLineElement>(vertices, global_indices);
  } else if (element_type == 2) {
    return std::make_shared<LinearTriangleElement>(vertices, global_indices);
  } else if (element_type == 3) {
    return std::make_shared<LinearRectangleElement>(vertices, global_indices);
  }
  throw std::runtime_error("Can't assemble element of this element type");
}