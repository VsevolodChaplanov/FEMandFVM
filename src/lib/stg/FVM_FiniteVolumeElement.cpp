#include "FVM_FiniteVolumeElement.hpp"

FVM::FiniteVolumeElement::FiniteVolumeElement(
  const std::vector<double>& vertices, 
  const std::vector<size_t>& global_indices
) : global_indices(global_indices),
    n_vert(global_indices.size()),
    area(calc_area(vertices)),
    center(calc_cell_center(vertices))
{ }

double FVM::FiniteVolumeElement::calc_area(const std::vector<double>& vertices) const {
  double summ = 0;
  for (size_t i = 0; i < n_vert - 1; i++)
  {
    size_t x_ind = i * 3;
    size_t y_ind = i * 3 + 1;
    size_t x_ind_next = x_ind + 3;
    size_t y_ind_next = y_ind + 3;
    summ += vertices[x_ind] * vertices[y_ind_next] - vertices[x_ind_next] * vertices[y_ind];
  }
  summ += vertices[n_vert - 3] * vertices[1] - vertices[n_vert - 3 + 1] * vertices[0];
  return summ / 2.;
}

const double* FVM::FiniteVolumeElement::calc_cell_center(const std::vector<double>& vertices) const {
  double* center_coordinates = new double[2] { 0., 0. };
  for (size_t i = 0; i < n_vert; i++)
  {
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

const double* FVM::FiniteVolumeElement::cell_center() const {
  return center;
}

const std::vector<size_t>& FVM::FiniteVolumeElement::get_global_indices() const {
  return global_indices;
}

FVM::FiniteVolumeElement::~FiniteVolumeElement() {
  delete[] center;
}

FVM::Face::Face(
  const std::vector<double>& vertices, 
  const std::vector<size_t>& global_indices
) : global_indices(global_indices),
    n_vert(global_indices.size()),
    lenght(calc_lenght(vertices)),
    center(calc_face_center(vertices))
{ }

double FVM::Face::calc_lenght(const std::vector<double>& vertices) const {
  double face_length = sqrt( (vertices[3] - vertices[0]) * (vertices[3] - vertices[0]) +
    (vertices[3 + 1] - vertices[1]) * (vertices[3 + 1] - vertices[1]) );
  return face_length;
}

const double* FVM::Face::calc_face_center(const std::vector<double>& vertices) const {
  double* center = new double[2] { 0., 0. };
  center[0] = (vertices[0] + vertices[3]) / 2;
  center[1] = (vertices[1] + vertices[3 + 1]) / 2;
  return center;
}

const std::vector<size_t>& FVM::Face::get_global_indices() const {
  return global_indices;
}

FVM::Face::~Face() {
  delete[] center;
}