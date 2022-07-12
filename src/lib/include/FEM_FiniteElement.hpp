#ifndef __FINITE_ELEMENT__
#define __FINITE_ELEMENT__

#include <vector>
#include <array>
#include <functional>
#include <cmath>
#include <memory>

class IFiniteElement {
protected:
  std::vector<size_t> global_indices;
public:
  IFiniteElement() = delete;
  IFiniteElement(const std::vector<double> &vertices, const std::vector<size_t> global_indices);
  // Returns mass matrix element with local indices [i,j]
  virtual double get_mass(size_t i, size_t j) const = 0; 
  // Returns stiffness matrix elements with local indices [i,j]
  virtual double get_stiffness(size_t i, size_t j) const = 0;
  // Returns [i,i] element of lumped mass matrix
  virtual double get_lumped(size_t i) const = 0;
  // Returns point in parametric space [x,y,z]
  virtual void phys_to_param(const double* phys_in, double* param_out) const = 0;
  //virtual double* phys_to_param(double* point) const = 0;
  // Returns point in
  virtual void param_to_phys(const double* param_in, double* phys_out) const = 0;
  // Center coordinates
  virtual const double* get_center_coordinates() const = 0;
  // Get volume of the element
  virtual double get_volume() const = 0;
  // Get number of basis functions 
  virtual size_t get_number_basis_func() const = 0;
  // Get type of element
  virtual size_t get_element_type() const = 0;
  // Get global indices
  virtual const std::vector<size_t>& get_global_indices() const = 0;
  virtual ~IFiniteElement();
};


class LinearLineElement : public IFiniteElement {
protected:
  const double det_j;
  const std::array<double, 2> lumped_matrix;
  const std::array<double, 4> mass_matrix;
  const std::array<double, 4> stiffness_matrix;
  const std::array<std::function<double(const double*)>, 2> basis_functions;
  const double* base_point;
  const std::array<double, 1> center_point;
  std::array<double, 1> calc_center_point(const std::vector<double>& vertices) const;
public:
  double get_mass(size_t i, size_t j) const override; 
  double get_stiffness(size_t i, size_t j) const override;
  double get_lumped(size_t i) const override;
  void phys_to_param(const double* phys_in, double* param_out) const override;
  void param_to_phys(const double* param_in, double* phys_out) const override;
  const double* get_center_coordinates() const override;
  double get_volume() const override;
  size_t get_number_basis_func() const override;
  size_t get_element_type() const override;
  const std::vector<size_t>& get_global_indices() const override;
  LinearLineElement(const std::vector<double> &vertices, const std::vector<size_t> global_indices);
  ~LinearLineElement();
};

class LinearTriangleElement : public IFiniteElement {
protected:
  const std::array<double, 4> jacobi_matrix;
  const double det_j;
  const std::array<double, 3> lumped_matrix;
  const std::array<double, 9> mass_matrix;
  const std::array<double, 9> stiffness_matrix;
  const std::array<std::function<double(const double*)>, 3> basis_functions;
  const double* base_point;
  const std::array<double, 2> center_point;
  std::array<double, 2> calc_center_point(const std::vector<double>& vertices) const;
public:
  double get_mass(size_t i, size_t j) const override; 
  double get_stiffness(size_t i, size_t j) const override;
  double get_lumped(size_t i) const override;
  void phys_to_param(const double* phys_in, double* param_out) const override;
  void param_to_phys(const double* param_in, double* phys_out) const override;
  const double* get_center_coordinates() const override;
  double get_volume() const override;
  size_t get_number_basis_func() const override;
  size_t get_element_type() const override;
  const std::vector<size_t>& get_global_indices() const override;
  LinearTriangleElement(const std::vector<double> &vertices, const std::vector<size_t> global_indices);
  ~LinearTriangleElement();
};

class LinearRectangleElement : public IFiniteElement {
protected:
  const double det_j;
  const std::array<double, 4> lumped_matrix;
  const std::array<double, 16> mass_matrix;
  const std::array<double, 16> stiffness_matrix;
  const std::array<std::function<double(const double*)>, 4> basis_functions;
  const double* base_point;
  const std::array<double, 2> center_point;
  std::array<double, 2> calc_center_point(const std::vector<double>& vertices) const;
public:
  double get_mass(size_t i, size_t j) const override; 
  double get_stiffness(size_t i, size_t j) const override;
  double get_lumped(size_t i) const override;
  void phys_to_param(const double* phys_in, double* param_out) const override;
  void param_to_phys(const double* param_in, double* phys_out) const override;
  const double* get_center_coordinates() const override;
  double get_volume() const override;
  size_t get_number_basis_func() const override;
  size_t get_element_type() const override;
  const std::vector<size_t>& get_global_indices() const override;
  LinearRectangleElement(const std::vector<double> &vertices, const std::vector<size_t> global_indices);
  ~LinearRectangleElement();
};

class FiniteElementFactory {
protected:
public:
  static std::shared_ptr<IFiniteElement> create_element(
    const std::vector<double>& vertices,
    const std::vector<size_t>& global_indices,
    size_t element_type
  );
};

#endif