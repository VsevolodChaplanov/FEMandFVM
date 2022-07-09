#ifndef __FINITE_VOLUME__
#define __FINITE_VOLUME__

#include <vector>
#include <cmath>

namespace FVM
{
  class Face;
  // Пока только двумерный элемент
  class FiniteVolumeElement {
  protected:
    const std::vector<size_t> global_indices;
    std::vector<Face> faces;
    size_t n_vert;
    double area;
    const double* center;
    // Calculate volume of cell
    // Function saves sign of the area
    double calc_area(const std::vector<double>& vertices) const;
    const double* calc_cell_center(const std::vector<double>& vertices) const;
  public:
    FiniteVolumeElement(const std::vector<double>& vertices, 
      const std::vector<size_t>& global_indices);
    // Calculates barycenter
    const double* cell_center() const;
    const std::vector<size_t>& get_global_indices() const;
    ~FiniteVolumeElement();
  };
  // Для двумерной задачи фэйсы всегда ребра
  // Пусть индексы и ребра сразу упорядочены
  class Face {
  protected: 
    const std::vector<size_t> global_indices;
    size_t n_vert;
    double lenght;
    const double* center;
    double calc_lenght(const std::vector<double>& vertices) const;
    const double* calc_face_center(const std::vector<double>& vertices) const;
  public:
    Face(const std::vector<double>& vertices, 
      const std::vector<size_t>& global_indices);
    const std::vector<size_t>& get_global_indices() const;
    ~Face();
  };
} // namespace FVM




#endif