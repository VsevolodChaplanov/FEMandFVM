#ifndef __FEM__BUILDER__
#define __FEM__BUILDER__

#include <string>
#include <vector>
#include <memory>
#include "FEM_Mesh.hpp"
#include "FEM_FiniteElement.hpp"
#include "FEM_Parser.hpp"
#include "FVM_FiniteVolumeElement.hpp"

// -------------------------------- FEM -------------------------------- //
class FiniteElementsMesh {
protected:
  const Mesh mesh;
  const std::vector<std::shared_ptr<IFiniteElement>> elements;
public:
  FiniteElementsMesh(const Mesh&& mesh, const std::vector<std::shared_ptr<IFiniteElement>>&& elements);
  const std::shared_ptr<IFiniteElement>& get_element(size_t i) const;
  const Mesh* get_mesh() const;
};

class FiniteElementsMeshBuilder {
protected:
public:
  static FiniteElementsMesh BuildFromFile(const std::string& filename);
  // Пока без указания границ
  static FiniteElementsMesh BuildFromMesh(const Mesh& mesh);
};

// -------------------------------- FEM -------------------------------- //
class FiniteVolumesMesh {
  protected:
  const Mesh mesh;
  const std::vector<std::shared_ptr<FiniteVolumeElement>> elements;
public:
  FiniteVolumesMesh(const Mesh&& mesh, const std::vector<std::shared_ptr<FiniteVolumeElement>>&& elements);
  const std::shared_ptr<FiniteVolumeElement>& get_element(size_t i) const;
  const Mesh* get_mesh() const;
};

class FiniteVolumesMeshBuilder {
protected:
public:
  static FiniteVolumesMesh BuildFromFile(const std::string& filename);
  // Пока без указания границ
  static FiniteVolumesMesh BuildFromMesh(const Mesh& mesh);
};


#endif