#include "FEM_MeshBuilder.hpp"

FiniteElementsMesh::FiniteElementsMesh(
  const Mesh&& mesh, 
  const std::vector<std::shared_ptr<IFiniteElement>>&& elements
) : mesh(mesh), elements(elements)
{ }
const std::shared_ptr<IFiniteElement>& FiniteElementsMesh::get_element(size_t i) const{
  return elements[i];
}
const Mesh* FiniteElementsMesh::get_mesh() const{
  return &mesh;
}

FiniteElementsMesh FiniteElementsMeshBuilder::BuildFromFile(const std::string& filename){
  std::shared_ptr<IFEMParser> parser = IFEMParser::Factory(filename);
  Mesh mesh = parser->load_mesh();
  return BuildFromMesh(mesh);
}

FiniteElementsMesh FiniteElementsMeshBuilder::BuildFromMesh(const Mesh& mesh){
  std::vector<std::shared_ptr<IFiniteElement>> elements;
  size_t N = mesh.get_elements_number();

  std::vector<size_t> cell_types = mesh.get_cell_types();
  std::vector<std::vector<size_t>> cells = mesh.get_cells();

  FiniteElementFactory elements_factory;

  for (size_t i = 0; i < N; i++) {
    if (cell_types[i] == 3) {
      std::vector<double> t_vert {
        mesh.get_vertex(cells[i][0])[0], mesh.get_vertex(cells[i][0])[1], mesh.get_vertex(cells[i][0])[2],
        mesh.get_vertex(cells[i][1])[0], mesh.get_vertex(cells[i][1])[1], mesh.get_vertex(cells[i][1])[2]
			};
      elements.push_back(
        elements_factory.create_element(
          t_vert,
          cells[i],
          1
        )
      );
    }
    if (cell_types[i] == 5) {
      std::vector<double> t_vert {
        mesh.get_vertex(cells[i][0])[0], mesh.get_vertex(cells[i][0])[1], mesh.get_vertex(cells[i][0])[2],
        mesh.get_vertex(cells[i][1])[0], mesh.get_vertex(cells[i][1])[1], mesh.get_vertex(cells[i][1])[2],
        mesh.get_vertex(cells[i][2])[0], mesh.get_vertex(cells[i][2])[1], mesh.get_vertex(cells[i][2])[2]
			};
      elements.push_back(
        elements_factory.create_element(
          t_vert,
          cells[i],
          2
        )
      );
    }
    if (cell_types[i] == 8)
    {
      std::vector<double> t_vert {
        mesh.get_vertex(cells[i][0])[0], mesh.get_vertex(cells[i][0])[1], mesh.get_vertex(cells[i][0])[2],
        mesh.get_vertex(cells[i][1])[0], mesh.get_vertex(cells[i][1])[1], mesh.get_vertex(cells[i][1])[2],
        mesh.get_vertex(cells[i][2])[0], mesh.get_vertex(cells[i][2])[1], mesh.get_vertex(cells[i][2])[2],
        mesh.get_vertex(cells[i][3])[0], mesh.get_vertex(cells[i][3])[1], mesh.get_vertex(cells[i][3])[2]
			};
      elements.push_back(
        elements_factory.create_element(
          t_vert,
          cells[i],
          3
        )
      );
    }
  }

  FiniteElementsMesh femesh(std::move(mesh), std::move(elements));
  return femesh;
}


FiniteVolumesMesh::FiniteVolumesMesh(
  const Mesh&& mesh, 
  const std::vector<std::shared_ptr<FiniteVolumeElement>>&& elements
) : mesh(mesh), elements(elements)
{ }
const std::shared_ptr<FiniteVolumeElement>& FiniteVolumesMesh::get_element(size_t i) const{
  return elements[i];
}
const Mesh* FiniteVolumesMesh::get_mesh() const{
  return &mesh;
}

FiniteVolumesMesh FiniteVolumesMeshBuilder::BuildFromFile(const std::string& filename) {
  std::shared_ptr<IFEMParser> parser = IFEMParser::Factory(filename);
  Mesh mesh = parser->load_mesh();
  return BuildFromMesh(mesh);
}

FiniteVolumesMesh FiniteVolumesMeshBuilder::BuildFromMesh(const Mesh& mesh) {
  std::vector<std::shared_ptr<FiniteVolumeElement>> elements;
  size_t N = mesh.get_elements_number();

  std::vector<size_t> cell_types = mesh.get_cell_types();
  std::vector<std::vector<size_t>> cells = mesh.get_cells();

  for (size_t i = 0; i < N; i++)
  {
    std::vector<double> t_vert;
    for (size_t vert_ind : cells[i])
    {
      t_vert.push_back(mesh.get_vertex(vert_ind)[0]);
      t_vert.push_back(mesh.get_vertex(vert_ind)[1]);
      t_vert.push_back(mesh.get_vertex(vert_ind)[2]);
    }
    elements.push_back(
      std::make_shared<FiniteVolumeElement>(t_vert, cells[i])
    );
  }
  FiniteVolumesMesh femesh(std::move(mesh), std::move(elements));
  return femesh;
}