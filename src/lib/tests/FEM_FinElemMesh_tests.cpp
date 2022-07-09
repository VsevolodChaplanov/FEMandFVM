#include "TestsHeader.hpp"

template<typename T>
bool compare_vectors(const std::vector<T>& first, const std::vector<T>& second) {
  if (first.size() != second.size()) {
    return false;
  }
  for (size_t i = 0; i < first.size(); i++) {
    if (first[i] != second[i]) {
      return false;
    }
  }
  return true;
}

TEST_CASE( "Finite elements mesh parser", "[ParserVTK][Parser]" )
{
  FEM::VtkFEMParser test_vtk("./test_resources/Rect.1.1mesh.1.vtk");

  FEM::Mesh test_mesh = test_vtk.load_mesh();

  REQUIRE( test_mesh.get_vertices_number() == 98 );
  REQUIRE( test_mesh.get_elements_number() == 198 );
  REQUIRE( test_mesh.get_cell_types().size() == 198);

  std::vector<double> vertices_test = test_mesh.get_vertices();
  std::vector<std::vector<size_t>> cells_test = test_mesh.get_cells();
  std::vector<size_t> cell_types_test = test_mesh.get_cell_types();

  // -------------- Vertices data check -------------- // 
  // First vertice
  for (size_t i = 0; i < 3; i++)
  {
    CHECK( vertices_test[i] == 0.0 );
  }

  // Random vertice from the centre
  CHECK( vertices_test[177] ==  0.31331768157997741);
  CHECK( vertices_test[178] ==  0.75237107335505982);
  CHECK( vertices_test[179] ==  0.0);

  // Last vertice
  CHECK( vertices_test[98 * 3 - 1] == 0.0);
  CHECK( vertices_test[98 * 3 - 2] == 0.8253550026152658);
  CHECK( vertices_test[98 * 3 - 3] == 0.7489644266062752);
  // -------------- Vertices data check -------------- // 

  // -------------- Cells data check -------------- // 
  // Nodes of the rectangle
  REQUIRE( cells_test[0].size() == 1);
  CHECK( cells_test[0][0] == 0 );
  // Edges of the rectangle
  REQUIRE( cells_test[4].size() == 2);
  CHECK( cells_test[4][0] == 0 );
  CHECK( cells_test[4][1] == 4 );
  // Internal triangles
  REQUIRE( cells_test[37].size() == 3);
  CHECK( cells_test[36][0] == 67);
  CHECK( cells_test[36][1] == 78);
  CHECK( cells_test[36][2] == 37);
  // Last triangle
  REQUIRE( cells_test[197].size() == 3);
  CHECK( cells_test[97][0] == 6);
  CHECK( cells_test[97][1] == 36);
  CHECK( cells_test[97][2] == 5);
  // Last triangle
  REQUIRE( cells_test[197].size() == 3);
  CHECK( cells_test[197][0] == 82);
  CHECK( cells_test[197][1] == 97);
  CHECK( cells_test[197][2] == 60);
  // -------------- Cells data check -------------- // 

  // -------------- Cell types data check -------------- // 

  // All cell types check
  for (std::vector<size_t>::iterator it = cell_types_test.begin(); it != cell_types_test.begin() + 4; it++)
  {
    CHECK( *it == 1 );
  }
  for (std::vector<size_t>::iterator it = cell_types_test.begin() + 4; it != cell_types_test.begin() + 36; it++)
  {
    CHECK( *it == 3 );
  }
  for (std::vector<size_t>::iterator it = cell_types_test.begin() + 36; it != cell_types_test.end(); it++)
  {
    CHECK( *it == 5 );
  }
  // -------------- Cell types data check -------------- // 
}

TEST_CASE( "FEM linear line element assemble, test properties", "[FEM][LinearLineElement]" ) {
  std::unique_ptr<FEM::IFiniteElement> test_linear_element = FEM::FiniteElementFactory::create_element(
    {0., 0., 0., 1., 0., 0.},
    {0, 1},
    1
  );

  double* test_center = new double[3] {.5, 0., 0.};
  std::vector<double> test_lumped {0.5, 0.5};
  std::vector<double> test_mass {1. / 3., 1. / 6., 1. / 6., 1. / 3.};
  std::vector<double> test_stiffness {1., -1., -1., 1.};


  CHECK( test_linear_element->get_center_coordinates()[0] == test_center[0] );

  CHECK( test_linear_element->get_element_type() == 1 );

  CHECK( test_linear_element->get_number_basis_func() == 2);

  CHECK( test_linear_element->get_volume() == 1.);

  CHECK( test_linear_element->get_lumped(0) == test_lumped[0] );
  CHECK( test_linear_element->get_lumped(1) == test_lumped[1] );

  CHECK( test_linear_element->get_mass(0, 0) == test_mass[0] );
  CHECK( test_linear_element->get_mass(0, 1) == test_mass[1] );
  CHECK( test_linear_element->get_mass(1, 0) == test_mass[2] );
  CHECK( test_linear_element->get_mass(1, 1) == test_mass[3] );

  CHECK( test_linear_element->get_stiffness(0, 0) == test_stiffness[0] );
  CHECK( test_linear_element->get_stiffness(0, 1) == test_stiffness[1] );
  CHECK( test_linear_element->get_stiffness(1, 0) == test_stiffness[2] );
  CHECK( test_linear_element->get_stiffness(1, 1) == test_stiffness[3] );
}

TEST_CASE( "FEM linear triange element assemble, test properties", "[FEM][LinearTriangleElement]" ) {
  std::unique_ptr<FEM::IFiniteElement> test_linear_element = FEM::FiniteElementFactory::create_element(
      {0., 0., 0., 1., 0., 0., 0., 1., 0.},
      {0, 1, 2},
      2
  );


  double* test_center = new double[3] {1. / 3., 1. / 3., 0.};
  std::vector<double> test_lumped {1. / 6., 1. / 6., 1. / 6.};
  std::vector<double> test_mass {1./12.,	1./24.,	1./24.,
    1./24.,	1./12.,	1./24.,
    1./24.,	1./24.,	1./12.
  };
  std::vector<double> test_stiffness {1.,	-(1./2.),	-(1./2.),
    -(1./2.),	1./2.,	0.,
    -(1./2.),	0.,	1./2.
  };


  CHECK( test_linear_element->get_center_coordinates()[0] == test_center[0] );
  CHECK( test_linear_element->get_center_coordinates()[1] == test_center[1] );

  CHECK( test_linear_element->get_element_type() == 2 );

  CHECK( test_linear_element->get_number_basis_func() == 3 );

  CHECK( test_linear_element->get_volume() == 0.5 );

  CHECK( test_linear_element->get_lumped(0) == test_lumped[0] );
  CHECK( test_linear_element->get_lumped(1) == test_lumped[1] );
  CHECK( test_linear_element->get_lumped(2) == test_lumped[2] );
  
  CHECK( test_linear_element->get_mass(0, 0) == test_mass[0] );
  CHECK( test_linear_element->get_mass(0, 1) == test_mass[1] );
  CHECK( test_linear_element->get_mass(0, 2) == test_mass[2] );
  CHECK( test_linear_element->get_mass(1, 0) == test_mass[3] );
  CHECK( test_linear_element->get_mass(1, 1) == test_mass[4] );
  CHECK( test_linear_element->get_mass(1, 2) == test_mass[5] );
  CHECK( test_linear_element->get_mass(2, 0) == test_mass[6] );
  CHECK( test_linear_element->get_mass(2, 1) == test_mass[7] );
  CHECK( test_linear_element->get_mass(2, 2) == test_mass[8] );

  CHECK( test_linear_element->get_stiffness(0, 0) == test_stiffness[0] );
  CHECK( test_linear_element->get_stiffness(0, 1) == test_stiffness[1] );
  CHECK( test_linear_element->get_stiffness(0, 2) == test_stiffness[2] );
  CHECK( test_linear_element->get_stiffness(1, 0) == test_stiffness[3] );
  CHECK( test_linear_element->get_stiffness(1, 1) == test_stiffness[4] );
  CHECK( test_linear_element->get_stiffness(1, 2) == test_stiffness[5] );
  CHECK( test_linear_element->get_stiffness(2, 0) == test_stiffness[6] );
  CHECK( test_linear_element->get_stiffness(2, 1) == test_stiffness[7] );
  CHECK( test_linear_element->get_stiffness(2, 2) == test_stiffness[8] );
}

TEST_CASE( "Save as vtk format tests", "[Savevtk]" )
{
  FEM::VtkFEMParser vtk_gmsh("./test_resources/Rect.1.1mesh.1.vtk");
  FEM::Mesh mesh_gmsh_vtk = vtk_gmsh.load_mesh();

  mesh_gmsh_vtk.save_mesh_vtk("./test_resources/Rect.1.1mesh.1_test_save.vtk");

  FEM::VtkFEMParser vtk_mymesh("./test_resources/Rect.1.1mesh.1_test_save.vtk");
  FEM::Mesh mesh_mymesh_vtk = vtk_mymesh.load_mesh();

  std::vector<size_t> element_types_gmsh = mesh_gmsh_vtk.get_cell_types();
  std::vector<size_t> element_types_mymesh = mesh_mymesh_vtk.get_cell_types();

  REQUIRE( element_types_gmsh.size() == element_types_mymesh.size() );

  for (size_t i = 0; i < element_types_gmsh.size(); i+= 10) {
    CHECK( element_types_gmsh[i] == element_types_mymesh[i] );
  }

  std::vector<std::vector<size_t>> element_vert_gmsh = mesh_gmsh_vtk.get_cells();
  std::vector<std::vector<size_t>> element_vert_mymesh = mesh_mymesh_vtk.get_cells();

  REQUIRE( element_types_gmsh.size() == element_types_mymesh.size() );

  for (size_t i = 0; i < element_vert_gmsh.size(); i+= 10) {
    CHECK( 
      compare_vectors(
        element_vert_gmsh[i],
        element_vert_mymesh[i]
      )
    );
  }

  CHECK( Approx(mesh_gmsh_vtk.get_vertex(0)[0]) == mesh_mymesh_vtk.get_vertex(0)[0] );
  CHECK( Approx(mesh_gmsh_vtk.get_vertex(0)[1]) == mesh_mymesh_vtk.get_vertex(0)[1] );
  CHECK( Approx(mesh_gmsh_vtk.get_vertex(0)[2]) == mesh_mymesh_vtk.get_vertex(0)[2] );

  CHECK( Approx(mesh_gmsh_vtk.get_vertex(10)[0]) == mesh_mymesh_vtk.get_vertex(10)[0] );
  CHECK( Approx(mesh_gmsh_vtk.get_vertex(10)[1]) == mesh_mymesh_vtk.get_vertex(10)[1] );
  CHECK( Approx(mesh_gmsh_vtk.get_vertex(10)[2]) == mesh_mymesh_vtk.get_vertex(10)[2] );
}

TEST_CASE( "FVM create a finite volume cell", "[FVM]" ) {
  FVM::FiniteVolumeElement test_element(
    {0., 0., 0., 1., 0., 0., 0., 1., 0.},
    {0, 1, 2}
  );

  CHECK(
    compare_vectors(
      test_element.get_global_indices(),
      {0, 1, 2}
    )
  );

  CHECK(test_element.cell_center()[0] == 1./3.);
  CHECK(test_element.cell_center()[1] == 1./3.);
}
