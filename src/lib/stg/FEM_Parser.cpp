#include "FEM_Parser.hpp"

std::shared_ptr<IFEMParser> IFEMParser::Factory(const std::string &filename)
{
  boost::char_separator<char> sep {"."};
  boost::tokenizer<boost::char_separator<char>> tokenizer(filename, sep);
  std::vector<std::string> filename_parsed;
  for (boost::tokenizer<boost::char_separator<char>>::iterator it = tokenizer.begin(); 
    it != tokenizer.end(); ++it)
  {
      filename_parsed.push_back(*it);
  }


  if ( *(filename_parsed.end() - 1) == "vtk")
  {
    return std::make_shared<VtkFEMParser>(filename);
  }
  throw std::runtime_error("This file format is not supported");
}

IFEMParser::IFEMParser(const std::string &filename) : filename(filename) { }

IFEMParser::~IFEMParser() { }

VtkFEMParser::VtkFEMParser(const std::string &filename) : IFEMParser(filename) { }

Mesh  VtkFEMParser::load_mesh()
{
  std::ifstream file(filename, std::ios_base::in);

  std::string line;

  boost::char_separator<char> sep {" "};

  // Continous storage of the vertices as:
  // 	{x0,y0,z0,x1,y1,z1, ... }
  size_t Nvert;
  size_t Nelem;
  std::vector<double> vertices;
  std::vector<std::vector<size_t>> cells;
  std::vector<size_t> cell_types;

  size_t k = 0;
  if (!file.is_open())
  {
    throw std::runtime_error("Can't open file " + filename);
  }
  
  while (getline(file, line))
  {
    boost::tokenizer<boost::char_separator<char>> tok(line, sep);

    if (line == "")
    {
      state = DataType::Null;
      continue;
    } else if (*tok.begin() == "POINTS")
    {
      state = DataType::POINTS;
      continue;
    } else if (*tok.begin() == "CELLS")
    {
      state = DataType::CELLS;
      continue;
    } else if (*tok.begin() == "CELL_TYPES")
    {
      state = DataType::CELL_TYPES;
      continue;
    }

    if (state == DataType::POINTS)
    {
      for (boost::tokenizer<boost::char_separator<char>>::iterator it = tok.begin(); it != tok.end(); ++it)
      {
        vertices.push_back(std::stod(*it));
      }
    } else if(state == DataType::CELLS)
    {
      std::vector<size_t> temp_v;
      for (boost::tokenizer<boost::char_separator<char>>::iterator it = ++tok.begin(); it != tok.end(); ++it)
      {
        temp_v.push_back(std::stol(*it));
      }
      cells.push_back(temp_v);
    } else if(state == DataType::CELL_TYPES)
    {
      cell_types.push_back(std::stol(*tok.begin()));
    }
    k++;
  }

  Nvert = vertices.size() / 3;
  Nelem = cells.size();

  file.close();

  Mesh mesh(
      Nvert,
      Nelem,
      vertices,
      cells,
      cell_types
  );

  return mesh;
}

VtkFEMParser::~VtkFEMParser() { }