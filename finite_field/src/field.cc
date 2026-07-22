#include "field.hh"
#include "field_element.hh"
#include "typedefs.hh"

#include <filesystem>
#include <fstream>

/* ************************************************************************* */
auto Field::get_element(uint64_t index) const -> FieldElement {
  FieldElement result(shared_from_this());
  uint64_t p = get_characteristic();

  for (int i = 0; i < get_dimension(); ++i) {
    result[i] = index % p;
    index /= p;
  }

  return result;
}

/* ************************************************************************* */
auto Field::get_all_elements() const -> std::vector<FieldElement> {
  std::vector<FieldElement> elems;

  for (uint64_t i = 0; i < get_order(); ++i)
    elems.push_back(get_element(i));

  return elems;
}

/* ************************************************************************* */
auto Field::parse_field_file(const filesystem::path &file) -> Field {
  ifstream in(file);

  if (!in)
    throw runtime_error("Cannot open field file: " + file.string());

  string line;

  if (!getline(in, line))
    throw runtime_error("Missing characteristic line.");

  const string characteristic_prefix = "characteristic=";

  if (!line.starts_with(characteristic_prefix))
    throw runtime_error("Expected characteristic line.");

  int characteristic = stoi(line.substr(characteristic_prefix.size()));

  if (!getline(in, line))
    throw runtime_error("Missing reduction polynomial line.");

  const string polynomial_prefix = "reduction_polynomial=";

  if (!line.starts_with(polynomial_prefix))
    throw runtime_error("Expected reduction polynomial line.");

  vector<Fint> coefficients;

  istringstream iss(line.substr(polynomial_prefix.size()));

  Fint coeff;

  while (iss >> coeff)
    coefficients.push_back(coeff);

  if (coefficients.empty())
    throw runtime_error("Empty reduction polynomial.");

  return Field(characteristic, coefficients.size(), coefficients);
}
