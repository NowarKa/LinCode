#include "field_element.hh"
#include <memory>
#include "field.hh"

auto Field::get_element(uint64_t index) const -> FieldElement
{
  FieldElement result(shared_from_this());
  uint64_t p = get_characteristic();

  for (int i = 0; i < get_dimension(); ++i)
  {
    result[i] = index % p;
    index /= p;
  }

  return result;
}

auto Field::get_all_elements() const -> std::vector<FieldElement>
{
  std::vector<FieldElement> elems;

  for (uint64_t i = 0; i < get_order(); ++i)
    elems.push_back(get_element(i));

  return elems;
}
