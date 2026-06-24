#include "field_vector.hh"
#include "field_element.hh"
#include <cassert>
#include <cstddef>
#include <vector>


/* ************************************************************************* */
auto normalize(const FieldVector& v) -> FieldVector
{
  FieldVector normalized_v;
  for (auto& x : v) 
  {
    if (!x.is_zero()) 
    {
      FieldElement inv_x = mult_inverse(x);
      for (auto& y : v)
        normalized_v.push_back(y * inv_x);
      break;
    } // end if
  } // end for
  return normalized_v;
}


/* ************************************************************************* */
auto weight(const FieldVector &v) -> int
{
  int w = 0;

  for (const auto& x : v)
    if (!x.is_zero())
      w++;

  return w;
}


/* ************************************************************************* */
auto hamming_distance(const FieldVector &left, const FieldVector &right) -> int
{
  assert(left.size() == right.size());

  int d = 0;

  for (size_t i = 0; i < left.size(); i++)
    if (left[i] != right[i])
      d++;

  return d;
}


/* ************************************************************************* */
auto get_unit_vector(const int& n, const int& i, 
    shared_ptr<const Field>& field) -> FieldVector
{
  assert(i < n);

  auto field_zero = FieldElement(field);
  auto field_one = FieldElement(field);

  field_zero.set_zero();
  field_one.set_one();

  auto unit_i = vector<FieldElement>(n, field_zero);
  unit_i[i] = field_one;

  return unit_i;
}


/* ************************************************************************* */
auto is_zero(const FieldVector & v) -> bool
{
  if (v.empty())
    return true;

  auto field_zero = v[0].get_field()->get_element(0);

  for (auto & x : v)
    if (x != field_zero) 
      return false;
  
  return true;
}


/* ************************************************************************* */
auto is_canonical(const FieldVector& v) -> bool
{
  if (v.empty())
    return false;

  auto field_zero = v[0].get_field()->get_element(0);
  auto field_one = v[0].get_field()->get_element(1);

  for (auto & x : v)
    if (x != field_zero)
      return x == field_one;

  return false;
}
