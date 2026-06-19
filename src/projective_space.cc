#include "projective_space.hh"

#include <cstddef>


/* ************************************************************************* */
auto ProjectiveSpace::get_index_unit_vector(size_t i) const -> size_t
{
  return (pow(field_->get_order(),i) - 1)/(field_->get_order() - 1);
}


/* ************************************************************************* */
auto ProjectiveSpace::get_all_points() const -> vector<ProjectivePoint>
{
  vector<ProjectivePoint> result;
  FieldVector v(dimension_ + 1, FieldElement(field_));

  generate_all_points(v, 0, result);

  return result;
}


/* ************************************************************************* */
// First attempt to generate all points.
// It is clearly not the optimal way.
auto ProjectiveSpace::generate_all_points(FieldVector & v, size_t pos, 
    vector<ProjectivePoint>& generated_points) const -> void
{
  if (pos == dimension_ + 1) 
  {
    if (!is_zero(v) && is_canonical(v))
      generated_points.push_back(ProjectivePoint(v));
    return;
  }

  for (auto& fe : field_->get_all_elements()) 
  {
    v[pos] = fe;
    generate_all_points(v, pos + 1, generated_points);
  }
}
