#include "hyperplane.hh"
#include "projective_space.hh"

#include <memory>
#include <vector>


/* ************************************************************************* */
auto Hyperplane::contains(const ProjectivePoint& p) const -> bool
{
  auto field_zero = p.get_coordinate(0).get_field()->get_element(0);
  FieldElement sum = field_zero;

  for (size_t i = 0; i < linear_form_.size(); i++)
    sum = sum + (linear_form_[i] * p.get_coordinate(i));

  return sum == field_zero;
}


/* ************************************************************************* */
auto get_all_hyperplanes(shared_ptr<const ProjectiveSpace> ps) 
  -> vector<Hyperplane>
{
  auto ps_points = ps->get_all_points();
  vector<Hyperplane> hyperplanes;

  for (auto p : ps_points)
    hyperplanes.push_back(Hyperplane(p.get_coordinates()));

  return hyperplanes;
}


/* ************************************************************************* */
auto get_all_hyperplanes(const ProjectiveSpace& ps) -> vector<Hyperplane>
{
  auto ps_points = ps.get_all_points();
  vector<Hyperplane> hyperplanes;

  for (auto p : ps_points)
    hyperplanes.push_back(Hyperplane(p.get_coordinates()));

  return hyperplanes;
}
