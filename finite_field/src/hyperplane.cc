#include "hyperplane.hh"
#include "projective_space.hh"

#include <memory>
#include <unordered_map>
#include <vector>

/* ************************************************************************* */
auto Hyperplane::contains(const ProjectivePoint &p) const -> bool {
  auto field_zero = p.get_coordinate(0).get_field()->get_element(0);
  FieldElement sum = field_zero;

  for (size_t i = 0; i < linear_form_.size(); i++)
    sum = sum + (linear_form_[i] * p.get_coordinate(i));

  return sum == field_zero;
}

/* ************************************************************************* */
auto get_all_hyperplanes(shared_ptr<const ProjectiveSpace> ps)
    -> vector<Hyperplane> {
  auto ps_points = ps->get_all_points();
  vector<Hyperplane> hyperplanes;

  for (auto p : ps_points)
    hyperplanes.push_back(Hyperplane(p.get_coordinates()));

  return hyperplanes;
}

/* ************************************************************************* */
auto Hyperplane::get_all_hyperplanes(const ProjectiveSpace &ps)
    -> vector<Hyperplane> {
  auto ps_points = ps.get_all_points();
  vector<Hyperplane> hyperplanes;

  for (auto p : ps_points)
    hyperplanes.push_back(Hyperplane(p.get_coordinates()));

  return hyperplanes;
}

/* ************************************************************************* */
auto Hyperplane::containing_map(shared_ptr<const ProjectiveSpace> ps_kp1,
                                shared_ptr<const vector<Hyperplane>> h_kp1)
    -> unordered_map<Hyperplane, vector<int>> {
  auto points = ps_kp1->get_all_points();
  unordered_map<Hyperplane, vector<int>> containing_map;

  for (const auto &h : *h_kp1)
    for (int i = 0; i < points.size(); i++)
      if (h.contains(points[i]))
        containing_map[h].push_back(i);

  return containing_map;
}

/* ************************************************************************* */
bool operator==(const Hyperplane &left, const Hyperplane &right) {
  return left.linear_form_ == right.linear_form_;
}
