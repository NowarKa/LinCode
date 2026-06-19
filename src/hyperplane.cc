#include "hyperplane.hh"
#include "projective_space.hh"

#include <memory>
#include <vector>

auto get_all_hyperplanes(shared_ptr<const ProjectiveSpace> ps) -> vector<Hyperplane>
{
  auto ps_points = ps->get_all_points();
  vector<Hyperplane> hyperplanes;
  for (auto p : ps_points)
    hyperplanes.push_back(Hyperplane(p.get_coordinates()));
  return hyperplanes;
}

auto get_all_hyperplanes(const ProjectiveSpace& ps) -> vector<Hyperplane>
{
  auto ps_points = ps.get_all_points();
  vector<Hyperplane> hyperplanes;
  for (auto p : ps_points)
    hyperplanes.push_back(Hyperplane(p.get_coordinates()));
  return hyperplanes;
}

