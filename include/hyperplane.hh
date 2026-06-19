#pragma once

#include "field.hh"
#include "field_element.hh"
#include "projective_point.hh"
#include "projective_space.hh"

#include <cstddef>
#include <vector>

using namespace std;

class Hyperplane 
{
public:

  Hyperplane(vector<FieldElement> u) : linear_form_(u) {};

  ~Hyperplane() = default;

  auto contains(const ProjectivePoint& p) const -> bool
  {
    auto field_zero = p.get_coordinate(0).get_field()->get_element(0);
    FieldElement sum = field_zero;
    for (size_t i = 0; i < linear_form_.size(); i++)
      sum = sum + (linear_form_[i] * p.get_coordinate(i));
    return sum == field_zero;
  }

private:
  vector<FieldElement> linear_form_;
};

auto get_all_hyperplanes(const ProjectiveSpace& ps) -> vector<Hyperplane>;
