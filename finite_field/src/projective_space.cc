#include "projective_space.hh"

#include <cassert>
#include <cstddef>

/* ************************************************************************* */
auto ProjectiveSpace::get_index_unit_vector(size_t i) const -> size_t {
  return (pow(field_->get_order(), i) - 1) / (field_->get_order() - 1);
}

/* ************************************************************************* */
auto ProjectiveSpace::get_all_points() const -> vector<ProjectivePoint> {
  return points_;
}

/* ************************************************************************* */
auto ProjectiveSpace::generate_all_points() const -> vector<ProjectivePoint> {
  assert(field_ != nullptr);
  vector<ProjectivePoint> result;
  FieldVector v(dimension_ + 1, FieldElement(field_));

  for (size_t pivot = 0; pivot <= dimension_; ++pivot) {
    fill(v.begin(), v.end(), FieldElement(field_));

    for (size_t i = 0; i < pivot; ++i)
      v.at(i) = field_->get_element(0);

    v[pivot] = field_->get_element(1);

    generate_tail(v, pivot + 1, result);
  }

  return result;
}

/* ************************************************************************* */
void ProjectiveSpace::generate_tail(FieldVector &v, size_t pos,
                                    vector<ProjectivePoint> &result) const {
  if (pos == v.size()) {
    result.emplace_back(v);
    return;
  }

  for (const auto &fe : field_->get_all_elements()) {
    v[pos] = fe;
    generate_tail(v, pos + 1, result);
  }
}
