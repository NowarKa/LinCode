#pragma once

#include "field.hh"
#include "field_element.hh"
#include "projective_point.hh"
#include "projective_space.hh"

#include <cstddef>
#include <memory>
#include <unordered_map>
#include <vector>

using namespace std;

/**
  Defines hyperplane in ProjectiveSpace and some operations.

  In mathematics, an hyperplane is a (k-1)-dimensional subspace of a
  k-dimensional vector space.

  It stores the following private data:
  @param linear_form_ The linear form that defines the subspace.
 */
class Hyperplane {
public:
  /**
    Constructs an Hyperplane and initializes its linear form.

    @param u The linear form.
   */
  Hyperplane(vector<FieldElement> u) : linear_form_(u) {};

  /**
    Destructs an Hyperplane.
   */
  ~Hyperplane() = default;

  /**
   * Returns linear_form_
   */
  auto get_linear_form() const -> vector<FieldElement>;

  /**
    Returns true if the given ProjectivePoint belongs to the Hyperplane, and
    false otherwise.

    @param p A specified ProjectivePoint.
   */
  auto contains(const ProjectivePoint &p) const -> bool;

  /**
    Returns the set of all Hyperplane of a given ProjectiveSpace.

    @param ps A specified ProjectiveSpace.
    */
  static auto get_all_hyperplanes(const ProjectiveSpace &ps)
      -> vector<Hyperplane>;

  static auto containing_map(shared_ptr<const ProjectiveSpace> ps_kp1,
                             shared_ptr<const vector<Hyperplane>> h_kp1)
      -> unordered_map<Hyperplane, vector<int>>;

  friend bool operator==(const Hyperplane &left, const Hyperplane &right);

private:
  vector<FieldElement> linear_form_;
};

inline auto Hyperplane::get_linear_form() const -> vector<FieldElement> {
  return linear_form_;
}

template <> struct std::hash<Hyperplane> {
  std::size_t operator()(Hyperplane const &H) const noexcept {
    size_t seed = 0;

    for (const auto &c : H.get_linear_form())
      boost::hash_combine(seed, hash<FieldElement>{}(c));

    return seed;
  }
};
