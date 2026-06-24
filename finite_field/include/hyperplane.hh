#pragma once

#include "field.hh"
#include "field_element.hh"
#include "projective_point.hh"
#include "projective_space.hh"

#include <cstddef>
#include <vector>

using namespace std;

/**
  Defines hyperplane in ProjectiveSpace and some operations.

  In mathematics, an hyperplane is a (k-1)-dimensional subspace of a 
  k-dimensional vector space.

  It stores the following private data:
  @param linear_form_ The linear form that defines the subspace.
 */
class Hyperplane 
{
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
    Returns true if the given ProjectivePoint belongs to the Hyperplane, and 
    false otherwise.

    @param p A specified ProjectivePoint.
   */
  auto contains(const ProjectivePoint& p) const -> bool;

private:
  vector<FieldElement> linear_form_;
};


/**
  Returns the set of all Hyperplane of a given ProjectiveSpace.

  @param ps A specified ProjectiveSpace.
 */
auto get_all_hyperplanes(const ProjectiveSpace& ps) -> vector<Hyperplane>;
