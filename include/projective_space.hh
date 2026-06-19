#pragma once

#include "field.hh"
#include "projective_point.hh"

#include <cstddef>
#include <cstdint>
#include <memory>
#include <vector>

using namespace std;

/**
  Defines a projective space of characteristic q and dimension k (PG(k-1,q)). 
  It stores the following (private) data:

  @param dimension_ the dimension of ProjectiveSpace.
  @param field_ a shared pointer to the field GF(q).
*/
class ProjectiveSpace 
{
public:
  /**
    Constructs a projective space.
    @param k ProjectiveSpace's dimension.
    @param field a shared pointer to ProjectiveSpace's field.
   */
  ProjectiveSpace(uint32_t k, shared_ptr<const Field> field)
    : dimension_(k-1), field_(field) {}


  /**
    Destructs a ProjectiveSpace.
   */
  ~ProjectiveSpace() = default;


  /**
    Returns all points in this ProjectiveSpace.
 
    The number of points is (q^k - 1) / (q - 1).
   */
  auto get_all_points() const -> vector<ProjectivePoint>;


  /**
    Returns the index of the i-th unit vector in the vector
    returned by get_all_points().
   */
  auto get_index_unit_vector(size_t i) const -> size_t;

private:
  uint32_t dimension_;
  shared_ptr<const Field> field_;


  /**
    Function used in order to generate all points in get_all_point().
   */
  auto generate_all_points(FieldVector & v, size_t pos, 
      vector<ProjectivePoint>& generated_points) const -> void;
};
