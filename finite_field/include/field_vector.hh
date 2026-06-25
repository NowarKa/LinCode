#pragma once

#include "field.hh"
#include "field_element.hh"

#include <memory>
#include <vector>

using namespace std;
using FieldVector = vector<FieldElement>;

/**
  Defines FieldVector and some various operations.

  FieldVector is a vector (in C++ meaning) of FieldElement.
 */

/**
  Returns a normalized copy of the given FieldVector.
 
  The normalized vector is obtained by multiplying the vector by the
  inverse of its first non-zero coefficient, making the first non-zero
  coefficient equal to 1.
 
  @param v The vector to normalize.
 */
auto normalize(const FieldVector& v) -> FieldVector;


/**
  Returns the weight of a given FieldVector.

  The weight of a vector is the number its non-zero coefficients.

  @param v The vector whose weight is returned.
 */
auto hamming_weight(const FieldVector &v) -> int;


/**
  Returns the hamming distance of two FieldVectors.

  The Hamming distance of two vectors is the weight of their difference.

  @param left The left vector.
  @param right The right vector.
 */
auto hamming_distance(const FieldVector& left, const FieldVector& right) 
  -> int;


/**
  Retruns the i-th unit vector of GF(q)^n.

  @param n The length of the unit vector.
  @param i The index of the unit vector.
  @param field The field.
 */
auto get_unit_vector(const int& n, const int& i, 
    shared_ptr<const Field>& field) -> FieldVector;


/**
  Returns true if the given vector is equal to the zero vector and false 
  otherwise.

  @param v The vector to be compared to zero.
 */
auto is_zero(const FieldVector & v) -> bool;


auto is_canonical(const FieldVector& v) -> bool;
