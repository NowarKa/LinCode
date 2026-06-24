#pragma once

#include "field.hh"
#include "field_element.hh"
#include "field_vector.hh"

#include <cstddef>
#include <cstdint>
#include <vector>

using namespace std;
using FieldVector = vector<FieldElement>;

/**
  Defines a projective point. Also defines some operation on projective points.

  In mathematics, a projective point is a one dimension subspace of the linear 
  space GF(q)^k.

  It stores the following (private) data:
  @param coordinates_ the coordinates of a generating vector (normalized).
  @param dimension_ the length of the generating vector (ie the k such that 
  ProjectivePoint is a subspace of GF(q)^k).
*/
class ProjectivePoint 
{
public:
  /**
    Constructs a projective point without initializing its coordinates.
    @param k ProjectivePoint's dimension.
   */
  ProjectivePoint(uint32_t k) : dimension_(k) {}


  /**
    Constructs a projective point, initializes its coordinates and normalize them.
    @param coordinates a vector of FieldElement.
    @param k ProjectivePoint's dimension (length of coordinates = dimension).
   */
  ProjectivePoint(const FieldVector& coordinates);


  /**
    Destructs a ProjectivePoint.
   */
  ~ProjectivePoint() = default;


  /**
    Returns the coordinates of the generating vector.
   */
  auto get_coordinates(void) const -> const FieldVector&;
  

  /**
    Returns the i-th coordinate of the generating vector.
   */
  auto get_coordinate(size_t i) const -> FieldElement;
  

  /**
    Sets the i-th coefficient of ProjectivePoint's coordinates to a specified 
    constant.
    @param i input the index to the i-th coordinate.
    @param coordinate the value of the i-th coordinate.
   */
  auto set_coordinate(size_t i, const FieldElement& coordinate) -> void;


  /**
    Sets ProjectivePoint's coordinates to a specified vector.
    @param coordinates the value of the coordinates.
   */
  auto set_coordinates(const FieldVector& coordinates) -> void;


  /**
    Returns true if the given vector belongs to ProjectivePoint, and false 
    otherwise.
    @param v the specified vector.
   */
  auto contains(const FieldVector& v) const -> bool;


  /**
    Returns the ProjectivePoint resulted from the concatenation of 
    ProjectivePoint and a specified FieldElement.
    @param fe the specified FieldElement.
    @returns ProjectivePoint(<coordinates_,fe>)
   */
  auto concatenate(const FieldElement& fe) const -> ProjectivePoint;
  

  /**
    Compares two ProjectivePoints, checks for equality
    @param right ProjectivePoint on the right of == sign
    @param left ProjectivePoint on the left of == sign
    @returns true if equal, false otherwise
   */
  friend auto operator==(const ProjectivePoint& left, 
      const ProjectivePoint& right) -> bool;


  /**
    Outputs a ProjectivePoint to the standard output
    @param right ProjectivePoint on the right of << sign
    @param output an output stream on the left of << sign
    @returns output stream
    */
  friend ostream& operator<<(ostream& output, const ProjectivePoint& right);

private:
  // I am not sure whether it is important to add the integers k and q.
  FieldVector coordinates_;
  uint64_t dimension_;
};


template<> struct std::hash<ProjectivePoint> {
  std::size_t operator()(ProjectivePoint const& p) const noexcept 
  {
    size_t seed = 0;
    // std::size_t h_field = std::hash<uint64_t>{}(p.get_dimension());
    // boost::hash_combine(seed, h_field);

    // auto normalized_p = normalize(p.get_coordinates());
    for (const auto& c : p.get_coordinates())
      boost::hash_combine(seed, hash<FieldElement>{}(c));

    return seed;
  }
};


inline auto ProjectivePoint::get_coordinates(void) const -> const FieldVector&
{
  return coordinates_;
}


inline auto ProjectivePoint::get_coordinate(size_t i) const -> FieldElement
{
  // TODO: check whether i < coordinates_.size()
  return coordinates_[i];
}


inline auto ProjectivePoint::set_coordinate(size_t i, 
    const FieldElement& coordinate) -> void
{
  coordinates_[i] = coordinate;
}


inline auto ProjectivePoint::set_coordinates(const FieldVector& coordinates) 
   -> void
{
  coordinates_ = coordinates;
}


inline auto ProjectivePoint::contains(const FieldVector& v) const -> bool
{
  return normalize(v) == coordinates_;
}
