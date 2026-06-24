#include "projective_point.hh"
#include "field_vector.hh"


/* ************************************************************************* */
ProjectivePoint::ProjectivePoint(const FieldVector& coordinates) : 
    dimension_(coordinates.size()), coordinates_(normalize(coordinates)) 
{}


/* ************************************************************************* */
auto ProjectivePoint::concatenate(const FieldElement& fe) const 
  -> ProjectivePoint
{
  auto coordinates = coordinates_;
  coordinates.push_back(fe);
  return ProjectivePoint(coordinates);
}


/* ************************************************************************* */
auto operator==(const ProjectivePoint& left, const ProjectivePoint& right) 
  -> bool
{
  return left.dimension_ == right.dimension_
    && normalize(left.coordinates_) == normalize(right.coordinates_);
}


/* ************************************************************************* */
ostream& operator<<(ostream& output, const ProjectivePoint& right)
{
  output << "( ";

  for (auto& c : right.coordinates_)
    cout << c.index() << " ";

  output << ")";

  return output;
}
