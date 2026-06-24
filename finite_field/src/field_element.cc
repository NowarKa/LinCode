#include <iostream>
#include <iterator>
#include <vector>

#include "field_element.hh"
#include "field.hh"
#include "polynomial_arith.hh"

using namespace std;

/* ************************************************************************* */
FieldElement::FieldElement(const std::shared_ptr<const Field>& gf) : 
  field(gf), modulus(field->get_characteristic()), 
  dimension(field->get_dimension()), 
  polynomial(vector<Fint>(dimension))
{}


/* ************************************************************************* */
FieldElement::FieldElement(const std::shared_ptr<const Field>& gf, 
    const vector<Fint> coeff) :
  field(gf), modulus(field->get_characteristic()), 
  dimension(field->get_dimension()), 
  polynomial(dimension, 0)
{
  int n = std::min<int>(dimension, coeff.size());

  for (int i = 0; i < n; ++i)
    polynomial[i] = coeff[i];

  if (polynomial_degree(coeff) >= dimension)
    throw std::runtime_error("Polynomial degree too large");
}


/* ************************************************************************* */
auto FieldElement::index() const -> uint64_t
{
  uint64_t result = 0;
  uint64_t power  = 1;

  for (int i = 0; i < dimension; ++i)
  {
    result += polynomial[i] * power;
    power *= modulus;
  }

  return result;
}


/* ************************************************************************* */
auto FieldElement::set_zero(void) -> void
{
  for (int i=0; i<dimension; i++)
    polynomial[i] = 0;
}


/* ************************************************************************* */
auto FieldElement::set_one(void) -> void
{
  polynomial[0] = 1;
  for (int i = 1; i < dimension; i++)
    polynomial[i] = 0;
}


/* ************************************************************************* */
Fint& FieldElement::operator[](int index)
{
  return polynomial[index];
}


Fint FieldElement::operator[](int index) const
{
  return polynomial[index];
}


/* ************************************************************************* */
FieldElement FieldElement::operator-()
{
  FieldElement result(field);
  for (int i = 0; i < dimension; i++)
    result.polynomial[i] = (polynomial[i]>0) ? (modulus - polynomial[i]) : 0;
  return result;
}


/* ************************************************************************* */
bool operator==(const FieldElement& left, const FieldElement& right)
{
  // TODO
  // Make sure to manipulate fixed-length polynomial
  if (*left.field != *right.field) throw ErrorIncompatibleFields;

  if (polynomial_degree(left.polynomial) != polynomial_degree(right.polynomial))
    return false;

  for (int i = 0; i <= polynomial_degree(left.polynomial); i++)
    if (left.polynomial[i] != right.polynomial[i])
      return false;

  return true;
  // return left.polynomial == right.polynomial;
}


/* ************************************************************************* */
bool operator!=(const FieldElement& left, const FieldElement& right)
{
  return !(left == right);
}

bool operator<(const FieldElement& left, const FieldElement& right)
{
  if (*left.field != *right.field) throw ErrorIncompatibleFields;
  bool isless = false;
  for (int i = left.dimension-1; i >= 0; i--)
  {
    if (left.polynomial[i] < right.polynomial[i])
    {
      isless = true;
      break;
    }
    else if (left.polynomial[i] > right.polynomial[i])
    {
      isless = false;
      break;
    }
  }
  return isless;
}


/* ************************************************************************* */
bool operator>(const FieldElement& left, const FieldElement& right)
{
  if (*left.field != *right.field) throw ErrorIncompatibleFields;
  bool isgreater = false;
  for (int i=left.dimension-1; i>=0; i--)
  {
    if (left.polynomial[i] > right.polynomial[i])
    {
      isgreater = true;
      break;
    }
    else if (left.polynomial[i] < right.polynomial[i])
    {
      isgreater = false;
      break;
    }
  }
  return isgreater;
}


/* ************************************************************************* */
bool operator<=(const FieldElement& left, const FieldElement& right)
{
  if (*left.field != *right.field) throw ErrorIncompatibleFields;
  return !(left > right);
}


/* ************************************************************************* */
bool operator>=(const FieldElement& left, const FieldElement& right)
{
  if (*left.field != *right.field) throw ErrorIncompatibleFields;
  return !(left < right);
}


/* ************************************************************************* */
FieldElement operator+(const FieldElement& left, const FieldElement& right)
{
  if (*left.field != *right.field) throw ErrorIncompatibleFields;
  FieldElement result(left.field);
  for (int i = 0; i < min(left.polynomial.size(), right.polynomial.size()); i++)
    result.polynomial[i] = add_modular(left.polynomial[i],
        right.polynomial[i],
        left.modulus);
  return result;
}


/* ************************************************************************* */
FieldElement operator-(const FieldElement& left, const FieldElement& right)
{
  if (*left.field != *right.field) throw ErrorIncompatibleFields;
  FieldElement result(left.field);
  for (int i = 0; i < min(left.polynomial.size(), right.polynomial.size()); i++)
    result.polynomial[i] = subtract_modular(left.polynomial[i],
        right.polynomial[i],
        left.modulus);
  return result;
}


/* ************************************************************************* */
FieldElement operator*(const FieldElement& left, const FieldElement& right)
{
  if (*left.field != *right.field) throw ErrorIncompatibleFields;

  int dim = left.dimension;
  Fint modulus = left.modulus;

  vector<Fint> reductionpolynomial(dim+1);

  // multiply two polynomials left and right
  auto product = polynomial_multiply(left.polynomial, right.polynomial, 
      modulus);

  // copy the reduction polynomial from GaloisField associated with this 
  // FieldElement
  for (int i = 0; i < dim+1; i++)
    reductionpolynomial[i] = left.field->reduction_polynomial(i);

  // reduce the product polynomial by the reduction polynomial
  auto [Q, R] = polynomial_divide(product, reductionpolynomial, modulus);

  // R polynomial is the result to be returned
  FieldElement result(left.field, R);

  return result;
}



/* ************************************************************************* */
FieldElement operator/(const FieldElement& left, const FieldElement& right)
{
  if (*left.field != *right.field)
    throw ErrorIncompatibleFields;

  return  left * mult_inverse( right );
}



/* ************************************************************************* */
FieldElement operator*(Fint left, const FieldElement& right)
{
  FieldElement result(right.field);
  for (int i = 0; i < right.dimension; i++)
    result.polynomial[i] = multiply_modular(left, right.polynomial[i], 
        right.modulus);
  return result;
}



/* ************************************************************************* */
FieldElement operator*(const FieldElement& left, Fint right)
{
  FieldElement result(left.field);

  for (int i = 0; i < left.dimension; i++)
    result.polynomial[i] = multiply_modular(right, left.polynomial[i], 
        left.modulus);
  return result;
}



/* ************************************************************************* */
FieldElement operator/(const FieldElement& left, Fint right)
{
  FieldElement result(left.field);
  for (int i = 0; i < left.dimension; i++)
    result.polynomial[i] = divide_modular(left.polynomial[i], right, 
        left.modulus);
  return result;
}



/* ************************************************************************* */
ostream& operator<<(ostream& output, const FieldElement& right)
{
  output << "(";
  for (int i = 0; i < right.dimension; i++)
  {
    if (i < right.dimension-1)
      output << right.polynomial[i] << " ";
    else
      output << right.polynomial[i] << ")\n";
  }

  return output;
}



/* ************************************************************************* */
FieldElement mult_inverse(const FieldElement& FieldE)
{
  int dim = FieldE.dimension;
  Fint modulus = FieldE.modulus;

  vector<Fint> reductionpolynomial;

  // get reduction polynomial from GF associated with this FieldElement
  for (int i = 0; i < dim + 1; i++)
    reductionpolynomial.push_back(FieldE.field->reduction_polynomial(i));

  // determine inverse of FieldE.polynomial (calling it inverse)
  auto [dummy, inversepoly, gcd] = polynomial_extended_euclid(reductionpolynomial, 
      FieldE.polynomial, modulus);

  // check to make sure if gcd is 1
  int degGcd = polynomial_degree(gcd);

  // if degree of gcd is not 0, or if it is the zero polynomial, then throw 
  // exception
  if ( degGcd > 0 || gcd[0]==0 )
    throw ErrorNoInverse;

  if (gcd[0] > 1) // if gcd is a scalar that is not equal to 1, then
                  // inversepoly must be divided by gcd.
  {
    Fint invgcd = inverse_modular(gcd[0],modulus);
    for (int i = 0; i < dim; i++)
      inversepoly[i] = multiply_modular(inversepoly[i], invgcd, modulus);
  }

  // polynomial R is the result to be returned.
  FieldElement result(FieldE.field,inversepoly);

  return result;
}



/* ************************************************************************* */
FieldElement field_exp(const FieldElement& G, Fint m)
{
  FieldElement S( G.get_field() ); // create S to be in the same field as G

  S.set_coefficients(0); S[0] = 1; // set to S(x)=1

  if (m == 0) return S;

  FieldElement Gtemp( G.get_field() ); // create G in the same field as g
  Gtemp = G; // set Gtemp(x) = G(x)

  Fint mask = 1;

  if (m & mask)
    S = G;

  for (int i = 0; i < static_cast<signed>(8*sizeof(Fint)); i++)
  {
    Gtemp = Gtemp * Gtemp;
    mask <<= 1;
    if (m & mask)
      S = Gtemp * S;
  }

  return S;
}
