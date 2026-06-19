/**********************************************************************
 *
 * The implementation of Galois field is inspired from GaloisCPP 
 * library. Refer to https://github.com/saiedhk/GaloisCPP
 *
**********************************************************************/
#pragma once

#include "field.hh"
#include "polynomial_arith.hh"

#include <iostream>
#include <memory>
#include <vector>

using namespace std;
/**
  Defines a Galois Field element. Also defines arithmetic operations on Galois Field elements.

  In mathematics, an element of the Galois Field GF(p^k) is a polynomial of degree less than k
  with coefficients chosen from the set {0,1,2,..., p-1}.

  In this program, a FieldElement contains the following private data:
  @param field a shared pointer to the Galois Field this FieldElement belongs to
  @param modulus the prime modulus of the Galois Field
  @param dimension the dminesion of the Galois Field
  @param polynomial the array representing the polynomial describing this FieldElement
  */
class FieldElement
{
public:
  /**
    Contructs FieldElement without initializing its coefficients
    @param gf pointer to a Galois Field
    */
  FieldElement(const std::shared_ptr<const Field>& gf);


  /**
    Contructs a FieldElement and initializes its coefficients
    @param gf pointer to a Galois Field
    @param poly an external array containing coefficients (length of array = dimension)
    */
  FieldElement(const std::shared_ptr<const Field>& gf, const vector<Fint> poly);


  /**
    Constructs and clones a FieldElement.
    @param FieldE reference to a FieldElement
    */
  FieldElement(const FieldElement& FieldE) = default;

  /**
    Destructs a FieldElement and frees its allocated memory
    */
  ~FieldElement() = default;

  /**
    Returns the index of the element in its field.
    The index is computed as a0 + a1*p + ... + a_(k-1) * p^(k-1)
    where p is the characteristic of the field and the a_i are 
    the coefficients of the polynomial.
   */
  auto index() const -> uint64_t;


  /**
    Returns the pointer to the GF this FieldElement is associated with.
    */
  auto get_field() const -> shared_ptr<const Field>;


  /**
    Returns the i-th coefficient of this FieldElement's polynomial
    @param i the index to the i-th coefficient
    */
  auto get_coefficient(int i) const -> Fint;

  /**
    Returns the FieldElement's polynomial
    */
  auto get_polynomial(void) const -> const vector<Fint>&;


  /**
   * Is the FieldElement equal to field ZERO?
   */
  auto is_zero() const -> bool;

  /**
    Sets FieldElement to the ZERO polynomial.
    It sets all coefficients of FieldElement polynomial to zero.
    */
  auto set_zero(void) -> void;


  /**
    Sets FieldElement to the ONE polynomial.
    It sets all coefficients of FieldElement polynomial to zero except the zero-th
    degree coefficient which is set to 1.
    */
  auto set_one(void) -> void;


  /**
    Sets the i-th coefficient of FieldElement's polynomial to a specified constant
    @param i input the index to the i-th coefficient
    @param index the value of the i-th coefficient
    */
  auto set_coefficient(int index, Fint value) -> void;


  /**
    Sets all coefficients of FieldElement's polynomial to the same specified constant
    @param value input, the value for all coefficients
    */
  auto set_coefficients(Fint value) -> void;


  /**
    Allows accessing the i-th coefficient of FieldElement polynomial.

    This can be used as an l-value. It sets the i-th coefficient to a specified value.
    @param i (input), index
    @returns Reference of the i-th coefficient of FieldElement polynomial
    */
  auto operator[](int i) -> Fint&;


  /**
    Returns the i-th coefficient of FieldElement polynomial
    @param i (input), index
    @returns the i-th coefficient of FieldElement polynomial
    */
  auto operator[](int i) const -> Fint;


  /**
    Assigns/copies the values of another FieldElement to this FieldElement
    @param right FieldElement on the right of = sign
    @returns FieldElement on the left of = sign (this FieldElement)
    */
  FieldElement& operator=(const FieldElement& right) = default;


  /**
    Returns the negated copy of this FieldElement
    */
  FieldElement operator-();


  /**
    Compares two FieldElements, checks for equality
    @param right FieldElement on the right of == sign
    @param left FieldElement on the left of == sign
    @returns true if equal, false otherwise
    */
  friend bool operator==(const FieldElement& left, const FieldElement& right);


  /**
    Compares two FieldElements, checks for inequality
    @param right FieldElement on the right of != sign
    @param left FieldElement on the left of != sign
    @returns true if equal, false otherwise
    */
  friend bool operator!=(const FieldElement& left, const FieldElement& right);


  /**
    Compares two FieldElements based on their lexicographical ordering.
    @param right FieldElement on the right of < sign
    @param left FieldElement on the left of < sign
    @returns true if left < right, false otherwise
    */
  friend bool operator<(const FieldElement& left, const FieldElement& right);


  /**
    Compares two FieldElements based on their lexicographical ordering.
    @param right FieldElement on the right of > sign
    @param left FieldElement on the left of > sign
    @returns true if left > right, false otherwise
    */
  friend bool operator>(const FieldElement& left, const FieldElement& right);


  /**
    Compares two FieldElements based on their lexicographical ordering.
    @param right FieldElement on the right of <= sign
    @param left FieldElement on the left of <= sign
    @returns true if left <= right, false otherwise
    */
  friend bool operator<=(const FieldElement& left, const FieldElement& right);


  /**
    Compares two FieldElements based on their lexicographical ordering.
    @param right FieldElement on the right of >= sign
    @param left FieldElement on the left of >= sign
    @returns true if left >= right, false otherwise
    */
  friend bool operator>=(const FieldElement& left, const FieldElement& right);


  /**
    Adds two FieldElements
    @param right FieldElement on the right of + sign
    @param left FieldElement on the left of + sign
    @returns sum of two operands
    */
  friend FieldElement operator+(const FieldElement& left, const FieldElement& right);


  /**
    Subtracts two FieldElements
    @param right FieldElement on the right of - sign
    @param left FieldElement on the left of - sign
    @returns difference of two operands
    */
  friend FieldElement operator-(const FieldElement& left, const FieldElement& right);


  /**
    Multiplies two FieldElements
    @param right FieldElement on the right of * sign
    @param left FieldElement on the left of * sign
    @returns product of two operands
    */
  friend FieldElement operator*(const FieldElement& left, const FieldElement& right);


  /**
    Divides two FieldElements
    @param right FieldElement on the right of / sign
    @param left FieldElement on the left of / sign
    @returns quotient of two operands
    */
  friend FieldElement operator/(const FieldElement& left, const FieldElement& right);


  /**
    Multiplies a FieldElement by a scalar (i.e. a Fint)
    @param right FieldElement on the right of * sign
    @param left Fint on the left of * sign
    @returns product of two operands
    */
  friend FieldElement operator*(Fint left, const FieldElement& right);


  /**
    Multiplies a FieldElement by a scalar (i.e. a Fint)
    @param right Fint on the right of * sign
    @param left FieldElement on the left of * sign
    @returns product of two operands
    */
  friend FieldElement operator*(const FieldElement& left, Fint right);


  /**
    Divides a FieldElement by a scalar (i.e. a Fint)
    @param right Fint on the right of / sign
    @param left FieldElement on the left of / sign
    @returns quotient of two operands
    @throws ErrorDivideByZero
    */
  friend FieldElement operator/(const FieldElement& left, Fint right);


  /**
    Outputs a FieldElement to the standard output
    @param right FieldElement on the right of << sign
    @param output an output stream on the left of << sign
    @returns output stream
    */
  friend ostream& operator<<(ostream& output, const FieldElement& right);


  /**
    Determines the multiplicative inverse of a FieldElement.
    @param FieldE input FieldElement
    @returns FieldElement (the multiplicative inverse of FieldE)
    @throws ErrorNoInverse, ErrorInvalidArraySize
    */
  friend FieldElement mult_inverse(const FieldElement& FieldE);


  /**
    Raises a FieldElement G to a positive integer power m

    It uses the repeated square-and-multiply algorithm for exponentiation in GF(p^k).

    (See Chapter 2, Algorithm 2.227 in Handbook of Applied Cryptography,
    by Menezes, van Oorschot, Vanstone, CRC Press, 1996.)

    @param G an input FieldElement
    @param m an input positive integer
    @returns G raised to the power m
    @throws ErrorInvalidArraySize, ErrorDivideByZero, ErrorIncompatibleFields
    */
  friend FieldElement field_exp(const FieldElement& G, Fint m);


private:
  shared_ptr<const Field> field;
  Fint modulus;
  int dimension;
  vector<Fint> polynomial;
};


template<> struct std::hash<FieldElement> {
  std::size_t operator()(FieldElement const& fe) const noexcept {
    size_t seed = 0;
    // size_t h_field = std::hash<Field>{}(*fe.get_field());
    // boost::hash_combine(seed, h_field);

    for (const auto& c : fe.get_polynomial())
      boost::hash_combine(seed, hash<Fint>{}(c));

    return seed;
  }
};


inline auto FieldElement::get_field() const -> shared_ptr<const Field>
{
  return field;
}


inline auto FieldElement::get_coefficient(int index) const -> Fint
{
  return polynomial[index];
}


inline auto FieldElement::set_coefficient(int index, Fint value) -> void
{
  polynomial[index] = value % modulus;
}


inline auto FieldElement::set_coefficients(Fint value) -> void
{
  for (int i=0; i<dimension; i++)
    polynomial[i] = value % modulus;
}


inline auto FieldElement::is_zero() const -> bool
{
  return polynomial_degree(polynomial) == -1;
}

inline auto FieldElement::get_polynomial(void) const -> const vector<Fint>&
{
  return polynomial;
}
