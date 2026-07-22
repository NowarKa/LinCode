#pragma once

#include <cassert>
#include <tuple>
#include <utility>
#include <vector>

#include "modular_arith.hh"

/**
   Implements addition of two polynomials over a finite ring

   It adds polynomial A to polynomial B.  (i.e., S = A + B)

   @param A input, addend polynomial (an array of Fint elements)
   @param B input, addend polynomial (an array of Fint elements)
   @param modulus input, a positive integer
   @throws ErrorInvalidArraySize
*/
auto polynomial_add(const vector<Fint> &A, const vector<Fint> &B, Fint modulus)
    -> vector<Fint>;

/**
   Implements subtraction of two polynomials over a finite ring.

   It subtracts polynomial B from polynomial A. (i.e., D = A - B))

   @param A input, minuend polynomial (an array of Fint elements)
   @param B input, subtrahend polynomial (an array of Fint elements)
   @param modulus input, a positive integer
   @throws ErrorInvalidArraySize
*/
auto polynomial_subtract(const vector<Fint> &A, const vector<Fint> &B,
                         Fint modulus) -> vector<Fint>;

/**
   Implements multiplication of two polynomials.

   It multiplies polynomial A by polynomial B.

   @param A input, multiplicand polynomial (an array of Fint elements)
   @param B input, multiplier polynomial (an array of Fint elements)
   @param modulus input, a positive integer
   @throws ErrorInvalidArraySize
*/
auto polynomial_multiply(const vector<Fint> &A, const vector<Fint> &B,
                         Fint modulus) -> vector<Fint>;

/**
   Implements long division of two polynomials over a prime field

   It divides polynomial A by polynomial B and determines the resulting
   quotient Q and remainder R (such that A = BQ + R).

   @param A input, dividend polynomial (an array of Fint elements)
   @param B input, divisor polynomial (an array of Fint elements)
   @param modulus input, a prime positive integer
   @throws ErrorInvalidArraySize or ErrorDivideByZero
*/
auto polynomial_divide(const vector<Fint> &A, const vector<Fint> &B,
                       Fint modulus) -> pair<vector<Fint>, vector<Fint>>;

/**
   Implements Extended Euclid algorithm on two polynomials over a prime field.

   (See Chapter 2, Algorithm 2.221 in Handbook of Applied Cryptography,
   by Menezes, van Oorschot, Vanstone, CRC Press, 1996.)

   @param A input polynomial (an array of Fint elements)
   @param B input polynomial (an array of Fint elements)
   @param modulus input, a prime positive integer
   @throws ErrorInvalidArraySize
*/
auto polynomial_extended_euclid(const vector<Fint> &A, const vector<Fint> &B,
                                Fint modulus)
    -> tuple<vector<Fint>, vector<Fint>, vector<Fint>>;

/**
  Prints a polynomial
  @param A input array containing the source polynomial
*/
auto polynomial_print(const vector<Fint> &A) -> void;

/**
   Determines degree of a polynomial (an input array)
   @param poly input array of type Fint containing polynomial coefficients
   @returns degree of polynomial (returns -1 for the all zero polynomial)
*/
inline auto polynomial_degree(const vector<Fint> &poly) -> int {
  for (int i = poly.size() - 1; i >= 0; i--)
    if (poly[i])
      return i;
  return -1;
}

/**
   Sets a polynomial to the all zero polynomial
   @param poly input array (polynomial)
*/
inline auto polynomial_set_zero(vector<Fint> &poly) -> void {
  for (int i = 0; i < poly.size(); i++)
    poly[i] = 0;
}
