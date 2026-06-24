#pragma once

#include <iostream>
#include <tuple>

#include "typedefs.hh"

using namespace std;

/**
  Implements the Extended Euclid Algorithm.

  It computes the greatest common divisor d of two positive integers a and b
  and determine two integers x and y such that d = ax + by.

  (See Chapter 2, Algorithm 2.107 in Handbook of Applied Cryptography,
  by Menezes, van Oorschot, Vanstone, CRC Press, 1996.)

  @param a input integer (a > 0)
  @param b input integer (b > 0)
*/
auto extended_euclid(Fint a, Fint b) -> tuple<Fint, Fint, Fint>;

/**
  Adds two integers modulo m
  @param a input integer
  @param b input integer
  @param mod input integer (modulus)
  @returns sum modulo mod
*/
inline auto add_modular(Fint a, Fint b, Fint mod) -> Fint
{
  return (a + b) % mod;
}


/**
  Subtracts two integers modulo m
  @param a input integer
  @param b input integer
  @param mod input integer (modulus)
  @returns difference modulo mod
  */
inline auto subtract_modular(Fint a, Fint b, Fint mod) -> Fint
{
  return (a + (mod - b)) % mod;
}

/**
  Multiplies two integers modulo m
  @param a input integer
  @param b input integer
  @param mod input integer (modulus)
  @returns product modulo mod
  */
inline auto multiply_modular(Fint a, Fint b, Fint mod) -> Fint
{
  return (a * b) % mod;
}


/**
  Determines the multiplicative inverse of an integer modulo mod
  @param a input integer (a>0)
  @param mod input integer (modulus, mod>0)
  @returns b the inverse of a such that a * b = 1 modulo mod
  @throws ErrorNoInverse if mod is not prime relative to a
  */
inline auto inverse_modular(Fint a, Fint mod) -> Fint
{
  auto [x, y, d] =  extended_euclid(a, mod);
  if (d != 1) {
    cout << a << " has no inverse mod " << mod << "\n";
    throw ErrorNoInverse;
  }
  return ( (x>0) ? x : (mod + x) );
}


/**
  Divides two integers modulo mod
  @param a input integer (dividend)
  @param b input integer (divisor)
  @param mod input integer (modulus)
  @returns q quotient such that a = q * b modulo mod
  @throws ErrorNoInverse if mod is not prime relative to a, or ErrorDivideByZero if b=0
  */
inline auto divide_modular(Fint a, Fint b, Fint mod) -> Fint
{
  if (b==0) {
    cout << "Dividing by zero!!!\n";
    throw ErrorDivideByZero;
  }
  return (a * inverse_modular(b,mod)) % mod;
}
