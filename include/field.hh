/**********************************************************************
 *
 * The implementation of Galois field is inspired from GaloisCPP 
 * library. Refer to https://github.com/saiedhk/GaloisCPP
 *
**********************************************************************/
#pragma once

#include "typedefs.hh"

#include <boost/functional/hash.hpp>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <memory>
#include <vector>
#include <cmath>

class FieldElement;

/**
   This class defines a Galois Field of characteristic p and dimesion k (for GF(p^k)).
   It stores the following (private) data about the field:
   @param the characteristics p, which MUST be a prime number.
   @param dimension The dimension (k) of Galois Field (a nonzero positive integer).
   @param order The order (p^k) of Galois Field.
   @param reductpoly A polynomial of degree k over the prime field GF(p). This polynomial
   MUST be irreducible. The program does not check this. YOU MUST correctly choose
   this polynomial.

   The attributes of a GaloisField are set when it is instantiated.
   Once instantiated, the attributes of cannot be changed.
*/
class Field : public std::enable_shared_from_this<Field> 
{
public:
  /**
    Constructs a Galois Field (GF)
    @param p a positive prime integer (the characteristic of GF)
    @param dim a positive integer (the dimension of GF)
    @param poly an array of length dim(the reduction polynomial of GF)
    (must be an irreducible polynomial over the prime field of GF)
    */
  Field(uint32_t p, uint32_t k, std::vector<Fint> reductpoly) : 
  characteristic(p), dimension(k), reductpoly(reductpoly)
  {
    order = std::pow(p, k);
  };


  /**
    Destructs Field.
    */
  ~Field() = default;


  /**
   * Returns characteristic of GF.
   */
  uint32_t get_characteristic() const;


  /**
   * Returns dimension of GF.
   */
  uint32_t get_dimension() const;


  /**
   * Returns order of GF.
   */
  uint64_t get_order() const;


  /**
    Returns i-th coefficient of reduction polynomial of GF
    */
  Fint reduction_polynomial(int i) const;


  /**
    Returns the reduction polynomial of GF
    */
  auto get_reduction_polynomial(void) const -> std::vector<Fint>;

  
  /**
    Returns i-th element of GF (lexicographical order).
    */
  auto get_element(uint64_t i) const
    -> FieldElement;


  /**
    Returns all elements of GF, sorted (lexicographical order).
    */
  auto get_all_elements() const -> std::vector<FieldElement>;


  /**
   * Constructs Field from a .txt file. The format is the following one:
   * p <- the characteristic of the Field
   * a0 a1 a2 ... an <- the coefficients of the reduction polynomial.
   *
   * @param file Path to the .txt file.
   */
  static auto parse_field_file(const std::filesystem::path& file) -> Field;


  // I might change this function.
  // By unicity of finite field, it is sufficient to test 
  // whether the two orders are equals.
  friend bool operator==(const Field& left, const Field& right)
  {
    return left.order == right.order && 
      left.reductpoly == right.reductpoly;
  };


  friend bool operator!=(const Field& left, const Field& right)
  {
    return !(left == right);
  };

private:
  uint32_t characteristic;
  uint32_t dimension;
  uint64_t order;
  std::vector<Fint> reductpoly;
};


template<> struct std::hash<Field> {
  std::size_t operator()(Field const& f) const noexcept {
    size_t seed = 0;

    std::size_t h_characteristic = std::hash<int>{}(f.get_characteristic());
    std::size_t h_dimension = std::hash<int>{}(f.get_dimension());

    boost::hash_combine(seed, h_characteristic);
    boost::hash_combine(seed, h_dimension);

    for (int i = 0; i < f.get_dimension(); i++)
      boost::hash_combine(seed, hash<Fint>{}(f.reduction_polynomial(i)));

    return seed;
  }
};


inline auto Field::get_characteristic() const -> uint32_t
{
  return characteristic;
}


inline auto Field::get_dimension() const -> uint32_t
{
  return dimension;
}


inline auto Field::get_order() const -> uint64_t
{
  return order;
}


inline Fint Field::reduction_polynomial(int i) const
{
  return reductpoly[i];
}

inline auto Field::get_reduction_polynomial() const -> std::vector<Fint>
{
  return reductpoly;
}
