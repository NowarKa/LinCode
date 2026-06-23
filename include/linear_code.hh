#pragma once

#include "field_element.hh"
#include "projective_point.hh"
#include "hyperplane.hh"

#include <cstddef>
#include <cstdint>
#include <unordered_map>
#include <vector>

using namespace std;

/**
  Defines a linear code. Also defines some operations on linear codes.  

  In coding theory, a linear code of length n is a k-dimensional linear 
  subspace of the vector space GF(q)^n, where GF(q) is the finite field 
  of q elements. A common representation of an [n, k]-code is as the row-span 
  of a k × n matrix over GF(q) - called generator matrix.

  In this program, an LCode contains the following private data:
  @param columns_ The set of the generator matrix's *COLUMNS*.
  @param rows_ The set of the generator matrix's *ROWS*.
  @param canonical_form_ The generator matrix's canonical form.
 */
class LCode 
{
public:
  /**
    Constructs an LCode without initialization.
   */
  LCode() = default;


  /**
    Constructs an LCode and initializes its generator matrix.

    @param rows A rectangular matrix represented by the set of its *ROWS*.
   */
  LCode(const vector<vector<FieldElement>>& rows);


  /**
    Destructs an LCode.
   */
  ~LCode() = default;


  /**
    Returns the number LCode's columns.
   */
  auto get_nb_columns() const -> size_t;


  /**
    Returns the number LCode's rows.
   */
  auto get_nb_rows() const -> size_t;


  /**
    Returns the i-th row of LCode's generator matrix.

    @param i A given integer i.
   */
  auto get_row(size_t i) const -> vector<FieldElement>;


  /**
    Returns the LCode's minimum column multiplicity.
   */
  auto get_minimum_column_multiplicity() -> int;


  /**
    Returns LCode's multiset.
    
    Starting from a generator matrix of a linear code we obtain a multiset 
    of points by considering the multiset of one-dimensional subspaces 
    spanned by the respective columns.
   */
  auto to_multiset() const -> unordered_map<ProjectivePoint, uint32_t>;


  /**
    Returns the minimal distance of LCode.
   */
  auto weight(const vector<Hyperplane>& hyperplanes) const -> uint32_t;


  /**
    Returns the minimal distance of LCode.
   */
  auto weight() const -> uint32_t;


  /**
   * Returns true if all weights of the codewords of LCode is divisible by 
   * delta (aka LCode is delta-divisible), false otherwise.
   *
   * @param delta A specified integer.
   */
  auto is_codeword_weight_divisible(int delta) -> bool;


  /**
    Returns LCode's canonical form.
   */
  auto canonical_form(const string& sage_binary = "sage") -> const string;


  /**
    Returns u.G where G is LCode's generator matrix and u a given vector.

    @param u The specified vector.
   */
  auto encode_column_vector(const vector<FieldElement>& u) const
    -> vector<FieldElement>;


  /**
    Returns the number of columns of LCode's generator matrix whose spanned 
    one-dimensional subspace is equal to p - a given ProjectivePoint.

    @param p The specified ProjectivePoint.
   */
  auto nb_columns_belong_to(const ProjectivePoint& p) const -> uint32_t;


  /**
   * Access for mutable matrix
   */
  auto operator()(size_t i, size_t j) -> FieldElement&;


  /**
   * Access for immutable matrix
   */
  auto operator()(size_t i, size_t j) const -> const FieldElement&;


  /**
   * Constructs LCode from canonical form.
   *
   * @param canonical_form A given canonical form.
   * @param field A given field.
   */
  static auto from_canonical_form(const string& canonical_form,
      shared_ptr<const Field> field) -> LCode;


  /**
    Outputs an LCode to the standard output.

    @param right The LCode on the right of << sign
    @param output An output stream on the left of << sign
    @returns output stream
   */
  friend ostream& operator<<(ostream& output, const LCode& right);

private:
  vector<vector<FieldElement>> columns_;
  vector<vector<FieldElement>> rows_;
  string canonical_form_ = "";
  int minimum_column_multiplicity_ = 0;
};


auto construct_from_columns(const vector<vector<FieldElement>>& columns) 
  -> LCode;

inline auto LCode::get_nb_columns() const -> size_t
{
  return columns_.size();
}


inline auto LCode::get_nb_rows() const -> size_t
{
  return rows_.size();
}


inline auto LCode::get_row(size_t i) const -> vector<FieldElement>
{
  return rows_[i];
}


inline auto LCode::operator()(size_t i, size_t j) const -> const FieldElement&
{
  return rows_[i][j];
}


inline auto LCode::operator()(size_t i, size_t j) -> FieldElement&
{
  return rows_[i][j];
}
