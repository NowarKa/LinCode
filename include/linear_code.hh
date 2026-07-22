#pragma once

#include "field_element.hh"
#include "hyperplane.hh"
#include "projective_point.hh"
#include "projective_space.hh"

#include <chrono>
#include <cstddef>
#include <cstdint>
#include <memory>
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
  LCode(const vector<vector<FieldElement>> &rows);

  /**
   * Constructs LCode from multiset.
   *
   * @param multiset A given multiset.
   */
  LCode(const unordered_map<ProjectivePoint, uint32_t> &multiset);

  LCode(
      shared_ptr<const unordered_map<Hyperplane, vector<int>>> &containing_map,
      const vector<int> &solution, shared_ptr<const ProjectiveSpace> ps);

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
   * Returns LCode's weight enumerator.
   */
  auto
  get_weight_enumerator(shared_ptr<const ProjectiveSpace> ps = nullptr) const
      -> vector<int>;

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
  auto weight(const vector<Hyperplane> &hyperplanes) const -> uint32_t;

  /**
    Returns the minimal distance of LCode.
   */
  auto minimum_distance(shared_ptr<const ProjectiveSpace> ps = nullptr) const
      -> uint32_t;

  /**
   * Returns true if all weights of the codewords of LCode is divisible by
   * delta (aka LCode is delta-divisible), false otherwise.
   *
   * @param delta A specified integer.
   */
  auto is_codeword_weight_divisible(int delta) -> bool;

  /**
    Returns u.G where G is LCode's generator matrix and u a given vector.

    @param u The specified vector.
   */
  auto encode_column_vector(const vector<FieldElement> &u) const
      -> vector<FieldElement>;

  /**
    Returns the number of columns of LCode's generator matrix whose spanned
    one-dimensional subspace is equal to p - a given ProjectivePoint.

    @param p The specified ProjectivePoint.
   */
  auto nb_columns_belong_to(const ProjectivePoint &p) const -> uint32_t;

  /**
   * Constructs an LCode resulting by deleting all columns equal
   * to a given ProjectivePoint.
   *
   * @param p A given ProjectivePoint.
   */
  auto remove_projective_point(const ProjectivePoint &p) const -> LCode;

  /**
   * Returns true if LCode is projective, i.e., m(P) \in {0,1} for
   * all P.
   */
  auto is_projective() const -> bool;

  /**
   * Access for mutable matrix
   */
  auto operator()(size_t i, size_t j) -> FieldElement &;

  /**
   * Access for immutable matrix
   */
  auto operator()(size_t i, size_t j) const -> const FieldElement &;

  /**
   * Constructs LCode from columns.
   *
   * @param columns The set of *COLUMNS*
   */
  static auto
  construct_from_columns(const vector<vector<FieldElement>> &columns) -> LCode;

  /**
    Outputs an LCode to the standard output.

    @param right The LCode on the right of << sign
    @param output An output stream on the left of << sign
    @returns output stream
   */
  friend ostream &operator<<(ostream &output, const LCode &right);

  friend bool operator==(const LCode &left, const LCode &right);

private:
  vector<vector<FieldElement>> columns_;
  vector<vector<FieldElement>> rows_;
  int minimum_column_multiplicity_ = -1;
  int minimum_distance_ = -1;
  mutable vector<int> weight_enumerator_;
  mutable unordered_map<ProjectivePoint, uint32_t> multiset_;
};

// Hash function for LCode in order to use unordered_map.
// Two LCodes are equals if they are equivalent.
template <> struct std::hash<LCode>
{
  size_t operator()(const LCode &code) const
  {
    std::size_t seed = 0;

    for (size_t i = 0; i < code.get_nb_rows(); i++)
    {
      std::size_t h = 0;

      for (size_t j = 0; j < code.get_nb_columns(); j++)
      {
        h ^= std::hash<int>{}(code(i, j).index()) + 0x9e3779b9 + (h << 6) +
             (h >> 2);
      }

      seed ^= h + 0x9e3779b9 + (seed << 6) + (seed >> 2);
    }

    return seed;
  }
};

inline auto LCode::get_nb_columns() const -> size_t { return columns_.size(); }

inline auto LCode::get_nb_rows() const -> size_t { return rows_.size(); }

inline auto LCode::get_row(size_t i) const -> vector<FieldElement>
{
  return rows_[i];
}

inline auto LCode::operator()(size_t i, size_t j) const -> const FieldElement &
{
  return rows_[i][j];
}

inline auto LCode::operator()(size_t i, size_t j) -> FieldElement &
{
  return rows_[i][j];
}
