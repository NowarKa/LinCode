#pragma once

#include "field.hh"
#include "linear_code.hh"
#include "utils.hh"

#include <climits>
#include <cstddef>
#include <filesystem>
#include <iostream>
#include <memory>
#include <queue>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

using namespace std;

/**
  Defines the structure that stores classified linear codes. For each cople
  (n, k) it assocites the set of classified [n, k] linear codes in their
  CANONICAL FORM.

  It stores the following private data:
  @table_ The unordered_map that stores the classified linear code.
 */
class ConstructedCodesTable
{
public:
  /**
    Constructs a ConstructedCodesTable.
   */
  ConstructedCodesTable() = default;

  /**
    Destructs a ConstructedCodesTable.
   */
  ~ConstructedCodesTable() = default;

  /**
    Inserts a given linear code.

    @param code A specified linear code.
   */
  auto insert_code(LCode &code) -> void;

  /**
    Returns true if ConstructedCodesTable contains a given linear code, and
    false otherwise.

    @param code A specified linear code.
   */
  auto contains_code(LCode &code, shared_ptr<long long int> time_sage = nullptr)
      -> bool;

  /**
    Returns true if ConstructedCodesTable contains a given linear code, and
    false otherwise. Update the field `minimum_weight_enumerator_extension_`
    of ConstructedCodesTable's LCodes in order to keep the smallest weight
    enumerator (lexicographical order).

    @param code A specified linear code.
   */
  // auto contains_code(LCode &basic_code, LCode &extended_code) -> bool;

  /**
    Saves ConstructedCodesTable's content to disk. By default, files are stored
    in `data/` folder.

    @param directory Directrory's name.
   */
  auto save(shared_ptr<const Field> field,
            const filesystem::path &directory = "data") const -> void;

  /**
    Loads ConstructedCodesTable's content from files.
   */
  auto load(shared_ptr<const Field> field,
      const filesystem::path& directory = "data") -> pair<size_t, size_t>;

  /**
   * Loads all codes of dimension [n', k] to queue to start classifcation.
   *
   * @param k A given integer.
   * @param extended_code A given queue.
   */
  auto load_queue(int k, queue<LCode> &extended_code,
                  shared_ptr<const Field> field) -> void;

  /**
   * Inserts all [n, k] (for all n) LCodes of a list of ConstructedCodesTables.
   *
   * @param results A given list of ConstructedCodesTables.
   * @param k A specified integer.
   */
  auto merge_list(const vector<ConstructedCodesTable> &results, int k) -> void;

  /**
   * Splits all [n, k] LCodes into `nb_threads`
   * groups according to their weight enumerators. Codes having the
   * same weight enumerator are assigned to the same group.
   *
   * @param k A specified integer.
   * @param nb_thread The number of groups (typically one per thread).
   */
  auto split_by_weight_enumerator(int k, int nb_thread)
      -> vector<vector<LCode>>;

  /**
   * Splits all [n, k] LCodes into `nb_threads` groups according to 
   * their number of columns.
   *
   * @param k A specified integer.
   * @param nb_thread The number of groups (typically one per thread).
   */
  auto split(int k, int nb_threads) -> vector<vector<LCode>>;

  /**
    Outputs a ConstructedCodesTable to the standard output
    @param right ConstructedCodesTable on the right of << sign
    @param output an output stream on the left of << sign
    @returns output stream
    */
  friend ostream &operator<<(ostream &output,
                             const ConstructedCodesTable &right);

private:
  unordered_map<pair<size_t, size_t>, unordered_set<LCode>> table_;
};
