#pragma once

#include "linear_code.hh"

#include <climits>
#include <cstddef>
#include <unordered_map>
#include <unordered_set>
#include <iostream>
#include <filesystem>

using namespace std;

// Hash function for pair<size_t, size_t> in order to use unordered_map.
template<> struct std::hash<pair<size_t, size_t>> {
  size_t operator()(const pair<size_t, size_t>& p) const noexcept 
  {
    size_t h1 = hash<size_t>{}(p.first);
    size_t h2 = hash<size_t>{}(p.second);
    return h1 * 31 + h2;
  }
};

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
  auto insert_code(LCode& code) -> void;
  

  /**
    Returns true if ConstructedCodesTable contains a given linear code, and 
    false otherwise.

    @param code A specified linear code.
   */
  auto contains_code(LCode& code) -> bool;


  /**
    Saves ConstructedCodesTable's content to disk. By default, files are stored 
    in `data/` folder.

    @param directory Directrory's name.
   */
  auto save(const filesystem::path& directory = "data") const -> void;


  /**
    Loads ConstructedCodesTable's content from files.
   */
  auto load(int upper_bound_n = INT_MAX, 
      const filesystem::path& directory = "data") -> void;


  /**
    Outputs a ConstructedCodesTable to the standard output
    @param right ConstructedCodesTable on the right of << sign
    @param output an output stream on the left of << sign
    @returns output stream
    */
  friend ostream& operator<<(ostream& output, 
      const ConstructedCodesTable& right);

private:
  unordered_map<pair<size_t, size_t>, unordered_set<string>> table_;
};
