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


class ConstructedCodesTable
{
public:
  ConstructedCodesTable() = default;

  ~ConstructedCodesTable() = default;

  auto insert_code(LCode& code) -> void;
  
  auto contains_code(LCode& code) -> bool;

  auto save(const filesystem::path& directory = "data") const -> void;

  auto load(int upper_bound_n = INT_MAX, 
      const filesystem::path& directory = "data") -> void;

  friend ostream& operator<<(ostream& output, 
      const ConstructedCodesTable& right);

private:
  unordered_map<pair<size_t, size_t>, unordered_set<string>> table_;
};
