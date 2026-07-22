#pragma once

#include "field.hh"
#include "hyperplane.hh"
#include "linear_code.hh"

#include <climits>
#include <cstddef>
#include <math.h>
#include <memory>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

using std::pair;
using std::shared_ptr;
using std::string;
using std::vector;

template <> struct std::hash<std::vector<int>>
{
  size_t operator()(const std::vector<int> &v) const noexcept
  {
    size_t h = 0;

    for (int i = 0; i < v.size(); i++)
      h += pow(2, i) * hash<int>{}(v[i]);

    return h;
  }
};

// Hash function for pair<size_t, size_t> in order to use unordered_map.
template <> struct std::hash<pair<size_t, size_t>>
{
  size_t operator()(const pair<size_t, size_t> &p) const noexcept
  {
    size_t h1 = hash<size_t>{}(p.first);
    size_t h2 = hash<size_t>{}(p.second);
    return h1 * 31 + h2;
  }
};

struct Equations
{
  vector<vector<int>> A;
  vector<int> d;
  vector<int> upper_bounds;
};

struct Params
{
  int delta = 1;
  int a = 0;
  string field_file = "";
  int minimum_weight = 1;
  int maximum_weight = INT_MAX;
  int upper_bound_n = 6;
  int k = 0;
  int nb_threads = 1;
  bool check_feasibility = false;
  bool save_results = false;
  shared_ptr<const Field> field;
  shared_ptr<const ProjectiveSpace> ps_k;
  shared_ptr<const ProjectiveSpace> ps_kp1;
  shared_ptr<const vector<Hyperplane>> h_kp1;
  shared_ptr<long long int> time_equivalence = nullptr;
  shared_ptr<long long int> time_sd = nullptr;
  shared_ptr<long long int> time_scip = nullptr;
  shared_ptr<vector<LCode>> candidates;
  shared_ptr<const unordered_map<Hyperplane, vector<int>>> containing_map;
  shared_ptr<const unordered_map<ProjectivePoint, int>> indexes;
};

struct ExtensionParams
{
  Params &params;
  int a;
  int b;
  int r;
  int thread_id;
  LCode &lcode;
};
