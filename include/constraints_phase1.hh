#pragma once

#include "constructed_codes_table.hh"
#include "linear_code.hh"

#include <vector>

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
};


struct ExtensionParams
{
  Params& params;
  int a; 
  int b; 
  int r; 
  int thread_id;
  LCode lcode; 
};


auto generate_equations_phase1(int delta, int a, int b, int r, 
    shared_ptr<const Field>& field, LCode g, vector<LCode>& extended_code) -> void;

auto extend_code(LCode& code, Params &params, int thread_id, 
    ConstructedCodesTable& constructed_codes) -> void;
