#include "code_extension.hh"
#include "field_vector.hh"
#include "linear_code.hh"
#include "phase0.hh"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <execution>
#include <fstream>
#include <iostream>
#include <iterator>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

/* ************************************************************************* */
auto passes_corollary_8(const LCode &parent, LCode &child,
                        ExtensionParams &ext_params) -> bool
{
  int r = child.get_minimum_column_multiplicity();

  auto parent_we = parent.get_weight_enumerator(ext_params.params.ps_k);
  auto min_we = parent_we;

  bool first = true;
  auto multiset = child.to_multiset();

  for (const auto &[point, multiplicity] : multiset)
  {
    if (multiplicity != r)
      continue;

    auto candidate_parent = child.remove_projective_point(point);
    auto we = candidate_parent.get_weight_enumerator(ext_params.params.ps_kp1);

    if (we < min_we)
      return false;
  }

  return parent_we == min_we;
}

/* ************************************************************************* */
auto check_solution(const vector<int> &solution, 
    ExtensionParams &ext_params, size_t np_kp1) -> bool
{
  // Check x_<e_i> >= 1
  for (auto &i : *ext_params.params.unit_vector_index)
    if (solution[i] < ext_params.r)
      return false;

  for (size_t i = 0; i < np_kp1; i++)
    if (solution[i] > 0 && solution[i] < ext_params.r)
      return false;

  /*
  for (const auto &[h, v] : *ext_params.params.containing_map) {
    int sum = 0;

    for (auto i : v)
      sum += solution[i];

    if (n + ext_params.r - sum < ext_params.params.minimum_weight ||
        n + ext_params.r - sum > ext_params.params.maximum_weight)
      return false;
  }
  */
  // Check y_H \in {0, ..., b - a}
  /*
  for (size_t i = np_kp1; i < solution.size(); i++)
    if (solution[i] > ext_params.b - ext_params.a)
      return false;
  */

  return true;
}

/* ************************************************************************* */
auto check_lcode(LCode &code, ExtensionParams &ext_params) -> bool
{
  /*
  if (code.get_minimum_column_multiplicity() != ext_params.r)
    return false;


  return true;
  */
  return code.minimum_distance() >= ext_params.params.minimum_weight;
}

/* ************************************************************************* */
auto parse_solution(ExtensionParams &ext_params) -> vector<vector<int>>
{
  auto solution_filename =
      "temp/job_" + to_string(ext_params.thread_id) + "/solutions.txt";

  vector<vector<int>> solutions;
  ifstream file(solution_filename);

  if (!file)
    throw runtime_error("Failed to open file " + solution_filename);

  string line;
  while (getline(file, line))
  {
    if (line.empty())
      continue;

    istringstream iss(line);
    vector<int> row;
    int value;

    while (iss >> value)
      row.push_back(value);

    if (!row.empty())
      solutions.push_back(std::move(row));
  }

  return solutions;
}

/* ************************************************************************* */
auto construct_lcode_from_solution(const vector<int> &solution,
                                   ExtensionParams &ext_params) -> LCode
{
  return LCode(ext_params.params.containing_map, solution,
               ext_params.params.ps_kp1);
}

/* ************************************************************************* */
auto call_solvediophant(ExtensionParams &ext_params) -> void
{
  auto solution_filename =
      "temp/job_" + to_string(ext_params.thread_id) + "/solutions.txt";

  auto problem_filename =
      "temp/job_" + to_string(ext_params.thread_id) + "/problem.txt";

  auto command = "./sd2 -o" + solution_filename + " " + problem_filename;

  system(command.c_str());
}

/* ************************************************************************* */
auto build_equations(ExtensionParams &ext_params) -> Equations
{
  auto field = ext_params.params.field;

  auto n = ext_params.lcode.get_nb_columns();
  auto k = ext_params.lcode.get_nb_rows();

  auto field_elements = field->get_all_elements();

  auto p_kp1 = ext_params.params.ps_kp1->get_all_points();
  auto p_k = ext_params.params.ps_k->get_all_points();

  auto np_k = p_k.size();
  auto np_kp1 = p_kp1.size();
  auto nh_kp1 = ext_params.params.h_kp1->size();

  vector<vector<int>> A;

  // Adding equation 2
  for (size_t i = 0; i < nh_kp1; i++)
  {
    auto line = vector<int>(np_kp1 + nh_kp1, 0);
    line[np_kp1 + i] = ext_params.params.delta;

    for (auto j :
         ext_params.params.containing_map->at(ext_params.params.h_kp1->at(i)))
      line[j] = 1;

    A.push_back(line);
  } // end for

  // Adding equation 3
  for (size_t i = 0; i < np_k; i++)
  {
    auto line = vector<int>(np_kp1 + nh_kp1, 0);

    for (const auto &fe : field_elements)
    {
      auto p = p_k[i].concatenate(fe);
      auto index = ext_params.params.indexes->at(p);
      line[index] = 1;
    }

    A.push_back(line);
  }

  auto line = vector<int>(np_kp1 + nh_kp1, 0);
  auto ekp1 = get_unit_vector(k + 1, k, field);
  auto p_ekp1 = ProjectivePoint(ekp1);
  auto idx = ext_params.params.indexes->at(p_ekp1);
  line[idx] = 1;
  A.push_back(line);

  // Constructing d (such that A*x = d)
  vector<int> d =
      vector<int>(nh_kp1 + np_k + 1,
                  n + ext_params.r - ext_params.a * ext_params.params.delta);
  auto m = ext_params.lcode.to_multiset();

  for (size_t i = 0; i < np_k; i++)
  {
    try
    {
      d[nh_kp1 + i] = m.at(p_k[i]);
    }
    catch (out_of_range)
    {
      d[nh_kp1 + i] = 0;
    }
  }

  // Adding the special line <u> = 0
  d[nh_kp1 + np_k] = ext_params.r;

  // Upper bounds
  vector<int> upper_bounds =
      vector<int>(nh_kp1 + np_kp1, ext_params.b - ext_params.a);
  for (size_t i = 0; i <= np_kp1; i++)
    // I can do better here
    upper_bounds[i] = n + ext_params.r;

  Equations equations_system = Equations({A, d, upper_bounds});

  return equations_system;
}

/* ************************************************************************* */
auto generate_equations(ExtensionParams &ext_params) -> void
{
  auto k = ext_params.params.k;
  auto p_kp1 = ext_params.params.ps_kp1->get_all_points();
  auto np_kp1 = p_kp1.size();

  auto system = build_equations(ext_params);

  if (ext_params.params.check_feasibility &&
      !check_feasibility(system, ext_params.r, ext_params.params.time_scip))
    return;

  // Writing problem to file
  ofstream problem_phase1;

  auto problem_filename =
      "temp/job_" + to_string(ext_params.thread_id) + "/problem.txt";

  problem_phase1.open(problem_filename);

  // Adding problem's dimensions
  problem_phase1 << system.A.size() << " " << system.A[0].size() << " " << 1
                 << endl;

  // Adding matrices A and d
  for (size_t i = 0; i < system.A.size(); i++)
  {
    for (size_t j = 0; j < system.A[0].size(); j++)
      problem_phase1 << system.A[i][j] << " ";

    problem_phase1 << system.d[i] << endl;
  }

  // Adding upper bounds
  problem_phase1 << "BOUNDS " << system.upper_bounds.size() << endl;

  for (size_t i = 0; i < system.upper_bounds.size(); i++)
    problem_phase1 << system.upper_bounds[i] << " ";

  problem_phase1.close();

  chrono::steady_clock::time_point begin = chrono::steady_clock::now();

  call_solvediophant(ext_params);
  auto solutions = parse_solution(ext_params);

  chrono::steady_clock::time_point end = chrono::steady_clock::now();

  *ext_params.params.time_sd +=
      chrono::duration_cast<chrono::milliseconds>(end - begin).count();

  vector<LCode> candidates;

  auto n = ext_params.lcode.get_nb_columns();

  cout << "Filtering solutions...\n";
  for (const auto &solution : solutions)
  {
    if (!check_solution(solution, ext_params, np_kp1))
      continue;

    auto m = n + ext_params.r;

    if (k == 1)
    {
      auto code = construct_lcode_from_solution(solution, ext_params);
      if (code.get_nb_columns() < n || 
          !check_lcode(code, ext_params))
        continue;
      m = code.get_nb_columns();
    }
    ext_params.params.candidates->insert({solution, m});
  } // end for
  cout << "Done.\n";

  return;
}

/* ************************************************************************* */
auto extend_code(LCode &code, Params &params, int thread_id) -> void
{
  if (code.get_nb_rows() >= params.upper_bound_n)
    return;

  for (int r = 1; r <= code.get_minimum_column_multiplicity(); r++)
  {
    int m = code.get_nb_columns() + r;

    if (m > params.upper_bound_n || (m <= params.load_n && code.get_nb_rows() <= params.load_k))
      continue;

    int b = min(m/params.delta, params.maximum_weight/params.delta);
    ExtensionParams ext_params = {
        params, max(params.minimum_weight, 1), b, r, thread_id, code};
    generate_equations(ext_params);
  }

  return;
}
