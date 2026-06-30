#include "constraints_phase1.hh"
#include "constructed_codes_table.hh"
#include "field_vector.hh"
#include "hyperplane.hh"
#include "phase0.hh"
#include "projective_space.hh"
#include "linear_code.hh"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <iterator>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>
#include <sstream>


/* ************************************************************************* */
auto get_position_unit_vector(size_t i, uint32_t q) -> size_t
{
  // The first is (0, ..., 1)
  return (pow(q, i) - 1)/(q - 1);
}


/* ************************************************************************* */
auto check_solution(const vector<int>& solution, 
    uint32_t q, uint32_t k, 
    size_t np_kp1, int a, int b) -> bool
{
  // Check x_<e_i> >= 1
  for (size_t i = 0; i < k; i++)
    if (solution[get_position_unit_vector(i, q)] == 0)
      return false;

  // Check y_H \in {0, ..., b - a}
  for (size_t i = np_kp1; i < solution.size(); i++)
    if (solution[i] > b - a)
      return false;

  return true;
}


/* ************************************************************************* */
auto check_lcode( 
    ConstructedCodesTable& constructed_codes,
    LCode& code, ExtensionParams &ext_params) 
  -> bool
{
  if (code.get_nb_columns() > ext_params.params.upper_bound_n)
    return false;

  if (code.get_minimum_column_multiplicity() != ext_params.r)
    return false;

  if (constructed_codes.contains_code_and_update(ext_params.lcode, code))
    return false;

  return code.minimum_distance() >= ext_params.params.minimum_weight && 
     code.minimum_distance() <= ext_params.params.maximum_weight && 
    code.is_codeword_weight_divisible(ext_params.params.delta);
}


/* ************************************************************************* */
auto parse_solution(ExtensionParams &ext_params) -> vector<vector<int>>
{
  auto solution_filename = 
    "temp/job_" + 
    to_string(ext_params.thread_id) + 
    "/solutions.txt";

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
auto construct_lcode_from_solution(const vector<int>& solution, 
    const vector<ProjectivePoint>& p_kp1) -> LCode
{
  vector<vector<FieldElement>> matrix;
  auto np_kp1 = p_kp1.size();

  for (size_t i = 0; i < np_kp1; i++) 
  {
    auto x_p = solution[i];

    for (size_t j = 0; j < x_p; j++) 
      matrix.push_back(p_kp1[i].get_coordinates());
  }

  return LCode::construct_from_columns(matrix);
}


/* ************************************************************************* */
auto call_solvediophant(ExtensionParams &ext_params) -> void
{
  auto solution_filename = 
    "temp/job_" + 
    to_string(ext_params.thread_id) + 
    "/solutions.txt";

  auto problem_filename = 
    "temp/job_" + 
    to_string(ext_params.thread_id) + 
    "/problem.txt";

  auto command = 
    "./sd2 -o" + 
    solution_filename + 
    " " + problem_filename;

  system(command.c_str());
}


/* ************************************************************************* */
auto build_equations(ExtensionParams &ext_params, 
    ConstructedCodesTable& constructed_codes) -> Equations
{
  auto field = ext_params.params.field;

  auto n = ext_params.lcode.get_nb_columns();
  auto k = ext_params.lcode.get_nb_rows();

  auto projective_space_k = ProjectiveSpace(k, field);
  auto projective_space_kp1 = ProjectiveSpace(k+1, field);
  auto field_elements = field->get_all_elements();

  auto h_kp1 = get_all_hyperplanes(projective_space_kp1);
  auto p_kp1 = projective_space_kp1.get_all_points();
  auto p_k = projective_space_k.get_all_points();

  auto np_k = p_k.size();
  auto np_kp1 = p_kp1.size();
  auto nh_kp1 = h_kp1.size();

  vector<vector<int>> A;

  cout << "Adding equation 2\n";
  // Adding equation 2
  for (size_t i = 0; i < nh_kp1; i++)
  {
    auto line = vector<int>(np_kp1 + nh_kp1, 0);
    line[np_kp1 + i] = ext_params.params.delta;

    for (size_t j = 0; j < np_kp1; j++)
      if (h_kp1[i].contains(p_kp1[j]))
        line[j] = 1;

    A.push_back(line);
  } // end for

  cout << "Adding equation 3\n";
  // Adding equation 3
  for(size_t i = 0; i < np_k; i++)
  {
    auto line = vector<int>(np_kp1 + nh_kp1, 0);

    for (const auto& fe : field_elements)
    {
      auto p = p_k[i].concatenate(fe);
      auto index = distance(p_kp1.begin(), find(p_kp1.begin(), p_kp1.end(), p));
      line[index] = 1;
    }

    A.push_back(line);
  }

  auto line = vector<int>(np_kp1 + nh_kp1, 0);
  auto ekp1 = get_unit_vector(k+1, k, field);
  auto p_ekp1 = ProjectivePoint(ekp1);
  auto idx = distance(
      p_kp1.begin(),
      find(p_kp1.begin(), p_kp1.end(), p_ekp1));
  line[idx] = 1;
  A.push_back(line);

  cout << "Constructing d\n";
  // Constructing d (such that A*x = d)
  vector<int> d = vector<int>(nh_kp1 + np_k + 1, n - ext_params.a*ext_params.params.delta);
  auto m = ext_params.lcode.to_multiset();


  for(size_t i = 0; i < np_k; i++)
  {
    try {
      d[nh_kp1 + i] = m.at(p_k[i]);
    } catch (out_of_range) {
      d[nh_kp1 + i] = 0;
    }
  }

  cout << "Adding special line to d\n";
  // Adding the special line <u> = 0
  d[nh_kp1 + np_k] = ext_params.r;

  // Upper bounds
  vector<int> upper_bounds = vector<int>(nh_kp1 + np_kp1, ext_params.b - ext_params.a);
  for (size_t i = 0; i <= np_kp1; i++)
    // I can do better here
    upper_bounds[i] = n + ext_params.r;

  Equations equations_system = Equations({A, d, upper_bounds});

  return equations_system;
}


/* ************************************************************************* */
auto generate_equations_phase1(ExtensionParams &ext_params, 
    ConstructedCodesTable& constructed_codes)
  -> void
{ 
  auto k = ext_params.lcode.get_nb_rows();
  auto projective_space_kp1 = ProjectiveSpace(k+1, ext_params.params.field);
  auto p_kp1 = projective_space_kp1.get_all_points();
  auto np_kp1 = p_kp1.size();

  auto system = build_equations(ext_params, constructed_codes);

  if (ext_params.params.check_feasibility && 
      !check_feasibility(system, ext_params.r))
    return;

  // Writing problem to file
  ofstream problem_phase1;

  auto problem_filename = 
    "temp/job_" + 
    to_string(ext_params.thread_id) + 
    "/problem.txt";

  problem_phase1.open(problem_filename);

  // Adding problem's dimensions
  problem_phase1 << system.A.size() << " " << system.A[0].size() << " " << 1 << endl;

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

  call_solvediophant(ext_params);
  auto solutions = parse_solution(ext_params);

  for (const auto& solution : solutions)
  {
    if (!check_solution(solution, ext_params.params.field->get_order(), k, 
          np_kp1, ext_params.a, ext_params.b))
      continue;

    auto code = construct_lcode_from_solution(solution, p_kp1);

    if (check_lcode(constructed_codes, code, ext_params))
    {
      constructed_codes.insert_code(code);
    }
  } // end for

  return;
}


/* ************************************************************************* */
auto extend_code(LCode& code, Params &params, int thread_id, 
    ConstructedCodesTable& constructed_codes) -> void
{
  for (int r = 1; r <= code.get_minimum_column_multiplicity() + 1; r++)
  {
    if (!code.should_extend(r))
      continue;

    int b = code.get_nb_columns() + r;
    ExtensionParams ext_params = {params, params.a, b, r, thread_id, code};
    generate_equations_phase1(ext_params, constructed_codes);
  }

  return;
}
