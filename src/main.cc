#include "code_extension.hh"
#include "constructed_codes_table.hh"
#include "equivalence_code.hh"
#include "field.hh"
#include "hyperplane.hh"
#include "linear_code.hh"
#include "projective_space.hh"

#include <CLI/CLI.hpp>
#include <cassert>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <execution>
#include <filesystem>
#include <iostream>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

int main(int argc, char *argv[])
{
  chrono::steady_clock::time_point begin = chrono::steady_clock::now();

  Params params;
  CLI::App app{"LinCode - Linear Code Classification Tool"};

  app.add_option(
      "--delta", params.delta,
      "Restrict the classification to linear codes whose codeword weights "
      "are divisible by delta (default: 1).");

  app.add_option("--field-file", params.field_file,
                 "Path to the file defining the finite field.")
      ->check(CLI::ExistingFile);

  app.add_option(
      "--minimum-weight", params.minimum_weight,
      "Minimum minimum distance of the linear codes to classify (default: 1).");

  app.add_option("--maximum-weight", params.maximum_weight,
                 "Maximum allowed minimum distance of the linear codes to "
                 "classify (default: INT_MAX).");

  app.add_flag(
      "--check-feasibility", params.check_feasibility,
      "Check the feasibility of solutions with SCIP before enumerating them.");

  app.add_flag("--save-results", params.save_results,
               "Save the classified codes to disk.");

  app.add_option(
      "--load", params.k,
      "Load previously saved results and "
      "initialize the queue with all saved codes of dimension [n', k], "
      "where k is fixed.");

  app.add_option("--ubn", params.upper_bound_n,
                 "Classify linear codes of "
                 "dimension [n, k] for all lengths n <= ubn (default: 6).");

  app.add_option("--nb-threads", params.nb_threads,
                 "The number of threads "
                 "used for parallelization.");

  CLI11_PARSE(app, argc, argv);

  const Field GF4 = Field(2, 2, {1, 1, 1});
  const Field GF2 = Field(2, 1, {1, 1});
  const Field GF3 = Field(3, 1, {1, 1});
  const Field GF8 = Field(2, 3, {1, 0, 1, 1});

  auto GF2_ptr = make_shared<const Field>(GF2);
  auto GF3_ptr = make_shared<const Field>(GF3);
  auto GF4_ptr = make_shared<const Field>(GF4);
  auto GF8_ptr = make_shared<const Field>(GF8);

  ConstructedCodesTable constructed_codes;

  // Adding the directory in which solvediophant will work and
  // the solutions will be added
  {
    if (!filesystem::exists("temp/"))
      system("mkdir temp");

    for (int i = 0; i < params.nb_threads; i++)
      if (!filesystem::exists("temp/job_" + to_string(i)))
        system(("mkdir temp/job_" + to_string(i)).c_str());
  }

  // TODO: To change later
  int ub_k = 8;
  int k = params.k;

  // Constructing Field
  if (params.field_file != "")
    params.field =
        make_shared<const Field>(Field::parse_field_file(params.field_file));

  else
    params.field = GF2_ptr;
  //  params.field = GF4_ptr;
  // params.field = GF8_ptr;
  // params.field = GF3_ptr;

  // Loading results and initializing the queue
  if (params.k != 0)
    cout << "Loading existing results...\n";

  else
  {
    cout << "Adding one-dimensional linear codes...\n";

    for (size_t n = 2; n < params.upper_bound_n; ++n)
    {
      for (size_t d = params.minimum_weight; d <= n; ++d)
      {
        if (d % params.delta != 0)
          continue;

        vector<FieldElement> g(n, params.field->get_element(0));

        for (size_t i = 0; i < d; ++i)
          g[i] = params.field->get_element(1);

        auto code = LCode({g});
        constructed_codes.insert_code(code);
      }
    }
    cout << "Done\n";
    k = 1;
  } // end else

  long long int time_equivalence = 0;
  long long int time_sd = 0;
  long long int time_scip = 0;
  params.time_equivalence = make_shared<long long int>(time_equivalence);
  params.time_sd = make_shared<long long int>(time_sd);
  params.time_scip = make_shared<long long int>(time_scip);

  while (k < ub_k)
  {
    // system("clear");

    cout << "=================================================================="
            "\n";
    cout << "=== ****************** Intermediate results ****************** ==="
            "\n";
    cout << "=================================================================="
            "\n";
    cout << constructed_codes;
    cout << "=================================================================="
            "\n";
    cout << "=== ************************** END *************************** ==="
            "\n";
    cout << "=================================================================="
            "\n";

    cout << "Splitting jobs...\n";
    auto jobs =
        constructed_codes.split_by_weight_enumerator(k, params.nb_threads);

    bool is_empty = true;

    for (int t = 0; t < params.nb_threads; ++t)
      if (!jobs[t].empty())
      {
        is_empty = false;
        break;
      }

    if (is_empty)
    {
      k++;
      continue;
    }

    cout << "Constructing projective space...\n";
    auto projective_space_k = ProjectiveSpace(k, params.field);
    auto projective_space_kp1 = ProjectiveSpace(k + 1, params.field);

    auto h_kp1 = Hyperplane::get_all_hyperplanes(projective_space_kp1);

    cout << "Constructing containing map...\n";
    params.ps_k = make_shared<const ProjectiveSpace>(projective_space_k);
    params.ps_kp1 = make_shared<const ProjectiveSpace>(projective_space_kp1);
    params.h_kp1 = make_shared<const vector<Hyperplane>>(h_kp1);
    params.containing_map =
        make_shared<const unordered_map<Hyperplane, vector<int>>>(
            Hyperplane::containing_map(params.ps_kp1, params.h_kp1));

    unordered_map<ProjectivePoint, int> indexes;

    for (int i = 0; i < params.ps_kp1->get_all_points().size(); i++)
      indexes.insert({params.ps_kp1->get_all_points()[i], i});

    params.indexes = make_shared<unordered_map<ProjectivePoint, int>>(indexes);

    vector<LCode> candidates;
    params.candidates = make_shared<vector<LCode>>(candidates);

    cout << "Extending...\n";
#pragma omp parallel for
    for (int t = 0; t < params.nb_threads; ++t)
    {
      for (auto &code : jobs[t])
      {
        extend_code(code, params, t);
      }
    }

    chrono::steady_clock::time_point begin_time_equivalence =
        chrono::steady_clock::now();

    auto reps = remove_equivalent_codes(*params.candidates, params);

    chrono::steady_clock::time_point end_time_equivalence =
        chrono::steady_clock::now();

    *params.time_equivalence +=
        chrono::duration_cast<chrono::milliseconds>(end_time_equivalence -
                                                    begin_time_equivalence)
            .count();

    for (auto &code : reps)
      constructed_codes.insert_code(code);

    k++;
  } // end while

  if (params.save_results)
    constructed_codes.save(params.field);

  chrono::steady_clock::time_point end = chrono::steady_clock::now();
  long long int execution_time =
      chrono::duration_cast<chrono::seconds>(end - begin).count();

  auto equivalence_time_sec =
      static_cast<float>(*params.time_equivalence) / 1000;
  auto sd_time_sec = static_cast<float>(*params.time_sd) / 1000;
  auto scip_time_sec = static_cast<float>(*params.time_scip) / 1000;
  auto ratio_equivalence_time =
      equivalence_time_sec / static_cast<float>(execution_time);
  auto ratio_sd_time = sd_time_sec / static_cast<float>(execution_time);
  auto ratio_scip_time = scip_time_sec / static_cast<float>(execution_time);

  auto execution_time_string =
      "Execution time = " + to_string(execution_time) + "[s]\n";

  auto equivalence_time_string =
      "Ellapsed time in equivalence test = " + to_string(equivalence_time_sec) +
      "[s]\n";

  auto equivalence_ratio_string =
      "Ratio = " + to_string(ratio_equivalence_time) + "\n";

  auto sd_time_string =
      "Ellapsed time in Solvediophant = " + to_string(sd_time_sec) + "[s]\n";

  auto sd_ratio_string = "Ratio = " + to_string(ratio_sd_time) + "\n";

  auto scip_time_string =
      "Ellapsed time in SCIP = " + to_string(scip_time_sec) + "[s]\n";

  auto scip_ratio_string = "Ratio = " + to_string(ratio_scip_time) + "\n";

  filesystem::path file = "data/general.txt";
  ofstream out(file);

  out << constructed_codes;
  out << "---\n\n";
  out << execution_time_string << equivalence_time_string
      << equivalence_ratio_string << sd_time_string << sd_ratio_string
      << scip_time_string << scip_ratio_string;

  cout << execution_time_string << equivalence_time_string
       << equivalence_ratio_string << sd_time_string << sd_ratio_string;

  return 0;
}
