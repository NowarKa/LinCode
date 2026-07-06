#include "constraints_phase1.hh"
#include "constructed_codes_table.hh"
#include "field.hh"
#include "field_vector.hh"
#include "linear_code.hh"
#include "projective_space.hh"

#include <cstdlib>
#include <filesystem>
#include <memory>
#include <string>
#include <vector>
#include <iostream>
#include <CLI/CLI.hpp>


int main (int argc, char *argv[]) 
{
  Params params;
  CLI::App app{"LinCode - Linear Code Classification Tool"};

  app.add_option("--delta", params.delta,
      "Restrict the classification to linear codes whose codeword weights "
      "are divisible by delta (default: 1).");

  app.add_option("--field-file", params.field_file,
      "Path to the file defining the finite field.")
    ->check(CLI::ExistingFile);

  app.add_option("--minimum-weight", params.minimum_weight,
      "Minimum minimum distance of the linear codes to classify (default: 1).");

  app.add_option("--maximum-weight", params.maximum_weight,
      "Maximum allowed minimum distance of the linear codes to "
      "classify (default: INT_MAX).");

  app.add_flag("--check-feasibility", params.check_feasibility,
      "Check the feasibility of solutions with SCIP before enumerating them.");

  app.add_flag("--save-results", params.save_results,
      "Save the classified codes to disk.");

  app.add_option("--load", params.k, "Load previously saved results and " 
      "initialize the queue with all saved codes of dimension [n', k], " 
      "where k is fixed.");

  app.add_option("--ubn", params.upper_bound_n, "Classify linear codes of "
      "dimension [n, k] for all lengths n <= ubn (default: 6).");

  app.add_option("--nb-threads", params.nb_threads, "The number of threads "
      "used for parallelization.");

  CLI11_PARSE(app, argc, argv);

  const Field GF4 = Field(2, 2, {1,1,1});
  const Field GF2 = Field(2, 1, {1,1});
  const Field GF3 = Field(3, 1, {1,1});

  auto GF2_ptr = make_shared<const Field>(GF2);
  auto GF3_ptr = make_shared<const Field>(GF3);
  auto GF4_ptr = make_shared<const Field>(GF4);

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

  // Constructing Field
  if (params.field_file != "")
    params.field = 
      make_shared<const Field>(Field::parse_field_file(params.field_file));
  else 
    params.field = GF2_ptr;

  // Loading results and initializing the queue
  if (params.k != 0)
  {
    cout << "Loading existing results...\n";
    constructed_codes.load(params.field, 2);
    // constructed_codes.load_queue(params.k, extended_code, params.field);
  }


  else 
  {
    cout << "Adding one-dimensional linear codes...\n";
    for (size_t n = 2; n < params.upper_bound_n; n++)
    {
      auto ps = ProjectiveSpace(n, params.field);

      // #pragma omp parallel for
      for (auto& p : ps.get_all_points())
      {
        auto w = hamming_weight(p.get_coordinates());
        if (w >= params.minimum_weight && w % params.delta == 0)
        {
          // cout << "n is now " << n << endl;
          // cout << "Code added with weight " << w << endl;
          auto code = LCode({p.get_coordinates()});
          // extended_code.push(code);
          constructed_codes.insert_code(code);
        } // end if
      } // end for
    } // end for
    constructed_codes.save(params.field);
    cout << "Done\n";
  } // end else

  /*
  // Engine of classification
  while (!extended_code.empty()) 
  {
    auto c = extended_code.front();
    extended_code.pop();

    size_t n = c.get_nb_columns();
    size_t k = c.get_nb_rows();

    if (n > params.upper_bound_n)
      continue;

    system("clear");

    cout << "==================================================================\n";
    cout << "=== ****************** Intermediate results ****************** ===\n";
    cout << "==================================================================\n";
    cout << constructed_codes;
    cout << "Extending code (" << n << ", " << k << ")" << endl;
    cout << c << endl;
    cout << c.canonical_form() << endl;
    cout << "==================================================================\n";
    cout << "=== ************************** END *************************** ===\n";
    cout << "==================================================================\n";

    extend_code(c, params, extended_code, constructed_codes);
  }
  */


  // TODO: To change later
  int ub_k = 7;
  int k = params.k;

  while (k < ub_k)
  {
    system("clear");

    cout << "==================================================================\n";
    cout << "=== ****************** Intermediate results ****************** ===\n";
    cout << "==================================================================\n";
    cout << constructed_codes;
    cout << "==================================================================\n";
    cout << "=== ************************** END *************************** ===\n";
    cout << "==================================================================\n";

    auto jobs = constructed_codes.split_by_weight_enumerator(k, params.nb_threads);
    vector<ConstructedCodesTable> local_tables(params.nb_threads);

    #pragma omp parallel for
    for (int t = 0; t < params.nb_threads; ++t)
    {
      for (auto& code : jobs[t])
        extend_code(code,
            params,
            t,
            local_tables[t]);
    }

    constructed_codes.merge_list(local_tables, k+1);
    k++;
  }

  cout << constructed_codes;

  if (params.save_results)
    constructed_codes.save(params.field);

  return 0;
}
