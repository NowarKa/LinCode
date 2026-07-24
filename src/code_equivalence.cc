#include "equivalence_code.hh"

#include <fstream>
#include <ostream>
#include <unordered_set>

/* ************************************************************************* */
auto write_lcequivalence_input(
    const unordered_set<pair<vector<int>, int>> &codes, 
    const Params &params,
    const string &filename) -> void
{
  std::ofstream file(filename);

  if (!file)
    throw std::runtime_error("Cannot open " + filename);

  size_t idx = 1;

  auto points = params.ps_kp1->get_all_points();

  for (const auto &[code, n] : codes)
  {
    auto q = params.field->get_order();

    file << "? " << params.k + 1 << " " << n << " " << q << " " << idx++ << "\n";

    bool use_commas = q > 10;

    for (size_t j = 0; j < params.k + 1; j++)
    {
      for (size_t i = 0; i < points.size(); i++)
      {
        if (code[i] == 0)
          continue;

        int value = points[i].get_coordinate(j).index();

        for (size_t l = 0; l < code[i]; l++)
        {
          if (use_commas)
          {
            if (i)
              file << ",";
            file << value;
          }
          else
          {
            file << value;
          }
        }
      }
      file << "\n";
    }
  }
}

/* ************************************************************************* */
auto call_lcequivalence() -> void
{
  std::string command = "printf '1\n10\n' | ./LCequivalence > /dev/null";
  auto rc = system(command.c_str());

  if (rc != 0)
    throw std::runtime_error("LCequivalence failed");
}

/* ************************************************************************* */
auto parse_generator_row(const std::string &line, uint32_t q)
    -> std::vector<int>
{
  std::vector<int> row;

  if (q <= 10)
  {
    for (char c : line)
      if (isdigit(c))
        row.push_back(c - '0');
  }
  else
  {
    std::stringstream ss(line);
    std::string token;

    while (std::getline(ss, token, ','))
      row.push_back(std::stoi(token));
  }

  return row;
}

/* ************************************************************************* */
auto parse_lcequivalence_output(const std::string &filename, Params &params)
    -> std::vector<LCode>
{
  std::ifstream file(filename);

  if (!file)
    throw std::runtime_error("Cannot open " + filename);

  std::vector<LCode> result;
  std::string line;

  while (std::getline(file, line))
  {
    if (line.empty())
      continue;

    if (line.rfind("END;", 0) == 0)
      break;

    if (line[0] != '?')
      continue;

    std::stringstream header(line);

    char question_mark;
    int k;
    int n;
    int q;
    int output_idx;
    int original_idx;

    header >> question_mark >> k >> n >> q >> output_idx >> original_idx;

    std::vector<std::vector<FieldElement>> matrix;

    for (int i = 0; i < k; i++)
    {
      std::getline(file, line);

      auto row_values = parse_generator_row(line, q);

      std::vector<FieldElement> row;

      for (auto value : row_values)
        row.push_back(params.field->get_element(value));

      matrix.push_back(std::move(row));
    }

    result.push_back(LCode(matrix));

    while (std::getline(file, line))
    {
      if (line.empty())
        break;

      if (line.rfind("==", 0) == 0)
        break;
    }
  }

  return result;
}

/* ************************************************************************* */
auto remove_equivalent_codes(
    const unordered_set<pair<vector<int>, int>> &candidates,
    Params &params) -> vector<LCode>
{
  system("mkdir -p RES_DIR0");

  write_lcequivalence_input(candidates, params, "RES_DIR0/EXAM");

  cout << "Calling LCequivalence...\n";
  call_lcequivalence();
  cout << "Done...\n";

  return parse_lcequivalence_output("RES_DIR0/EXAM_r", params);
}
