#include "constructed_codes_table.hh"
#include "field.hh"
#include "linear_code.hh"

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <execution>
#include <fstream>
#include <iomanip>
#include <queue>
#include <string>
#include <utility>
#include <vector>

/* ************************************************************************* */
auto ConstructedCodesTable::insert_code(LCode &code) -> void
{
  size_t n = code.get_nb_columns();
  size_t k = code.get_nb_rows();

  table_[{n, k}].insert(code);

  return;
}

/* ************************************************************************* */
auto ConstructedCodesTable::contains_code(LCode &code,
                                          shared_ptr<long long int> time_sage)
    -> bool
{
  size_t n = code.get_nb_columns();
  size_t k = code.get_nb_rows();

  if (table_.find({n, k}) == table_.end())
    return false;

  for (auto &s : table_.at({n, k}))
    if (s == code)
      return true;

  return false;
}

/* ************************************************************************* */
ostream &operator<<(ostream &output, const ConstructedCodesTable &right)
{
  size_t max_i = 0;
  size_t max_j = 0;

  for (const auto &[key, value] : right.table_)
  {
    max_i = max(max_i, key.first);
    max_j = max(max_j, key.second);
  }

  constexpr int W = 6;

  output << setw(W) << "";
  for (size_t j = 0; j <= max_j; ++j)
    output << setw(W) << j;
  output << '\n';

  for (size_t i = 0; i <= max_i; ++i)
  {
    output << setw(W) << i;

    for (size_t j = 0; j <= max_j; ++j)
    {
      auto it = right.table_.find({i, j});
      size_t count = (it == right.table_.end()) ? 0 : it->second.size();
      output << setw(W) << count;
    }

    output << '\n';
  }

  return output;
}

/* ************************************************************************* */
auto ConstructedCodesTable::save(shared_ptr<const Field> field,
                                 const std::filesystem::path &directory) const
    -> void
{
  namespace fs = filesystem;

  fs::create_directories(directory);

  for (const auto &[nk, codes] : table_)
  {
    auto [n, k] = nk;

    fs::path file = directory / (to_string(n) + "_" + to_string(k) + ".txt");

    ofstream out(file);

    if (!out)
      throw std::runtime_error("Cannot open " + file.string());

    out << "characteristic=" << field->get_characteristic() << '\n';
    out << "reduction_polynomial=";

    const auto &poly = field->get_reduction_polynomial();

    for (size_t i = 0; i < poly.size(); ++i)
    {
      if (i)
        out << ' ';

      out << poly[i];
    }

    out << "\n\n";

    for (const auto &code : codes)
    {
      out << code.get_nb_columns() << " " << code.get_nb_rows() << " "
          << code.minimum_distance() << endl;
      out << "Is projective: " << code.is_projective() << endl;
      out << code.to_string();
      out << "---\n";
    }
  } // end for
}

/* ************************************************************************* */
auto ConstructedCodesTable::load_queue(int k, queue<LCode> &extended_code,
                                       shared_ptr<const Field> field) -> void
{
  for (auto &[p, codes] : table_)
  {
    auto [n, l] = p;

    if (l != k)
      continue;

    for (auto &code : codes)
      extended_code.push(code);
  }

  return;
}

/* ************************************************************************* */
auto ConstructedCodesTable::merge_list(
    const vector<ConstructedCodesTable> &results, int k0) -> void
{
  for (const auto &result : results)
    for (const auto &[nk, codes] : result.table_)
    {
      auto [n, k] = nk;

      if (k != k0)
        continue;

      for (const auto &code : codes)
        table_[nk].insert(code);
    }

  return;
}

/* ************************************************************************* */
auto ConstructedCodesTable::split_by_weight_enumerator(int k0, int nb_thread)
    -> vector<vector<LCode>>
{
  unordered_map<string, vector<LCode>> groups;

  for (const auto &[nk, codes] : table_)
  {
    auto [n, k] = nk;

    if (k != k0)
      continue;

    for (const auto &code : codes)
    {
      string key;

      for (auto x : code.get_weight_enumerator())
      {
        key += to_string(x);
        key += ',';
      }

      groups[key].push_back(code);
    } // end for
  } // end for

  vector<vector<LCode>> jobs(nb_thread);
  vector<size_t> loads(nb_thread, 0);

  for (auto &[key, group] : groups)
  {
    auto it = min_element(loads.begin(), loads.end());
    size_t thread_id = distance(loads.begin(), it);

    jobs[thread_id].insert(jobs[thread_id].end(), group.begin(), group.end());

    loads[thread_id] += group.size();
  }

  return jobs;
}

/* ************************************************************************* */
auto ConstructedCodesTable::split(int k, int nb_threads) 
  -> vector<vector<LCode>>
{
  vector<vector<LCode>> result(nb_threads);

  for (const auto& [key, value] : table_) {
    if (key.second != k)
      continue;

    size_t i = 0;
    size_t nb_added_codes = 0;
    size_t partition_size = value.size()/nb_threads;

    for (const auto &code : value)
    {
      if (nb_added_codes > partition_size)
      {
        nb_added_codes = 0;
        i = (i + 1) % nb_threads;
      }
      result[i].push_back(code);
      nb_added_codes++;
    }

  }

  return result;
}

/* ************************************************************************* */
struct ParsedFieldInfo
{
  int characteristic;
  vector<Fint> reduction_polynomial;
};

auto parse_field_header(ifstream &in) -> ParsedFieldInfo
{
  ParsedFieldInfo result;

  string line;

  if (!getline(in, line))
    throw runtime_error("Missing characteristic");

  const string characteristic_prefix = "characteristic=";

  if (!line.starts_with(characteristic_prefix))
    throw runtime_error("Invalid characteristic line");

  result.characteristic = stoi(line.substr(characteristic_prefix.size()));

  if (!getline(in, line))
    throw runtime_error("Missing reduction polynomial");

  const string polynomial_prefix = "reduction_polynomial=";

  if (!line.starts_with(polynomial_prefix))
    throw runtime_error("Invalid reduction polynomial line");

  istringstream iss(line.substr(polynomial_prefix.size()));

  Fint coeff;

  while (iss >> coeff)
    result.reduction_polynomial.push_back(coeff);

  getline(in, line);

  return result;
}

auto ConstructedCodesTable::load(
    shared_ptr<const Field> field,
    const filesystem::path& directory) -> pair<size_t, size_t>
{
  namespace fs = filesystem;

  size_t max_n = 0;
  size_t max_k = 0;

  for (const auto& entry : fs::directory_iterator(directory))
  {
    if (!entry.is_regular_file())
      continue;

    auto stem = entry.path().stem().string();

    size_t pos = stem.find('_');

    if (pos == string::npos)
      continue;

    size_t n = stoull(stem.substr(0, pos));
    size_t k = stoull(stem.substr(pos + 1));

    max_n = max(max_n, n);
    max_k = max(max_k, k);

    /*
    if (n > upper_bound_n || k > upper_bound_k)
      continue;
    */

    ifstream in(entry.path());

    if (!in)
      throw runtime_error("Cannot open " + entry.path().string());

    /*
    ParsedFieldInfo field_info =
        parse_field_header(in);
    */

    string line;
    string current;

    if (!in)
        throw std::runtime_error("Cannot open file.");

    getline(in, line); // characteristic=...
    getline(in, line); // reduction_polynomial=...
    getline(in, line); // empty line

    while (true)
    {
      // Reading line "n k d"
      if (!std::getline(in, line))
        break;

      if (line.empty())
        continue;

      std::istringstream iss(line);

      int n, k, d;
      if (!(iss >> n >> k >> d))
        throw std::runtime_error("Invalid format for line n k d.");


      // Reading "Is projective: ..."
      if (!std::getline(in, line))
        throw std::runtime_error("Inexepected end of line.");

      vector<vector<FieldElement>> G;
      G.reserve(k);

      // Reading generator matrix
      for (int i = 0; i < k; ++i)
      {
        if (!std::getline(in, line))
          throw std::runtime_error("Inexepected end of file.");

        std::istringstream row(line);
        std::vector<FieldElement> r;
        r.reserve(n);

        int x;
        while (row >> x)
          r.push_back(field->get_element(x));

        if ((int)r.size() != n)
          throw std::runtime_error("Incorrect number of columns.");

        G.push_back(std::move(r));
      }

      LCode code = LCode(G);
      table_[{n,k}].insert(code);

      // Reading ---
      if (!std::getline(in, line))
        break;

      if (line != "---")
      {
        cout << "In " << stem << endl;
        throw std::runtime_error("Expected ---.");
      }
    }

  }

  return {max_n, max_k};
}
