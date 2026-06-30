#include "constructed_codes_table.hh"
#include "field.hh"
#include "linear_code.hh"

#include <algorithm>
#include <fstream>
#include <iomanip>
#include <memory>
#include <queue>
#include <string>


/* ************************************************************************* */
auto ConstructedCodesTable::insert_code(LCode& code) -> void
{
  size_t n = code.get_nb_columns();
  size_t k = code.get_nb_rows();

  table_[{n, k}].insert(code);

  return;
}


/* ************************************************************************* */
auto ConstructedCodesTable::contains_code(LCode& code) -> bool
{
  size_t n = code.get_nb_columns();
  size_t k = code.get_nb_rows();

  if (table_.find({n, k}) == table_.end())
    return false;

  auto cn_code = code.canonical_form();

  for (auto& s : table_.at({n, k})) 
    if (s.canonical_form() == cn_code)
      return true;

  return false;
}


/* ************************************************************************* */
auto ConstructedCodesTable::contains_code_and_update(LCode &basic_code, 
    LCode &extended_code) -> bool
{
  size_t n = extended_code.get_nb_columns();
  size_t k = extended_code.get_nb_rows();

  if (table_.find({n, k}) == table_.end())
    return false;

  auto cn_code = extended_code.canonical_form();

  for (auto& s : table_.at({n, k})) 
    if (s.canonical_form() == cn_code)
    {
      s.updates_minimum_weight_enumerator_extension(basic_code);
      return true;
    }

  return false;
}


/* ************************************************************************* */
ostream& operator<<(ostream& output, const ConstructedCodesTable& right)
{
  size_t max_i = 0;
  size_t max_j = 0;

  for (const auto& [key, value] : right.table_)
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
      size_t count = (it == right.table_.end())
        ? 0
        : it->second.size();
      output << setw(W) << count;
    }

    output << '\n';
  }

  return output;
}


/* ************************************************************************* */
auto ConstructedCodesTable::save(shared_ptr<const Field> field, 
    const std::filesystem::path& directory) const -> void
{
  namespace fs = filesystem;

  fs::create_directories(directory);

  for (const auto& [nk, codes] : table_)
  {
    auto [n, k] = nk;

    fs::path file =
      directory /
      (to_string(n)
       + "_"
       + to_string(k)
       + ".txt");

    ofstream out(file);

    if (!out)
      throw std::runtime_error("Cannot open " + file.string());

    out << "characteristic="
        << field->get_characteristic()
        << '\n';
    out << "reduction_polynomial=";

    const auto& poly = field->get_reduction_polynomial();

    for (size_t i = 0; i < poly.size(); ++i)
    {
      if (i)
        out << ' ';

      out << poly[i];
    }

    out << "\n\n";

    for (const auto& code : codes)
    {
      out << code;
      out << "\n---\n";
    }
  } // end for
}


/* ************************************************************************* */
auto ConstructedCodesTable::load_queue(int k, queue<LCode> &extended_code, 
    shared_ptr<const Field> field) 
  -> void
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
    const vector<ConstructedCodesTable> &results, 
    int k0) 
  -> void
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
auto ConstructedCodesTable::split_by_weight_enumerator(
    int k0,
    int nb_thread)
  -> vector<vector<LCode>>
{
  unordered_map<string, vector<LCode>> groups;

  for (const auto& [nk, codes] : table_)
  {
    auto [n, k] = nk;

    if (k != k0)
      continue;


    for (const auto& code : codes)
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

  for (auto& [key, group] : groups)
  {
    auto it = min_element(loads.begin(), loads.end());
    size_t thread_id = distance(loads.begin(), it);

    jobs[thread_id].insert(
        jobs[thread_id].end(),
        group.begin(),
        group.end());

    loads[thread_id] += group.size();
  }

  return jobs;
}


/* ************************************************************************* */
struct ParsedFieldInfo
{
  int characteristic;
  vector<Fint> reduction_polynomial;
};

auto parse_field_header(ifstream& in) -> ParsedFieldInfo
{
  ParsedFieldInfo result;

  string line;

  if (!getline(in, line))
    throw runtime_error("Missing characteristic");

  const string characteristic_prefix = "characteristic=";

  if (!line.starts_with(characteristic_prefix))
    throw runtime_error("Invalid characteristic line");

  result.characteristic =
      stoi(line.substr(characteristic_prefix.size()));

  if (!getline(in, line))
    throw runtime_error("Missing reduction polynomial");

  const string polynomial_prefix = "reduction_polynomial=";

  if (!line.starts_with(polynomial_prefix))
    throw runtime_error("Invalid reduction polynomial line");

  istringstream iss(
      line.substr(polynomial_prefix.size()));

  Fint coeff;

  while (iss >> coeff)
    result.reduction_polynomial.push_back(coeff);

  getline(in, line);

  return result;
}


auto ConstructedCodesTable::load(
    shared_ptr<const Field> field,
    int upper_bound_n,
    const filesystem::path& directory) -> void
{
  namespace fs = filesystem;

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

    if (n > upper_bound_n)
      continue;

    ifstream in(entry.path());

    if (!in)
      throw runtime_error("Cannot open " + entry.path().string());

    ParsedFieldInfo field_info =
        parse_field_header(in);

    string line;
    string current;

    while (getline(in, line))
    {
      if (line == "---")
      {
        if (!current.empty())
        {
          table_[{n, k}].insert(LCode::from_canonical_form(current, field));
          current.clear();
        }
      }
      else
      {
        current += line;
        current += '\n';
      }
    }

    if (!current.empty())
      table_[{n, k}].insert(LCode::from_canonical_form(current, field));

    /*
    auto field = Field(field_info.characteristic, 
        field_info.reduction_polynomial.size(),
        field_info.reduction_polynomial);
    */
  }
}
