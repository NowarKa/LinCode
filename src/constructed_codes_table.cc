#include "constructed_codes_table.hh"
#include "linear_code.hh"

#include <fstream>
#include <iomanip>


/* ************************************************************************* */
auto ConstructedCodesTable::insert_code(LCode& code) -> void
{
  size_t n = code.get_nb_columns();
  size_t k = code.get_nb_rows();

  table_[{n, k}].insert(code.canonical_form());

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

  for (const auto& s : table_.at({n, k})) 
    if (s == cn_code)
      return true;

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
auto ConstructedCodesTable::save( const std::filesystem::path& directory) 
  const -> void
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

    for (const auto& code : codes)
    {
      out << code;
      out << "\n---\n";
    }
  } // end for
}


/* ************************************************************************* */
auto ConstructedCodesTable::load(int upper_bound_n, 
    const filesystem::path& directory) -> void
{
  namespace fs = std::filesystem;

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
    string line;
    string current;

    while (getline(in, line))
    {
      if (line == "---")
      {
        if (!current.empty())
        {
          table_[{n,k}].insert(current);
          current.clear();
        }
      }
      else
      {
        current += line;
        current += '\n';
      }
    } // end while

    if (!current.empty())
      table_[{n,k}].insert(current);
  } // end for

  // return result;
}
