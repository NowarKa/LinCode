#include "linear_code.hh"
#include "field_element.hh"
#include "field_vector.hh"
#include "hyperplane.hh"
#include "projective_space.hh"

#include <cstdint>
#include <iostream>
#include <string>
#include <iomanip>
#include <algorithm>
#include <vector>
#include <fstream>
#include <sstream>
#include <filesystem>
#include <cstdlib>


/* ************************************************************************* */
auto transpose(const vector<vector<FieldElement>>& matrix) 
  -> vector<vector<FieldElement>>
{
  if (matrix.empty())
    return {};

  size_t nb_rows = matrix.size();
  size_t nb_columns = matrix[0].size();

  vector<vector<FieldElement>> t_matrix;

  for (size_t j = 0; j < nb_columns; j++)
  {
    vector<FieldElement> line;

    for (size_t i = 0; i < nb_rows; i++)
      line.push_back(matrix[i][j]);

    t_matrix.push_back(line);
  }

  return t_matrix;
}


/* ************************************************************************* */
LCode::LCode(const vector<vector<FieldElement>>& rows) : 
  rows_(rows), 
  columns_(transpose(rows))
{}


/* ************************************************************************* */
auto construct_from_columns(const vector<vector<FieldElement>>& columns) 
  -> LCode
{
  return LCode(transpose(columns));
}


/* ************************************************************************* */
auto LCode::get_minimum_column_multiplicity() -> int
{
  if (minimum_column_multiplicity_ != 0)
    return minimum_column_multiplicity_;

  uint32_t min_mult = get_nb_columns();

  for (auto& [p, m] : to_multiset())
    min_mult = min(min_mult, m);

  minimum_column_multiplicity_ = min_mult;
  return min_mult;
}


/* ************************************************************************* */
auto LCode::nb_columns_belong_to(const ProjectivePoint& p) const -> uint32_t
{
  uint32_t nb = 0;
  for (auto& v : columns_)
    if (p.contains(v))
      nb++;
  return nb;
}


/* ************************************************************************* */
auto LCode::to_multiset() const -> unordered_map<ProjectivePoint, uint32_t>
{
  unordered_map<ProjectivePoint, uint32_t> multiset;
  for (const auto& column : columns_)
  {
    auto p = ProjectivePoint(column);
    if (multiset.find(p) == multiset.end())
      multiset.insert({p,1});
    else
      multiset.at(p)++;
  }
  return multiset;
}


/* ************************************************************************* */
auto LCode::weight(const vector<Hyperplane>& hyperplanes) const -> uint32_t
{
  int d = get_nb_rows();
  int n = get_nb_columns();

  auto M = to_multiset();

  for (const auto& H : hyperplanes)
  {
    int MH = 0;

    for (auto const& [P, mult] : M)
      if (H.contains(P))
        MH += mult;

    d = min(d, n - MH);
  }

  return d;
}


/* ************************************************************************* */
auto LCode::encode_column_vector(const vector<FieldElement>& u) const
  -> vector<FieldElement>
{
  if (rows_.empty()) 
    return {};

  size_t n = u.size();
  size_t m = rows_[0].size();

  FieldElement zero = rows_[0][0].get_field()->get_element(0);

  if (rows_.size() != n)
    throw invalid_argument("Incompatible dimension");

  /*
  for (const auto& row : G)
    if (row.size() != m)
      throw std::invalid_argument("The matrix is not rectangular");
  */

  vector<FieldElement> result(m, zero);

  for (size_t j = 0; j < m; ++j)
    for (size_t i = 0; i < n; ++i)
      result[j] = result[j] + (u[i] * rows_[i][j]);

  return result;
}


/* ************************************************************************* */
auto hamming_weight(const FieldVector& u) -> int
{
  int d = 0;

  if (u.empty())
    return d;

  FieldElement zero = u[0].get_field()->get_element(0);

  for (auto x : u)
    if (x != zero)
      d++;
  return d;
}


/* ************************************************************************* */
auto LCode::weight() const -> uint32_t
{
  auto ps = ProjectiveSpace(get_nb_rows(), rows_[0][0].get_field());
  auto points = ps.get_all_points();

  int d = get_nb_columns();

  for (const auto& p : points)
  {
    auto encoded = encode_column_vector(p.get_coordinates());
    int h = hamming_weight(encoded);
    
    if (h > 0)
      d = min(d, h);

  }
  return d;
}


/* ************************************************************************* */
ostream& operator<<(ostream& output, const LCode& right)
{
  int n = right.get_nb_rows();
  int m = right.get_nb_columns();

  // Width of each colomn
  vector<int> width(m, 0);

  for (int j = 0; j < m; ++j)
    for (int i = 0; i < n; ++i)
      width[j] = max(width[j], (int)to_string(right(i,j).index()).size());

  auto print_line = [&](int i, const string& left, const string& r)
  {
    output << left << " ";

    for (int j = 0; j < m; ++j)
    {
      output << setw(width[j]) << right(i, j).index();
      if (j + 1 < m) cout << " ";
    }

    output << " " << r << "\n";
  };

  for (int i = 0; i < n; ++i)
    if (i == 0)
      print_line(i, "⎛", "⎞");
    else if (i == n - 1)
      print_line(i, "⎝", "⎠");
    else
      print_line(i, "⎜", "⎟");

  return output;
}


/* ************************************************************************* */
auto LCode::canonical_form(const string& sage_binary) -> const string
{
  namespace fs = filesystem;

  if (columns_.empty())
    return "";

  if (canonical_form_ != "")
    return canonical_form_;

  auto q = columns_[0][0].get_field()->get_order();

  auto tmpdir = fs::temp_directory_path();

  auto input_file  = tmpdir / "codecan_input.txt";
  auto output_file = tmpdir / "codecan_output.txt";
  auto script_file = tmpdir / "codecan_script.sage";


  // Exporting matrix.
  {
    ofstream out(input_file);

    const size_t k = get_nb_rows();
    const size_t n = get_nb_columns();

    out << q << " " << k << " " << n << "\n";

    for (size_t i = 0; i < k; ++i)
    {
      for (size_t j = 0; j < n; ++j)
      {
        if (j)
          out << " ";
        out << (*this)(i,j).index();
      }

      out << "\n";
    }
  }


  // Writing Sage script.
  {
    ofstream out(script_file);

    out << R"(
from sage.all import *

with open(r')"
    << input_file.string()
    << R"(', 'r') as f:

    q,k,n = map(int, f.readline().split())

    F = GF(q)

    rows = []

    for _ in range(k):
        rows.append([F(int(x)) for x in f.readline().split()])

G = Matrix(F, rows)

C = LinearCode(G)

Can, T = C.canonical_representative(
    equivalence='linear'
)

M = Can.generator_matrix()

with open(r')"
    << output_file.string()
    << R"(', 'w') as f:

    f.write(str(M.nrows()) + " "
            + str(M.ncols()) + "\n")

    for row in M.rows():
        f.write(
            " ".join(
                str(x)
                for x in row
            )
        )
        f.write("\n")
)";
  }


  // Executing Sage.
  string command =
    sage_binary + " " + script_file.string();

  int ret = system(command.c_str());

  if (ret != 0)
    throw std::runtime_error(
        "Sage failed while computing canonical form."
        );

  // Parsing output.
  ifstream in(output_file);
  ostringstream result;
  result << in.rdbuf();
  canonical_form_ = result.str();

  return result.str();
}
