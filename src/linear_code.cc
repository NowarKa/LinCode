#include "linear_code.hh"
#include "field_element.hh"
#include "field_vector.hh"
#include "hyperplane.hh"
#include "projective_space.hh"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

/* ************************************************************************* */
auto transpose(const vector<vector<FieldElement>> &matrix)
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
LCode::LCode(const vector<vector<FieldElement>> &rows)
    : rows_(rows), columns_(transpose(rows)), weight_enumerator_({})
{
}

/* ************************************************************************* */
LCode::LCode(const unordered_map<ProjectivePoint, uint32_t> &multiset)
    : multiset_(multiset)
{
  for (const auto &[p, m] : multiset)
    for (size_t i = 0; i < m; i++)
      columns_.push_back(p.get_coordinates());

  rows_ = transpose(columns_);
}

/* ************************************************************************* */
LCode::LCode(
    shared_ptr<const unordered_map<Hyperplane, vector<int>>> &containing_map,
    const vector<int> &solution, shared_ptr<const ProjectiveSpace> ps)
{
  auto points = ps->get_all_points();
  auto np_kp1 = points.size();
  unordered_map<ProjectivePoint, uint32_t> multiset;

  for (size_t i = 0; i < np_kp1; i++)
    if (solution[i] > 0)
    {
      multiset.insert({points[i].get_coordinates(), solution[i]});
      for (int j = 0; j < solution[i]; j++)
        columns_.push_back(points[i].get_coordinates());
    }

  auto n = columns_.size();

  rows_ = transpose(columns_);
  multiset_ = multiset;
  weight_enumerator_ = vector<int>(n + 1, 0);

  for (const auto &[h, v] : *containing_map)
  {
    int sum = 0;

    for (auto i : v)
      sum += solution[i];

    weight_enumerator_[n - sum]++;

    if (n - sum > minimum_distance_)
      minimum_distance_ = n - sum;
  }
}

/* ************************************************************************* */
auto LCode::construct_from_columns(const vector<vector<FieldElement>> &columns)
    -> LCode
{
  return LCode(transpose(columns));
}

/* ************************************************************************* */
auto LCode::remove_projective_point(const ProjectivePoint &p) const -> LCode
{
  vector<vector<FieldElement>> new_columns;

  for (const auto &column : columns_)
    if (ProjectivePoint(column) != p)
      new_columns.push_back(column);

  return LCode::construct_from_columns(new_columns);
}

/* ************************************************************************* */
auto LCode::get_minimum_column_multiplicity() -> int
{
  if (minimum_column_multiplicity_ != -1)
    return minimum_column_multiplicity_;

  uint32_t min_mult = get_nb_columns();

  for (auto &[p, m] : to_multiset())
    min_mult = min(min_mult, m);

  minimum_column_multiplicity_ = min_mult;
  return min_mult;
}

/* ************************************************************************* */
auto LCode::nb_columns_belong_to(const ProjectivePoint &p) const -> uint32_t
{
  uint32_t nb = 0;
  for (auto &v : columns_)
    if (p.contains(v))
      nb++;
  return nb;
}

/* ************************************************************************* */
auto LCode::to_multiset() const -> unordered_map<ProjectivePoint, uint32_t>
{
  if (!multiset_.empty())
    return multiset_;

  unordered_map<ProjectivePoint, uint32_t> multiset;
  for (const auto &column : columns_)
  {
    auto p = ProjectivePoint(column);
    if (multiset.find(p) == multiset.end())
      multiset.insert({p, 1});
    else
      multiset.at(p)++;
  }

  multiset_ = multiset;

  return multiset;
}

/* ************************************************************************* */
auto LCode::get_weight_enumerator(shared_ptr<const ProjectiveSpace> ps) const
    -> vector<int>
{
  if (!weight_enumerator_.empty())
    return weight_enumerator_;

  weight_enumerator_ = vector<int>(get_nb_columns() + 1, 0);

  vector<ProjectivePoint> points;

  if (ps != nullptr)
  {
    points = ps->get_all_points();
  }

  else
  {
    auto psn = ProjectiveSpace(get_nb_rows(), rows_[0][0].get_field());
    points = psn.get_all_points();
  }

  for (const auto &p : points)
  {
    auto encoded = encode_column_vector(p.get_coordinates());
    int h = hamming_weight(encoded);

    if (h > 0)
      weight_enumerator_[h]++;
  }

  return weight_enumerator_;
}

/* ************************************************************************* */
auto LCode::weight(const vector<Hyperplane> &hyperplanes) const -> uint32_t
{
  int d = get_nb_rows();
  int n = get_nb_columns();

  auto M = to_multiset();

  for (const auto &H : hyperplanes)
  {
    int MH = 0;

    for (auto const &[P, mult] : M)
      if (H.contains(P))
        MH += mult;

    d = min(d, n - MH);
  }

  return d;
}

/* ************************************************************************* */
auto LCode::encode_column_vector(const vector<FieldElement> &u) const
    -> vector<FieldElement>
{
  if (rows_.empty())
    return {};

  size_t n = u.size();
  size_t m = rows_[0].size();

  FieldElement zero = rows_[0][0].get_field()->get_element(0);

  if (rows_.size() != n)
    throw invalid_argument("Incompatible dimension");

  vector<FieldElement> result(m, zero);

  for (size_t j = 0; j < m; ++j)
    for (size_t i = 0; i < n; ++i)
      result[j] = result[j] + (u[i] * rows_[i][j]);

  return result;
}

/* ************************************************************************* */
auto LCode::minimum_distance(shared_ptr<const ProjectiveSpace> ps) const
    -> uint32_t
{
  if (minimum_distance_ != -1)
    return minimum_distance_;

  auto we = get_weight_enumerator(ps);

  for (size_t i = 1; i < we.size(); i++)
    if (we[i] > 0)
      return i;

  return 0;
}

/* ************************************************************************* */
auto LCode::is_codeword_weight_divisible(int delta) -> bool
{
  if (delta <= 1)
    return true;

  auto we = get_weight_enumerator();

  for (size_t i = 1; i < we.size(); i++)
    if (we[i] > 0 && i % delta != 0)
      return false;

  return true;
}

/* ************************************************************************* */
ostream &operator<<(ostream &output, const LCode &right)
{
  int n = right.get_nb_rows();
  int m = right.get_nb_columns();

  // Width of each colomn
  vector<int> width(m, 0);

  for (int j = 0; j < m; ++j)
    for (int i = 0; i < n; ++i)
      width[j] = max(width[j], (int)to_string(right(i, j).index()).size());

  auto print_line = [&](int i, const string &left, const string &r)
  {
    output << left << " ";

    for (int j = 0; j < m; ++j)
    {
      output << setw(width[j]) << right(i, j).index();
      if (j + 1 < m)
        cout << " ";
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
bool operator==(const LCode &left, const LCode &right)
{
  return left.rows_ == right.rows_;
}

/* ************************************************************************* */
auto LCode::is_projective() const -> bool
{
  for (const auto &[p, m] : to_multiset())
    if (m > 1)
      return false;

  return true;
}
