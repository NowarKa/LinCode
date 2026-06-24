#include "field.hh"
#include "field_element.hh"
#include "projective_space.hh"
#include "hyperplane.hh"

#include <memory>

static void test_projective_point_equality(const shared_ptr<Field> GF) 
{
  cout << "Testing basic equalities...\n";
  auto v1 = {GF->get_element(1), GF->get_element(0)};
  auto v2 = {GF->get_element(1), GF->get_element(0)};

  ProjectivePoint p1(v1);
  ProjectivePoint p2(v2);

  assert(p1 == p2);

  cout << "Basic equalities test passed\n";
}


static void test_projective_belonging(const shared_ptr<Field> GF) 
{
  cout << "Testing belonging properties...\n";

  FieldVector v1 = {GF->get_element(1), GF->get_element(1)};
  FieldVector v2 = {GF->get_element(1), GF->get_element(1)};

  ProjectivePoint p(v1);

  assert(p.contains(v2));

  cout << "Belonging properties test passed\n";
}


static void test_hyperplane_contains(shared_ptr<const Field> GF) 
{
  cout << "Testing belonging properties...\n";

  const ProjectiveSpace ps(2, GF);

  FieldVector u = {GF->get_element(1), GF->get_element(1)};
  Hyperplane H(u);

  FieldVector v = {GF->get_element(1), GF->get_element(1)};
  ProjectivePoint p(v);

  assert(H.contains(p) == true);

  cout << "Belonging properties test passed\n";
}


static void test_hyperplane_not_contains(shared_ptr<const Field> GF) 
{
  cout << "Testing hyperplane properties...\n";

  ProjectiveSpace ps(2, GF);

  FieldVector u = {GF->get_element(1), GF->get_element(1)};
  Hyperplane H(u);

  FieldVector v = {GF->get_element(1), GF->get_element(0)};
  ProjectivePoint p(v);

  assert(H.contains(p) == false);

  cout << "Hyperplane properties test passed\n";
}


static void test_pg14_size(shared_ptr<const Field> GF) 
{
  cout << "Testing PG(1,4)'s size...\n";

  ProjectiveSpace ps(2, GF);
  auto pts = ps.get_all_points();

  assert(pts.size() == 5);

  cout << "PG(1,4)'s size test passed\n";
}


static void test_no_zero_vector(shared_ptr<const Field> GF) 
{
  ProjectiveSpace ps(3, GF);
  auto pts = ps.get_all_points();

  for (auto& p : pts) 
  {
    bool all_zero = true;

    for (auto& x : p.get_coordinates()) 
      if (!(x == GF->get_element(0))) 
      {
        all_zero = false;
        break;
      }

    assert(!all_zero);
  }
}


static void test_concatenate(const shared_ptr<Field> GF) 
{
  FieldVector v = {GF->get_element(1), GF->get_element(1)};
  ProjectivePoint p(v);

  auto p2 = p.concatenate(GF->get_element(1));

  assert(p2.get_coordinates().size() == 3);
  assert(p2.get_coordinate(2) == GF->get_element(1));
}

int main() {
  cout << "================= PROJECTIVE SPACE TESTS START =================\n";

  Field GF = Field(2, 2, {1,1,1});
  auto GF_ptr = make_shared<Field>(GF);

  test_projective_point_equality(GF_ptr);
  test_projective_belonging(GF_ptr);
  test_hyperplane_contains(GF_ptr);
  test_hyperplane_not_contains(GF_ptr);
  test_pg14_size(GF_ptr);
  test_no_zero_vector(GF_ptr);
  test_concatenate(GF_ptr);

  cout << "====================== ALL TESTS PASSED ========================\n";
}
