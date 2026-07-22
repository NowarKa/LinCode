#include <cassert>
#include <iostream>
#include <memory>
#include <vector>

#include "field.hh"
#include "field_element.hh"
#include "polynomial_arith.hh"

using namespace std;

static void test_basic_field_properties(const shared_ptr<Field> GF) {
  cout << "Testing basic field properties...\n";

  auto zero = FieldElement(GF);
  zero.set_zero();

  auto one = FieldElement(GF);
  one.set_one();

  assert(zero.is_zero());
  assert(!(one.is_zero()));

  cout << "Basic field properties test passed\n";
}

static void test_additive_group(const shared_ptr<Field> GF) {
  cout << "Testing additive group...\n";

  auto elems = vector<FieldElement>();

  for (uint64_t i = 0; i < GF->get_order(); ++i)
    elems.push_back(GF->get_element(i));

  for (auto &a : elems)
    for (auto &b : elems)
      for (auto &c : elems) {
        assert((a + b) + c == a + (b + c));
        assert(a + b == b + a);

        auto zero = FieldElement(GF);
        zero.set_zero();

        assert(a + zero == a);
        assert(a + (-a) == zero);
      }

  cout << "Additive group test passed\n";
}

static void test_multiplicative_group(const shared_ptr<Field> GF) {
  cout << "Testing multiplicative group...\n";

  auto elems = vector<FieldElement>();

  for (uint64_t i = 0; i < GF->get_order(); ++i)
    elems.push_back(GF->get_element(i));

  auto one = FieldElement(GF);
  one.set_one();

  for (auto &a : elems)
    for (auto &b : elems)
      for (auto &c : elems) {
        assert((a * b) * c == a * (b * c));
        assert(a * b == b * a);
        assert(a * one == a);
      }

  cout << "Multiplicative group test passed\n";
}

static void test_distributivity(const shared_ptr<Field> GF) {
  cout << "Testing distributivity...\n";

  auto elems = vector<FieldElement>();

  for (uint64_t i = 0; i < GF->get_order(); ++i)
    elems.push_back(GF->get_element(i));

  for (auto &a : elems)
    for (auto &b : elems)
      for (auto &c : elems)
        assert(a * (b + c) == a * b + a * c);

  cout << "Distributivity test passed\n";
}

static void test_inverses(const shared_ptr<Field> GF) {
  cout << "Testing inverses...\n";

  auto one = FieldElement(GF);
  one.set_one();

  for (uint64_t i = 1; i < GF->get_order(); ++i) {
    auto a = GF->get_element(i);

    auto inv = mult_inverse(a);

    assert(a * inv == one);
  }

  cout << "Inverses test passed\n";
}

static void test_exponentiation(const shared_ptr<Field> GF) {
  cout << "Testing exponentiation...\n";

  auto one = FieldElement(GF);
  one.set_one();

  for (uint64_t i = 0; i < GF->get_order(); ++i) {
    auto a = GF->get_element(i);

    assert(field_exp(a, 0) == one);
    assert(field_exp(a, 1) == a);
    assert(field_exp(a, 2) == a * a);
  }

  cout << "Exponentiation test passed\n";
}

static void test_polynomials() {
  cout << "Testing polynomial arithmetic...\n";

  Fint mod = 2;

  vector<Fint> A = {1, 1};    // x + 1
  vector<Fint> B = {1, 1, 1}; // x^2 + x + 1

  auto P = polynomial_multiply(A, B, mod);

  // (x+1)(x^2+x+1) = x^3 + 1 over F2
  vector<Fint> expected_result = {1, 0, 0, 1};

  assert(P == expected_result);

  cout << "Polynomial arithmetic test passed\n";
}

static void test_extended_euclid_polynomials() {
  cout << "Testing extended Euclid (polynomials)...\n";

  Fint mod = 2;

  vector<Fint> A = {1, 0, 1}; // x^2 + 1
  vector<Fint> B = {1, 1};    // x + 1

  auto [X, Y, G] = polynomial_extended_euclid(A, B, mod);

  auto left = polynomial_add(polynomial_multiply(X, A, mod),
                             polynomial_multiply(Y, B, mod), mod);

  assert(left == G);

  cout << "Extended Euclid test passed\n";
}

int main() {
  cout << "========== FIELD TESTS START ==========\n";

  auto GF = std::make_shared<Field>(2, 3, vector<Fint>{1, 0, 1, 1}); // GF(2^3)

  test_basic_field_properties(GF);
  test_additive_group(GF);
  test_multiplicative_group(GF);
  test_distributivity(GF);
  test_inverses(GF);
  test_exponentiation(GF);
  test_polynomials();
  test_extended_euclid_polynomials();

  cout << "========== ALL TESTS PASSED ==========\n";
}
