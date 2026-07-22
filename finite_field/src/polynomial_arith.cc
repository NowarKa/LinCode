#include <algorithm>
#include <iostream>
#include <tuple>
#include <utility>
#include <vector>

#include "polynomial_arith.hh"

using namespace std;

/* ************************************************************************* */
auto polynomial_add(const vector<Fint> &A, const vector<Fint> &B, Fint modulus)
    -> vector<Fint> {
  int min_degree = min(A.size(), B.size());
  int sum_degree = max(A.size(), B.size());
  vector<Fint> sum = vector<Fint>(sum_degree, 0);
  for (int i = 0; i < min_degree; i++)
    sum[i] = add_modular(A[i], B[i], modulus);

  for (int i = min_degree; i < sum_degree; i++)
    if (sum_degree == A.size())
      sum[i] = A[i];
    else
      sum[i] = B[i];
  return sum;
}

/* ************************************************************************* */
auto polynomial_subtract(const vector<Fint> &A, const vector<Fint> &B,
                         Fint modulus) -> vector<Fint> {
  int min_degree = min(A.size(), B.size());
  int sub_degree = max(A.size(), B.size());
  vector<Fint> sub = vector<Fint>(sub_degree, 0);
  for (int i = 0; i < min_degree; i++)
    sub[i] = subtract_modular(A[i], B[i], modulus);

  for (int i = min_degree; i < sub_degree; i++)
    if (sub_degree == A.size())
      sub[i] = A[i];
    else
      sub[i] = ((B[i] > 0) ? (modulus - B[i]) : 0);
  return sub;
}

/* ************************************************************************* */
auto polynomial_multiply(const vector<Fint> &A, const vector<Fint> &B,
                         Fint modulus) -> vector<Fint> {
  int Adeg = polynomial_degree(A);
  int Bdeg = polynomial_degree(B);
  if (Adeg == -1 || Bdeg == -1)
    return vector<Fint>(A.size(), 0);
  vector<Fint> P = vector<Fint>(Adeg + Bdeg + 1, 0);
  for (int i = 0; i <= Adeg + Bdeg; i++) {
    for (int j = 0; j <= i; j++) {
      if ((i - j <= Bdeg) && (j <= Adeg)) {
        Fint temp = multiply_modular(A[j], B[i - j], modulus);
        P[i] = add_modular(P[i], temp, modulus);
      }
    }
  }
  return P;
}

/* ************************************************************************* */
auto polynomial_divide(const vector<Fint> &A, const vector<Fint> &B,
                       Fint modulus) -> pair<vector<Fint>, vector<Fint>> {
  int Adeg = polynomial_degree(A);
  int Bdeg = polynomial_degree(B);

  if (Adeg < Bdeg)
    return {{0}, A};

  auto Q = vector<Fint>(Adeg - Bdeg + 1, 0);
  auto R = vector<Fint>(Bdeg, 0);

  if (Bdeg == -1)
    throw ErrorDivideByZero;

  if (Adeg == -1)
    return {Q, R};

  // allocate array to hold result of B multplied by q
  auto qB = vector<Fint>(Adeg + 1, 0);

  // copy A to R
  int Rdeg = Adeg;
  R = A;

  while (Rdeg >= Bdeg) {
    Fint q = divide_modular(R[Rdeg], B[Bdeg],
                            modulus); // divide leading coefficients
    int diff = Rdeg - Bdeg;
    // update Q
    Q[diff] = q;

    // do this: qB = q * B
    for (int i = Bdeg; i >= 0; i--)
      qB[i + diff] = multiply_modular(B[i], q, modulus);

    for (int i = 0; i < diff; i++)
      qB[i] = 0;

    // now do this: R = R - q*B
    for (int i = 0; i <= Rdeg; i++)
      R[i] = subtract_modular(R[i], qB[i], modulus);

    // determine degree of R
    Rdeg = polynomial_degree(R);

  } // end while
  return {Q, R};
}

/* ************************************************************************* */
auto polynomial_extended_euclid(const vector<Fint> &A, const vector<Fint> &B,
                                Fint modulus)
    -> tuple<vector<Fint>, vector<Fint>, vector<Fint>> {
  vector<Fint> r0 = A;
  vector<Fint> r1 = B;

  vector<Fint> x0 = {1};
  vector<Fint> x1 = {0};

  vector<Fint> y0 = {0};
  vector<Fint> y1 = {1};

  while (polynomial_degree(r1) != -1) {
    auto [q, r2] = polynomial_divide(r0, r1, modulus);

    auto x2 =
        polynomial_subtract(x0, polynomial_multiply(q, x1, modulus), modulus);
    auto y2 =
        polynomial_subtract(y0, polynomial_multiply(q, y1, modulus), modulus);

    r0 = r1;
    r1 = r2;

    x0 = x1;
    x1 = x2;

    y0 = y1;
    y1 = y2;
  }

  return {x0, y0, r0};
}

/* ************************************************************************* */
void polynomial_print(const vector<Fint> &A) {
  int Adeg = polynomial_degree(A);
  if (Adeg == -1)
    cout << "( 0 )";
  else {
    cout << "( ";
    for (int i = Adeg; i >= 0; i--) {
      if (i > 0)
        cout << A[i] << " x^" << i << " + ";
      else
        cout << A[i] << " )";
    }
  }
}
