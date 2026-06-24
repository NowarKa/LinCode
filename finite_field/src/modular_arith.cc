#include <cassert>
#include <tuple>

#include "modular_arith.hh"

// Computes greatest common divisor of two positive integers
auto extended_euclid(Fint a, Fint b) -> tuple<Fint, Fint, Fint>
{
  assert(a > 0);
  assert(b > 0);
  Fint x, y, d;
  Fint q, r;
  Fint x1, x2, y1, y2;

  if (b == 0)
  {
    d = a;
    x = 1;
    y = 0;
  }
  else
  {
    x1 = 0;
    x2 = 1;
    y1 = 1;
    y2 = 0;

    while (b > 0)
    {
      q = a / b;
      r = a % b;
      a = b;
      b = r;
      x = x2 - q * x1;
      y = y2 - q * y1;
      x2 = x1;
      x1 = x;
      y2 = y1;
      y1 = y;
    }

    d = a;
    x = x2;
    y = y2;
  }

  return make_tuple(x, y, d);
}
