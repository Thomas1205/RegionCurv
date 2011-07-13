/*** written by Thomas Schoenemann as a private person without employment, November 2009 ***/

#include "combinatoric.hh"

uint fac(uint n) {

  uint r = 1;
  for (uint k=2; k <= n; k++)
    r *= k;

  return r;
}

long double ldfac(uint n) {

  long double r = 1.0;

  for (uint k=2; k <= n; k++)
    r *= k;

  return r;  
}

uint choose(uint n, uint k) {

  assert(k <= n);

  if (k == n || k == 0)
    return 1;

  uint r = n;
  for (uint i=n-1; i > (n-k); i--)
    r *= i;

  for (uint i=2; i <= k; i++)
    r /= i;

  return r;
}

long double ldchoose(uint n, uint k) {

  assert(k <= n);

  if (k == n || k == 0)
    return 1;

  long double r = n;
  for (uint i=n-1; i > (n-k); i--)
    r *= i;

  for (uint i=2; i <= k; i++)
    r /= i;

  return r;
}
