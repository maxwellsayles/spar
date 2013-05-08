#include <iostream>

#include <gmp.h>

#include "spar.h"

using namespace std;

int main(int argc, char** argv) {
  mpz_t N;
  mpz_t d;
  mpz_init_set_str(N, "4722366490463147393501", 10);
  mpz_init(d);

  Spar spar;
  spar.factor(d, N);

  gmp_printf("d = %Zd\n", d);

  mpz_clear(N);
  mpz_clear(d);
  return 0;
}
