#include <iostream>

#include <gmp.h>

#include "spar.h"

using namespace std;

int main(int argc, char** argv) {
  if (argc == 1) {
    cout << "Usage: " << argv[0] << " <integer>" << endl;
    return 0;
  }
  mpz_t N;
  mpz_t d;

  mpz_init_set_str(N, argv[1], 10);
  mpz_init(d);

  Spar spar;
  spar.factor(d, N);

  gmp_printf("d = %Zd\n", d);

  mpz_clear(N);
  mpz_clear(d);
  return 0;
}
