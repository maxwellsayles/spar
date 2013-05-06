/**
 * @file main.c
 * 
 * This program generates C code for 2,3 representations of
 * the first N power primorials that are the product of all prime powers
 * ${p_i}^{e_i}$ for $e_i = \floor{\log_{p_i} {p_t}^2}$ and $p_i \le p_t$.
 */

#include <gmp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "liboptarith/closest_23.h"
#include "liboptarith/math_mpz.h"
#include "liboptarith/primes.h"
#include "liboptarith/primorial.h"
#include "libqform/mpz_qform.h"
#include "libqform/s64_qform.h"
#include "libqform/s128_qform.h"

// Largest prime encountered
#define N 1024

char* header = \
"#include <stdint.h>\n" \
"\n" \
"/**\n" \
" * Representation is given as an array 'A' such that\n" \
" * N = ((((A[0].b*A[0].a) \\pm A[1].b)*A[1].a) \\pm A[2].b)*A[2].a ....\n" \
" *\n" \
" * where A[i].b is added if the high bit of 'b' is clear,\n" \
" * and subtracted if the high bit of 'b' is set.\n" \
" *\n" \
" * We use the high bit so that if 'b' is 0, we can still distinguish between\n" \
" * a 'b' that should be added and a 'b' that should be subtracted.\n" \
" */\n" \
"typedef struct {\n" \
"    uint16_t a;\n" \
"    uint16_t b;\n" \
"} factored_two_three_term16_t;\n\n";

void write_two_three(const int i,
		     const int k,
		     const group_cost_t* costs,
		     const mpz_t primorial,
		     int* term_counts,
		     FILE* f) {
  int j;
 
  // Compute the factored 2,3 representation for s64
  int term_count;
  factored_two_three_term16_t* terms =
      factored_rep_prune_closest(&term_count,
				 primorial,
				 costs,
				 k);
  term_counts[i] = term_count;

  // Write to file
  fprintf(f,
	  "static const factored_two_three_term16_t primorial%d[%d] = {",
	  i,
	  term_count);
  for (j = 0; j < term_count; j ++) {
    if (j != 0) {
      fprintf(f, ", ");
    }
    fprintf(f, "{%d, %d}", terms[j].a, terms[j].b);
  }
  fprintf(f, "};\n");
  fflush(f);
  
  // Release the 2,3 form
  free(terms);
}

int main(int argc, char** argv) {
  if (argc != 5 && argc != 6) {
    printf("Usage: %s <s64|s128|mpz> <start_primorial> <end_primorial> <k> [<out_file>]\n\n", argv[0]);
    exit(0);
  }

  char* type = argv[1];
  int start = atoi(argv[2]);
  int end = atoi(argv[3]);
  int k = atoi(argv[4]);
  const char* filename = argc == 6 ? argv[5] : 0;

  int term_count[N];

  // Based on type, open the file and set the costs.
  FILE* file;
  const group_cost_t* costs;
  if (strcmp(type, "s64") == 0) {
    if (!filename) filename = "s64_primorials.c";
    costs = &s64_qform_costs;
  } else if (strcmp(type, "s128") == 0) {
    if (!filename) filename = "s128_primorials.c";
    costs = &s128_qform_costs;
  } else if (strcmp(type, "mpz") == 0) {
    if (!filename) filename = "mpz_primorials.c";
    costs = &mpz_qform_costs;
  } else {
    printf("Unrecognized type \"%s\"\n", type);
    exit(0);
  }
  file = fopen(filename, "a");

  // Write headers.
  if (start == 0) {
    fprintf(file, "%s", header);
    fflush(file);
  }

  mpz_t primorial;
  mpz_init(primorial);
  int i;
  for (i = start; i < end; i ++) {
    unsigned int p = prime_list[i - 1];
    mpz_power_primorial(primorial, i, p * p);
    // Remove powers of two from the primorial.
    mpz_tdiv_q_2exp(primorial, primorial, mpz_scan1(primorial, 0));

    // Print the primorial
    gmp_printf("P_%d %Zd\n", i, primorial);
    
    // Compute 2,3 reps and write them out.
    write_two_three(i, k, costs, primorial, term_count, file);
  }
  mpz_clear(primorial);
  fprintf(file, "\n");

  if (end == N) {
    // Write out indexes
    fprintf(file,
	    "const factored_two_three_term16_t* primorial_terms[%d] = {", N);
    for (i = 0; i < N; i ++) {
      if (i != 0) {
	fprintf(file, ", ");
      }
      fprintf(file, "primorial%d", i);
    }
    fprintf(file, "};\n");
    fflush(file);
    
    // Write out sizes
    fprintf(file,
	    "const int primorial_term_counts[%d] = {", N);
    for (i = 0; i < N; i ++) {
      if (i != 0) {
	fprintf(file, ", ");
      }
      fprintf(file, "%d", term_count[i]);
    }
    fprintf(file, "};\n");
    fflush(file);
  }

  // close all the files
  fclose(file);

  return 0;
}


