/**
 * @file spar.h
 * Interface for the SPAR factoring algorithm.
 */
#pragma once

#include <gmp.h>
#include <stdint.h>

extern "C" {
#include "liboptarith/closest_23.h"
#include "liboptarith/group.h"
#include "liboptarith/group_pow.h"
#include "libqform/mpz_qform.h"
#include "libqform/qform_group.h"
#include "libqform/s64_qform.h"
#include "libqform/s128_qform.h"
#include "libsspar/open_addr_hash.h"
}

const int spar_sylow_form_count = 16;
const int spar_random_form_count = 16;
const int spar_hist_form_count = 7;

/**
 * The structure that maintains the data for the SuperSPAR algorithm.
 */
class Spar {
 public:
  Spar();
  ~Spar();
  bool factor(mpz_t d, const mpz_t N);
  
 private:
  void setup_exponentiation_stage();
  void setup_search_stage();
  void setup_discriminant(const unsigned int k);
  int next_primeform(qform_t* form, int prime_index);
  void factor_using_group();
  void exponentiation_stage();
  uint32_t search_stage();
  bool in_history(uint32_t* out_exp, const qform_t* form);
  void record_in_history(int steps_taken,
			 const qform_t* form, const uint32_t exp);
  void repeatedly_square(qform_t* form);
  void add_to_sylow(const qform_t* ambig_form, const qform_t* form);
  void test_sylow();
  
  mpz_t N;   // Input to be factored.
  mpz_t d;   // divisor of N
  double r;  // sqrt(log(N) / log(log(N)))
  double t;  // N^(1/(2*r))
  int w;     // The number of primes used for the exponentiation stage.

  mpz_t D;   // Discriminant

  // Active group elements.
  // These are polymorphic based on discriminant.
  qform_group_t* qform_group;
  group_t* group;
  qform_t* initform;
  qform_t* search;
  qform_t* current;
  qform_t* random_forms[spar_random_form_count];
  uint32_t random_exps[spar_random_form_count];
  qform_t* hist_forms[spar_hist_form_count];
  uint32_t hist_exps[spar_hist_form_count];
  int hist_size;
  double next_history;
  qform_t* sylow_forms[spar_sylow_form_count];
  qform_t* ambig_forms[spar_sylow_form_count];
  int sylow_size;

  // Concrete group elements for each discriminant group.
  s64_qform_group_t qform_group_s64;
  s64_qform_t initform_s64;
  s64_qform_t search_s64;
  s64_qform_t current_s64;
  s64_qform_t ambiguous_forms_s64;
  s64_qform_t random_forms_s64[spar_random_form_count];
  s64_qform_t hist_forms_s64[spar_hist_form_count];
  s64_qform_t sylow_forms_s64[spar_sylow_form_count];
  s64_qform_t ambig_forms_s64[spar_sylow_form_count];
  s128_qform_group_t qform_group_s128;
  s128_qform_t initform_s128;
  s128_qform_t search_s128;
  s128_qform_t current_s128;
  s128_qform_t ambiguous_form_s128;
  s128_qform_t random_forms_s128[spar_random_form_count];
  s128_qform_t hist_forms_s128[spar_hist_form_count];
  s128_qform_t sylow_forms_s128[spar_sylow_form_count];
  s128_qform_t ambig_forms_s128[spar_sylow_form_count];
  mpz_qform_group_t qform_group_mpz;
  mpz_qform_t initform_mpz;
  mpz_qform_t search_mpz;
  mpz_qform_t current_mpz;
  mpz_qform_t ambiguous_form_mpz;
  mpz_qform_t random_forms_mpz[spar_random_form_count];
  mpz_qform_t hist_forms_mpz[spar_hist_form_count];
  mpz_qform_t sylow_forms_mpz[spar_sylow_form_count];
  mpz_qform_t ambig_forms_mpz[spar_sylow_form_count];

  // Powering
  group_pow_t* pow;  // Polymorphic
  group_pow_t pow_s64;
  group_pow_t pow_s128;
  group_pow_t pow_mpz;

  // The order of the ideal once found.
  uint32_t order;

  int logh;

  // Primorial index used for big exponentiation.
  const factored_two_three_term16_t* primorial_terms;
  int primorial_term_count;
  const factored_two_three_term16_t* primorial_terms_s64;
  int primorial_term_count_s64;
  const factored_two_three_term16_t* primorial_terms_s128;
  int primorial_term_count_s128;

  // The number of coprime baby-steps \le step_bound not including the
  // identity.
  int step_count;
};

