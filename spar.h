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

/**
 * The structure that maintains the data for the SuperSPAR algorithm.
 */
class Spar {
 public:
  Spar();
  ~Spar();
  bool factor(mpz_t d, const mpz_t N);
  
 private:
  void setup_exponentiation_stage(const mpz_t N, const double t);
  void setup_discriminant(const mpz_t N, const unsigned int k);
  int next_primeform(qform_t* form, int prime_index);
  void factor_using_group();
  void exponentiation_stage();
  void repeatedly_square(qform_t* form);
  
  mpz_t N;  // Input to be factored.
  mpz_t d;  // divisor of N

  mpz_t D;  // Discriminant
  mpz_t t;  // Temporary

  // Active group elements.
  // These are polymorphic based on discriminant.
  qform_group_t* qform_group;
  group_t* group;
  qform_t* initform;
  qform_t* search;
  qform_t* current;
  qform_t* temp_form;

  // Concrete group elements for each discriminant group.
  s64_qform_group_t qform_group_s64;
  s64_qform_t initform_s64;
  s64_qform_t search_s64;
  s64_qform_t giant_s64;
  s64_qform_t current_s64;
  s64_qform_t temp_form_s64;
  s128_qform_group_t qform_group_s128;
  s128_qform_t initform_s128;
  s128_qform_t search_s128;
  s128_qform_t giant_s128;
  s128_qform_t current_s128;
  s128_qform_t temp_form_s128;
  mpz_qform_group_t qform_group_mpz;
  mpz_qform_t initform_mpz;
  mpz_qform_t search_mpz;
  mpz_qform_t giant_mpz;
  mpz_qform_t current_mpz;
  mpz_qform_t temp_form_mpz;

  // Powering
  group_pow_t* pow;  // Polymorphic
  group_pow_t pow_s64;
  group_pow_t pow_s128;
  group_pow_t pow_mpz;

  // The order of the ideal once found.
  // NOTE: This is limited to 32bits as the hashtable used has 32bit keys.
  uint32_t order;
  int found_order;

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

