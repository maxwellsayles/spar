#include "spar/spar.h"

#include <iomanip>
#include <iostream>

#include <assert.h>
#include <gmp.h>
#include <inttypes.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

extern "C" {
#include "liboptarith/closest_23.h"
#include "liboptarith/group.h"
#include "liboptarith/math32.h"
#include "liboptarith/math_mpz.h"
#include "liboptarith/primes.h"
#include "liboptarith/primorial.h"
#include "liboptarith/square_free.h"
#include "libqform/mpz_qform.h"
#include "libqform/qform_group.h"
#include "libqform/s64_qform.h"
#include "libqform/s128_qform.h"
#include "spar/primorials/s64_theory_opt.h"
#include "spar/primorials/s128_theory_opt.h"
}

using namespace std;

/// Default values for an instance of SuperSPAR.
Spar::Spar() {
  mpz_init(this->N);
  mpz_init(this->d);
  mpz_init(this->D);

  // Initialize groups.
  s64_qform_group_init(&this->qform_group_s64);
  s64_qform_init(&this->qform_group_s64, &this->initform_s64);
  s64_qform_init(&this->qform_group_s64, &this->search_s64);
  s64_qform_init(&this->qform_group_s64, &this->current_s64);
  s64_qform_init(&this->qform_group_s64, &this->temp_form_s64);
  for (int i = 0; i < spar_random_form_count; i++) {
    s64_qform_init(&this->qform_group_s64, &this->random_forms_s64[i]);
  }
  for (int i = 0; i < spar_hist_form_count; i++) {
    s64_qform_init(&this->qform_group_s64, &this->hist_forms_s64[i]);
  }
  group_pow_init(&this->pow_s64, &this->qform_group_s64.desc.group);

  s128_qform_group_init(&this->qform_group_s128);
  s128_qform_init(&this->qform_group_s128, &this->initform_s128);
  s128_qform_init(&this->qform_group_s128, &this->search_s128);
  s128_qform_init(&this->qform_group_s128, &this->current_s128);
  s128_qform_init(&this->qform_group_s128, &this->temp_form_s128);
  for (int i = 0; i < spar_random_form_count; i++) {
    s128_qform_init(&this->qform_group_s128, &this->random_forms_s128[i]);
  }
  for (int i = 0; i < spar_hist_form_count; i++) {
    s128_qform_init(&this->qform_group_s128, &this->hist_forms_s128[i]);
  }
  group_pow_init(&this->pow_s128, &this->qform_group_s128.desc.group);

  mpz_qform_group_init(&this->qform_group_mpz);
  mpz_qform_init(&this->qform_group_mpz, &this->initform_mpz);
  mpz_qform_init(&this->qform_group_mpz, &this->search_mpz);
  mpz_qform_init(&this->qform_group_mpz, &this->current_mpz);
  mpz_qform_init(&this->qform_group_mpz, &this->temp_form_mpz);
  for (int i = 0; i < spar_random_form_count; i++) {
    mpz_qform_init(&this->qform_group_mpz, &this->random_forms_mpz[i]);
  }
  for (int i = 0; i < spar_hist_form_count; i++) {
    mpz_qform_init(&this->qform_group_mpz, &this->hist_forms_mpz[i]);
  }
  group_pow_init(&this->pow_mpz, &this->qform_group_mpz.desc.group);
    
  // Primorial.
  this->primorial_terms = nullptr;
  this->primorial_term_count = 0;
}

/// Release SSPAR resources.
Spar::~Spar() {
  group_pow_clear(&this->pow_s64);
  s64_qform_clear(&this->qform_group_s64, &this->initform_s64);
  s64_qform_clear(&this->qform_group_s64, &this->search_s64);
  s64_qform_clear(&this->qform_group_s64, &this->current_s64);
  s64_qform_clear(&this->qform_group_s64, &this->temp_form_s64);
  for (int i = 0; i < spar_random_form_count; i++) {
    s64_qform_clear(&this->qform_group_s64, &this->random_forms_s64[i]);
  }
  for (int i = 0; i < spar_hist_form_count; i++) {
    s64_qform_clear(&this->qform_group_s64, &this->hist_forms_s64[i]);
  }
  s64_qform_group_clear(&this->qform_group_s64);

  group_pow_clear(&this->pow_s128);
  s128_qform_clear(&this->qform_group_s128, &this->initform_s128);
  s128_qform_clear(&this->qform_group_s128, &this->search_s128);
  s128_qform_clear(&this->qform_group_s128, &this->current_s128);
  s128_qform_clear(&this->qform_group_s128, &this->temp_form_s128);
  for (int i = 0; i < spar_random_form_count; i++) {
    s128_qform_clear(&this->qform_group_s128, &this->random_forms_s128[i]);
  }
  for (int i = 0; i < spar_hist_form_count; i++) {
    s128_qform_clear(&this->qform_group_s128, &this->hist_forms_s128[i]);
  }
  s128_qform_group_clear(&this->qform_group_s128);

  group_pow_clear(&this->pow_mpz);
  mpz_qform_clear(&this->qform_group_mpz, &this->initform_mpz);
  mpz_qform_clear(&this->qform_group_mpz, &this->search_mpz);
  mpz_qform_clear(&this->qform_group_mpz, &this->current_mpz);
  mpz_qform_clear(&this->qform_group_mpz, &this->temp_form_mpz);
  for (int i = 0; i < spar_random_form_count; i++) {
    mpz_qform_clear(&this->qform_group_mpz, &this->random_forms_mpz[i]);
  }
  for (int i = 0; i < spar_hist_form_count; i++) {
    mpz_qform_clear(&this->qform_group_mpz, &this->hist_forms_mpz[i]);
  }
  mpz_qform_group_clear(&this->qform_group_mpz);

  mpz_clear(this->D);
  mpz_clear(this->d);
  mpz_clear(this->N);
}

/**
 * Find the next prime form for prime_index.
 * init_form is the set to the new form.
 */
int Spar::next_primeform(qform_t* form, int prime_index) {
  prime_index = qform_next_primeform(qform_group, form, prime_index);
  assert(prime_index < static_cast<int>(prime_list_count));
  return prime_index;
}

/**
 * Exponentiates initform by pow.
 * Result is in initform.
 */
void Spar::exponentiation_stage() {
  group_pow_factored23(this->pow, this->initform, this->initform,
		       this->primorial_terms, this->primorial_term_count);
}

/**
 * Set all polymorphic variables and the group discriminant.
 * NOTE: Requires this->D and this->primorial_index to be set.
 */
void Spar::setup_discriminant(const unsigned int k) {
  // Set up discriminant for this multiplier
  mpz_mul_ui(this->D, N, k);
  if ((this->D->_mp_d[0] & 3) != 3) {
    // D != 3 (mod 4)
    mpz_mul_2exp(this->D, this->D, 2);
  }
  mpz_neg(this->D, this->D);

  size_t logD = mpz_sizeinbase(this->D, 2);
  if (logD <= s64_qform_group_max_bits) {
    // Use s64 implementations.
    this->qform_group          = &this->qform_group_s64.desc;
    this->initform             = &this->initform_s64;
    this->search               = &this->search_s64;
    this->current              = &this->current_s64;
    this->temp_form            = &this->temp_form_s64;
    for (int i = 0; i < spar_random_form_count; i++) {
      this->random_forms[i] = &this->random_forms_s64[i];
    }
    for (int i = 0; i < spar_hist_form_count; i++) {
      this->hist_forms[i] = &this->hist_forms_s64[i];
    }
    this->pow                  = &this->pow_s64;
    this->primorial_terms      = this->primorial_terms_s64;
    this->primorial_term_count = this->primorial_term_count_s64;
  } else if (logD <= s128_qform_group_max_bits) {
    // Use s128 implementations.
    this->qform_group          = &this->qform_group_s128.desc;
    this->initform             = &this->initform_s128;
    this->search               = &this->search_s128;
    this->current              = &this->current_s128;
    this->temp_form            = &this->temp_form_s128;
    for (int i = 0; i < spar_random_form_count; i++) {
      this->random_forms[i] = &this->random_forms_s128[i];
    }
    for (int i = 0; i < spar_hist_form_count; i++) {
      this->hist_forms[i] = &this->hist_forms_s128[i];
    }
    this->pow                  = &this->pow_s128;
    this->primorial_terms      = this->primorial_terms_s128;
    this->primorial_term_count = this->primorial_term_count_s128;
  } else {
    // Use MPZ implementations.
    this->qform_group          = &this->qform_group_mpz.desc;
    this->initform             = &this->initform_mpz;
    this->search               = &this->search_mpz;
    this->current              = &this->current_mpz;
    this->temp_form            = &this->temp_form_mpz;
    for (int i = 0; i < spar_random_form_count; i++) {
      this->random_forms[i] = &this->random_forms_mpz[i];
    }
    for (int i = 0; i < spar_hist_form_count; i++) {
      this->hist_forms[i] = &this->hist_forms_mpz[i];
    }
    this->pow                  = &this->pow_mpz;
    this->primorial_terms      = this->primorial_terms_s128;
    this->primorial_term_count = this->primorial_term_count_s128;
  }
  this->qform_group->set_discriminant(this->qform_group, this->D);
  this->group = &this->qform_group->group;

  this->logh = (logD + 1) / 2;
  if (this->logh < 1) this->logh = 1;
}

void Spar::setup_exponentiation_stage() {
  w = prime_index_ge(t);
  primorial_terms_s64 = s64_primorial_terms[w-1];
  primorial_term_count_s64 = s64_primorial_term_counts[w-1];
  primorial_terms_s128 = s128_primorial_terms[w-1];
  primorial_term_count_s128 = s128_primorial_term_counts[w-1];

  // We take at most p_t steps.
  step_count = prime_list[w-1];

  /*
  // This is how the precomputed reps were computed.
  const group_cost_t* costs = &mpz_qform_costs;
  int logN = mpz_sizeinbase(N, 2);
  if (logN <= s64_qform_group_max_bits) {
    costs = &s64_qform_costs;
  } else if (logN <= s128_qform_group_max_bits) {
    costs = &s128_qform_costs;
  }
  
  mpz_t E;
  mpz_init(E);

  int w = prime_index_ge(t);
  unsigned int pt = prime_list[w];
  cout << "w=" << w << " prime=" << pt << endl;

  mpz_power_primorial(E, w, pt * pt);
  mpz_tdiv_q_2exp(E, E, mpz_scan1(E, 0));

  // Step count is what it would be if we used binary exponentiation.
  step_count = mpz_sizeinbase(E, 2) + mpz_popcount(E);
  cout << "Step count = " << step_count << endl;

  // compute 2,3 representation of E
  cout << "Computing 2,3 representation for exponentiation stage." << endl;
  primorial_terms = factored_rep_prune_closest(&primorial_term_count,
					       E, costs, 4);

  mpz_clear(E);
  */
}

// Pick 16 random exponents in the range [{p_t}^2, 2{p_t}^2]
// and exponentiate search to each exponent.
// Then perform a Pollard-Brenth random walk.
void Spar::setup_search_stage() {
  int pt2 = prime_list[w-1] * prime_list[w-1];
  for (int q = 0; q < spar_random_form_count; q++) {
    int e = (rand_u32() % pt2) + pt2;
    qform_pow_u32(pow, random_forms[q], search, e);
    random_exps[q] = e;
  }
}

/// d will be non-zero if this found a factor of N
void Spar::repeatedly_square(qform_t* form) {
  int i = 0;
  while (i < logh) {
    if (qform_group->is_ambiguous(qform_group, form))
      break;
    group->square(group, form, form);
    i++;
  }

  if (qform_group->is_ambiguous(qform_group, form)) {
    qform_group->split_ambiguous(qform_group, d, N, form);
  }
}

/// True if form is in the history.  If the form is in the history,
/// the exponent of the matching form is returned in out_exp.
bool Spar::in_history(uint32_t* out_exp, const qform_t* form) {
  for (int i = 0;
       i < hist_size && i < spar_hist_form_count;
       i++) {
    if (group->equal(group, hist_forms[i], form)) {
      *out_exp = hist_exps[i];
      return true;
    }
  }
  return false;
}

// Record to history when the current step is 1.1^8 times larger than
// than the step of the head of history.
void Spar::record_in_history(int steps_taken,
			     const qform_t* form, const uint32_t exp) {
  int m = min(spar_hist_form_count, hist_size) + 1;
  double x = ::pow(1.1, m);
  if (steps_taken >= next_history * x) {
    int i = hist_size % spar_hist_form_count;
    group->set(group, hist_forms[i], form);
    hist_exps[i] = exp;
    hist_size ++;
    if (hist_size >= spar_hist_form_count) {
      next_history *= 1.1;
    } else {
      next_history = 1;
    }
  }
}

/// Return the found order or 1.
uint32_t Spar::search_stage() {
  group->set(group, current, search);
  uint32_t exp = 0;
  hist_size = 0;
  next_history = 0;
  for (int steps_taken = 0;
       steps_taken < step_count;
       steps_taken++) {
    int i = steps_taken % spar_random_form_count;
    group->compose(group, current, current, random_forms[i]);
    exp += random_exps[i];
    uint32_t match_exp;
    if (in_history(&match_exp, current)) {
      return exp - match_exp;
    }
    record_in_history(steps_taken, current, exp);
  }
  return 1;
}

// The class group is determined, now we iterate on various prime ideals.
void Spar::factor_using_group() {
  int prime_index = 0;
  order = 1;
  while (prime_index < static_cast<int>(prime_list_count)) {
    prime_index = next_primeform(initform, 0);
    exponentiation_stage();
    if (order > 1) {  // Exponentiate by known order
      qform_pow_u32(pow, initform, initform, order);
    }

    if (group->is_id(group, initform)) {
      // Very likely this group has odd order.
      break;
    }

    // Repeatedly square looking for an ambiguous form.
    group->set(group, search, initform);
    repeatedly_square(search);
    if (mpz_cmp_ui(d, 0) != 0) {
      // Found a divisor of N.
      return;
    }

    // Begin search stage.
    setup_search_stage();
    uint32_t new_order = search_stage();

    // Check if we found the order.
    if (new_order > 1) {
      // Exponentiate initform by the order and then repeatedly square.
      qform_pow_u32(pow, search, initform, new_order);
      repeatedly_square(search);
      if (mpz_cmp_ui(d, 0) != 0) {
	// Found a divisor of N.
	return;
      }
      // We can reuse the order.
      order *= new_order;
    }
  }
}

/**
 * Attempts to factor the integer N.
 * @param N The integer to be factored.
 * @param d A divisor of N.
 * @return 1 if a non-trivial factor of N is found, 0 otherwise.
 */
bool Spar::factor(mpz_t d, const mpz_t N) {
  mpz_set_ui(this->d, 0);
  mpz_set(this->N, N);
  
  double Nd = mpz_get_d(N);
  r = sqrt(log(Nd) / log(log(Nd)));
  t = ::pow(Nd, 1 / (2 * r));

  setup_exponentiation_stage();

  for (size_t multiplier_index = 0;
       multiplier_index < square_free_count;
       multiplier_index++) {
    unsigned int k = square_free[multiplier_index];
    cout << "k = " << k << endl;

    // Verify that k is not a divisor of N
    if ((k > 1) && (mpz_cmp_ui(N, k) != 0) && mpz_divisible_ui_p(N, k)) {
      mpz_set_ui(d, k);
      return true;
    }

    // Setup discriminant and members
    setup_discriminant(k);

    factor_using_group();
    if (mpz_cmp_ui(this->d, 0)) {
      mpz_set(d, this->d);
      return true;
    }

    // TODO: remove this line. for now we only want one execution.
    //    return false;
  }
  return false;
}
