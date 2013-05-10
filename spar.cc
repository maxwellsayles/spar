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
  mpz_init(this->t);

  // Initialize groups.
  s64_qform_group_init(&this->qform_group_s64);
  s64_qform_init(&this->qform_group_s64, &this->initform_s64);
  s64_qform_init(&this->qform_group_s64, &this->search_s64);
  s64_qform_init(&this->qform_group_s64, &this->current_s64);
  s64_qform_init(&this->qform_group_s64, &this->temp_form_s64);
  group_pow_init(&this->pow_s64, &this->qform_group_s64.desc.group);

  s128_qform_group_init(&this->qform_group_s128);
  s128_qform_init(&this->qform_group_s128, &this->initform_s128);
  s128_qform_init(&this->qform_group_s128, &this->search_s128);
  s128_qform_init(&this->qform_group_s128, &this->current_s128);
  s128_qform_init(&this->qform_group_s128, &this->temp_form_s128);
  group_pow_init(&this->pow_s128, &this->qform_group_s128.desc.group);

  mpz_qform_group_init(&this->qform_group_mpz);
  mpz_qform_init(&this->qform_group_mpz, &this->initform_mpz);
  mpz_qform_init(&this->qform_group_mpz, &this->search_mpz);
  mpz_qform_init(&this->qform_group_mpz, &this->current_mpz);
  mpz_qform_init(&this->qform_group_mpz, &this->temp_form_mpz);
  group_pow_init(&this->pow_mpz, &this->qform_group_mpz.desc.group);
    
  // Orders.
  this->order = 0;
    
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
  s64_qform_group_clear(&this->qform_group_s64);

  group_pow_clear(&this->pow_s128);
  s128_qform_clear(&this->qform_group_s128, &this->initform_s128);
  s128_qform_clear(&this->qform_group_s128, &this->search_s128);
  s128_qform_clear(&this->qform_group_s128, &this->current_s128);
  s128_qform_clear(&this->qform_group_s128, &this->temp_form_s128);
  s128_qform_group_clear(&this->qform_group_s128);

  group_pow_clear(&this->pow_mpz);
  mpz_qform_clear(&this->qform_group_mpz, &this->initform_mpz);
  mpz_qform_clear(&this->qform_group_mpz, &this->search_mpz);
  mpz_qform_clear(&this->qform_group_mpz, &this->current_mpz);
  mpz_qform_clear(&this->qform_group_mpz, &this->temp_form_mpz);
  mpz_qform_group_clear(&this->qform_group_mpz);

  mpz_clear(this->D);
  mpz_clear(this->t);
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
 * Takes this->initform and exponentiates it by this->order.
 * Then squares this until an ambigous form is found and attempts to factor N.
 * @param d The non-trivial factor if one is found.
 * @return 1 if a non-trivial factor was found.
 */
/*
static int test_ambiguous_form(sspar_t* this, mpz_t d, const mpz_t N) {
  int i;
  qform_group_t* qgroup = this->qform_group;
  group_t* group = &qgroup->group;
    
  // Exponentiate initform by the odd part of order.
  group_pow_naf_r2l_u32(this->pow,
			this->temp_form, this->initform, this->order);

  // If this is the identity, then the initial form has odd order,
  // and the identity will be the only ambiguous form.
  if (group->is_id(group, this->temp_form)) {
    cprintf(2, "The primeform has an odd order and so only a trivial factorization.\n");
    return 0;
  }
    
  // Since we change primeforms, we aren't guaranteed that we
  // know the odd part of the order anymore.  Therefore, we must
  // limit the number of squares performed, since the remaining
  // order might be odd.
  for (i = 0;
       i < this->logh && !qgroup->is_ambiguous(qgroup, this->temp_form);
       i ++) {
    group->square(group, this->temp_form, this->temp_form);
  }

  // Try to split the ambiguous form.
  if (qgroup->split_ambiguous(qgroup, d, N, this->temp_form)) {
    if (verbose_level >= 2) {
      cprintf(2, "Successful with ", N);
      group->print(group, this->temp_form);
      cprintf(2, "\n");
    }
    return 1;
  } else {
    if (verbose_level >= 2) {
      cprintf(2, "Failed with ", N);
      group->print(group, this->temp_form);
      cprintf(2, "\n");
    }
    return 0;
  }
}
*/
/**
 * Set all polymorphic variables and the group discriminant.
 * NOTE: Requires this->D and this->primorial_index to be set.
 */
void Spar::setup_discriminant(const mpz_t N, const unsigned int k) {
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
    this->pow                  = &this->pow_mpz;
    this->primorial_terms      = this->primorial_terms_s128;
    this->primorial_term_count = this->primorial_term_count_s128;
  }
  this->qform_group->set_discriminant(this->qform_group, this->D);
  this->group = &this->qform_group->group;

  this->logh = (logD + 1) / 2;
  if (this->logh < 1) this->logh = 1;
}

void Spar::setup_exponentiation_stage(const mpz_t N, const double t) {
  int w = prime_index_ge(t);
  primorial_terms_s64 = s64_primorial_terms[w-1];
  primorial_term_count_s64 = s64_primorial_term_counts[w-1];
  primorial_terms_s128 = s128_primorial_terms[w-1];
  primorial_term_count_s128 = s128_primorial_term_counts[w-1];

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

// The class group is determined, now we iterate on various prime ideals.
void Spar::factor_using_group() {
  next_primeform(initform, 0);
  cout << "initform = ";
  group->print(group, initform);
  cout << endl;
  
  exponentiation_stage();
  cout << "After exponentiation stage ";
  group->print(group, initform);
  cout << endl;

  group->set(group, search, initform);
  repeatedly_square(search);

  cout << "Search = ";
  group->print(group, search);
  cout << endl;
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
  double r = sqrt(log(Nd) / log(log(Nd)));
  double t = ::pow(Nd, 1 / (2 * r));
  cout << setprecision(5) << fixed << "t=" << t << endl;

  setup_exponentiation_stage(N, t);

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
    setup_discriminant(N, k);

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
