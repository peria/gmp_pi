#include "drm/chudnovsky.h"

#include <glog/logging.h>
#include <gmp.h>
#include <algorithm>
#include "base/base.h"

namespace peria {

namespace {

const double kDigsPerTerm = 14.181647462725477;

const int64 kConstA = 13591409;
const int64 kConstB = 545140134;
const int64 kConstC = 640320;

}  // namespace

Chudnovsky::Chudnovsky() {}

Chudnovsky::~Chudnovsky() {}

void Chudnovsky::Init(int64 digits) {
  num_terms_ = digits / kDigsPerTerm + 5;
}

void Chudnovsky::Compute(mpf_t pi) {
  mpz_t a;
  mpz_t b;
  mpz_t c;
  mpz_inits(a, b, c, NULL);

  BinarySplit(0, num_terms_, a, b, c);

  mpz_clear(c);

  mpz_mul_ui(b, b, 12);

  int64 bits_a = mpz_sizeinbase(a, 2);
  int64 bits_b = mpz_sizeinbase(b, 2);
  LOG(INFO) << "Size of a: " << bits_a << " bits";
  LOG(INFO) << "Size of b: " << bits_b << " bits";
  mpf_set_default_prec(std::max(bits_a, bits_b));

  mpf_t p, q;
  mpf_inits(p, q, NULL);

  mpf_set_z(p, a);
  mpf_set_z(q, b);

  mpz_clears(a, b, NULL);

  mpf_div(p, p, q);

  mpf_set_ui(q, kConstC * kConstC * kConstC);
  mpf_sqrt(q, q);
  mpf_mul(p, p, q);

  mpf_set(pi, p);

  mpf_clears(p, q, NULL);
}

void Chudnovsky::BinarySplit(int64 low, int64 up, mpz_t a0, mpz_t b0, mpz_t c0) {
  if (low + 1 >= up) {
    SetValues(low, a0, b0, c0);
    return;
  }

  mpz_t a1;
  mpz_t b1;
  mpz_t c1;
  mpz_inits(a1, b1, c1, NULL);

  int64 mid = (low + up) / 2;
  BinarySplit(low, mid, a0, b0, c0);
  BinarySplit(mid, up, a1, b1, c1);

  mpz_mul(b0, b0, a1);
  mpz_mul(b1, b1, c0);
  mpz_add(b0, b0, b1);
  mpz_mul(a0, a0, a1);
  mpz_mul(c0, c0, c1);

  mpz_clears(a1, b1, c1, NULL);
}

void Chudnovsky::SetValues(int64 k, mpz_t a, mpz_t b, mpz_t c) {
  // a[k] = k^3 * C^3 / 24
  mpz_set_ui(a, kConstC / 24 * k);
  mpz_mul_ui(a, a, kConstC * k);
  mpz_mul_ui(a, a, kConstC * k);

  // b[k] = A * k + B;
  mpz_set_ui(b, kConstA * k);
  mpz_add_ui(b, b, kConstB);

  // c[k] = -(6k+1)(2k+1)(6k+5)
  mpz_set_ui(c, -(6 * k + 1));
  mpz_mul_ui(c, c, 2 * k + 1);
  mpz_mul_ui(c, c, 6 * k + 5);
}

}  // namespace peria

