#include "drm/ramanujan.h"

#include <glog/logging.h>
#include <gmp.h>
#include <sys/time.h>
#include <algorithm>
#include <cstdio>
#include <ctime>

#include "base/base.h"

namespace pi {

namespace {

const double kDigsPerTerm = 7.98254077839;  // log10(99^4)

const int64 kConstA = 26390;
const int64 kConstB = 1103;
const int64 kConstC = 99 * 99 * 99 * 99;

double GetTime() {
  timeval t;
  gettimeofday(&t, NULL);
  return t.tv_sec + t.tv_usec * 1e-6;
}
  
}  // namespace

void Ramanujan::Compute(int64 digits, mpf_t pi) {
  double all_start = GetTime();

  int64 num_terms = digits / kDigsPerTerm + 5;
  LOG(INFO) << "Computing terms: " << num_terms;
  LOG(INFO) << "Target digits: " << digits;

  ComputeCore(num_terms, pi);

  double all_end = GetTime();
  LOG(INFO) << "Time of computing: " << (all_end - all_start) << " sec.";
}

void Ramanujan::ComputeCore(int64 num_terms, mpf_t pi) {
  mpz_t a, b, c;
  mpz_inits(a, b, c, NULL);

  double bs_start = GetTime();
  BinarySplit(0, num_terms, a, b, c);
  double bs_end = GetTime();
  mpz_clear(c);
  LOG(INFO) << "Time of BS: " << (bs_end - bs_start) << " sec.";

  mpz_mul_ui(a, a, 99 * 99);
  mpz_mul_ui(b, b, 4);
  int64 bits_a = mpz_sizeinbase(a, 2);
  int64 bits_b = mpz_sizeinbase(b, 2);
  LOG(INFO) << "Size of a: " << bits_a << " bits";
  LOG(INFO) << "Size of b: " << bits_b << " bits";
  int64 bits = std::max(bits_a, bits_b);

  mpf_t q;
  mpf_init(q);

  mpf_set_prec(pi, bits);
  mpf_set_prec(q, bits);
  mpf_set_z(pi, a);
  mpf_set_z(q, b);

  mpz_clears(a, b, NULL);

  mpf_div(pi, pi, q);

  mpf_set_ui(q, 2);
  mpf_sqrt(q, q);
  mpf_mul(pi, pi, q);

  mpf_clear(q);
}

void Ramanujan::BinarySplit(int64 low, int64 up,
                            mpz_t a0, mpz_t b0, mpz_t c0) {
  if (low + 1 == up) {
    SetValues(low, a0, b0, c0);
    return;
  }

  int64 mid = (low + up) / 2;

  mpz_t a1, b1, c1;
  mpz_inits(a1, b1, c1, NULL);
  BinarySplit(low, mid, a0, b0, c0);
  BinarySplit(mid, up, a1, b1, c1);

  mpz_mul(b0, b0, a1);
  mpz_mul(b1, b1, c0);
  mpz_add(b0, b0, b1);
  mpz_mul(a0, a0, a1);
  mpz_mul(c0, c0, c1);

  mpz_clears(a1, b1, c1, NULL);
}

void Ramanujan::SetValues(int64 k, mpz_t a, mpz_t b, mpz_t c) {
  // a[k] = k^3 * C * 32 (for k > 0)
  if (k == 0) {
    mpz_set_ui(a, 1);
  } else {
    mpz_set_ui(a, 32 * kConstC * k);
    mpz_mul_ui(a, a, k * k);
  }

  // b[k] = A * k + B;
  mpz_set_ui(b, kConstA * k);
  mpz_add_ui(b, b, kConstB);

  // c[k] = (4k+1)(2k+1)(4k+3)
  mpz_set_ui(c, 4 * k + 1);
  mpz_mul_ui(c, c, 2 * k + 1);
  mpz_mul_ui(c, c, 4 * k + 3);
}

}  // namespace pi

