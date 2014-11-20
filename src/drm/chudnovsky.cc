#include "drm/chudnovsky.h"

#include <glog/logging.h>
#include <gmp.h>
#include <sys/time.h>
#include <algorithm>
#include <cstdio>
#include <ctime>
#include <thread>

#include "base/base.h"

namespace pi {

namespace {

const double kDigsPerTerm = 14.181647462725477;

const int64 kConstA = 545140134;
const int64 kConstB = 13591409;
const int64 kConstC = 640320;

double GetTime() {
  timeval t;
  gettimeofday(&t, NULL);
  return t.tv_sec + t.tv_usec * 1e-6;
}
  
}  // namespace

void Chudnovsky::Compute(int64 digits, mpf_t pi) {
  double all_start = GetTime();

  int64 num_terms = digits / kDigsPerTerm + 5;
  LOG(INFO) << "Computing terms: " << num_terms;
  LOG(INFO) << "Target digits: " << digits;
  int num_threads = std::thread::hardware_concurrency() * 2;
  LOG(INFO) << "Number of threads: " << num_threads;

  mpz_t a, b, c;
  mpz_inits(a, b, c, NULL);

  double bs_start = GetTime();
  BinarySplit(0, num_terms, a, b, c, num_threads);
  double bs_end = GetTime();
  mpz_clear(c);
  LOG(INFO) << "Time of BS: " << (bs_end - bs_start) << " sec.";

  mpz_mul_ui(b, b, 12);

  int64 bits_a = mpz_sizeinbase(a, 2);
  int64 bits_b = mpz_sizeinbase(b, 2);
  LOG(INFO) << "Size of a: " << bits_a << " bits";
  LOG(INFO) << "Size of b: " << bits_b << " bits";
  int64 bits = std::max(bits_a, bits_b);
  mpf_set_default_prec(bits);

  mpf_t p, q;
  mpf_inits(p, q, NULL);

  mpf_set_z(p, a);
  mpf_set_z(q, b);

  mpz_clears(a, b, NULL);

  mpf_div(p, p, q);

  mpf_set_ui(q, kConstC * kConstC);
  mpf_mul_ui(q, q, kConstC);
  mpf_sqrt(q, q);
  mpf_mul(p, p, q);

  mpf_set_prec(pi, bits);
  mpf_set(pi, p);

  mpf_clears(p, q, NULL);

  double all_end = GetTime();
  LOG(INFO) << "Time of computing: " << (all_end - all_start) << " sec.";
}

void Chudnovsky::BinarySplit(int64 low, int64 up,
                             mpz_t a0, mpz_t b0, mpz_t c0, int num_cpu) {
  if (low + 1 == up) {
    SetValues(low, a0, b0, c0);
    return;
  }

  int64 mid = (low + up) / 2;
  int cpu0 = num_cpu / 2;
  int cpu1 = num_cpu - cpu0;  // cpu0 <= cpu1

  mpz_t a1, b1, c1;
  mpz_inits(a1, b1, c1, NULL);
  if (cpu0) {
    std::thread th0([&a0, &b0, &c0, &low, &mid, &cpu0] {
	BinarySplit(low, mid, a0, b0, c0, cpu0);
    });
    BinarySplit(mid, up, a1, b1, c1, cpu1);
    th0.join();
  } else {
    BinarySplit(low, mid, a0, b0, c0, cpu0);
    BinarySplit(mid, up, a1, b1, c1, cpu1);
  }

  mpz_mul(b0, b0, a1);
  mpz_mul(b1, b1, c0);
  mpz_add(b0, b0, b1);
  mpz_mul(a0, a0, a1);
  mpz_mul(c0, c0, c1);

  mpz_clears(a1, b1, c1, NULL);
}

void Chudnovsky::SetValues(int64 k, mpz_t a, mpz_t b, mpz_t c) {
  // a[k] = k^3 * C^3 / 24
  if (k == 0) {
    mpz_set_ui(a, 1);
  } else {
    int64 base = kConstC * k;
    mpz_set_ui(a, base / 24);
    mpz_mul_ui(a, a, base);
    mpz_mul_ui(a, a, base);
  }

  // b[k] = A * k + B;
  mpz_set_ui(b, kConstA * k);
  mpz_add_ui(b, b, kConstB);

  // c[k] = -(6k+1)(2k+1)(6k+5)
  mpz_set_si(c, -(6 * k + 1));
  mpz_mul_ui(c, c, 2 * k + 1);
  mpz_mul_ui(c, c, 6 * k + 5);
}

}  // namespace pi

