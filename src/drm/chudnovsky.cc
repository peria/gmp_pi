#include "drm/chudnovsky.h"

#include <glog/logging.h>
#include <gmp.h>
#include <sys/time.h>
#include <algorithm>
#include <cstdio>
#include <ctime>

#include "base/base.h"
#include "base/prime.h"
#include "base/time.h"

namespace pi {

namespace {

const double kDigsPerTerm = 14.181647462725477;

const int64 kConstA = 545140134;
const int64 kConstB = 13591409;
const int64 kConstC = 640320;
  
}  // namespace

void Chudnovsky::Compute(int64 digits, mpf_t pi) {
  double all_start = base::GetTime();

  int64 num_terms = digits / kDigsPerTerm + 5;
  LOG(INFO) << "Computing terms: " << num_terms;
  LOG(INFO) << "Target digits: " << digits;

  ComputeCore(num_terms, pi);

  double all_end = base::GetTime();
  LOG(INFO) << "Time of computing: " << (all_end - all_start) << " sec.";
}

namespace {
const int64 kTournamentLevel = 10;
const int64 kTournamentWidth = 1 << kTournamentLevel;
mpz_t g_factor[kTournamentLevel];
}  // namespace

void Chudnovsky::ComputeCore(int64 num_terms, mpf_t pi) {
  mpz_t a, b, c;
  mpz_inits(a, b, c, NULL);

  int64 num_tournament = (num_terms + kTournamentWidth - 1) / kTournamentWidth;
  LOG(INFO) << "Computing tournaments: " << num_tournament;
  InitTournament(kTournamentLevel);

  double bs_start = base::GetTime();
  if (num_tournament == 1) {
    LOG(INFO) << "Direct PT";
    PerfectTournament(0, kTournamentWidth, kTournamentLevel, a, b, c);
  } else {
    BinarySplit(0, num_tournament, a, b, c);
  }
  double bs_end = base::GetTime();

  mpz_set_ui(c, 1);
  for (int64 i = 0; i < kTournamentLevel; ++i) {
    mpz_mul(c, c, g_factor[i]);
    mpz_clear(g_factor[i]);
  }
  mpz_mul(a, a, c);
  
  mpz_mul_ui(c, a, kConstB);
  mpz_mul_ui(b, b, 5);
  mpz_sub(b, c, b);

  mpz_clear(c);
  
  LOG(INFO) << "Time of BS: " << (bs_end - bs_start) << " sec.";

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

  mpf_set_ui(q, (kConstC / 12) * (kConstC / 12) * kConstC);
  mpf_sqrt(q, q);
  mpf_mul(pi, pi, q);

  mpf_clear(q);
}

void Chudnovsky::BinarySplit(int64 low, int64 up,
                             mpz_t a0, mpz_t b0, mpz_t c0) {
  if (low + 1 == up) {
    PerfectTournament(low * kTournamentWidth, kTournamentWidth, kTournamentLevel,
                      a0, b0, c0);
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

void Chudnovsky::InitTournament(int64 level) {
  for (int64 i = 0; i < kTournamentLevel; ++i)
    mpz_init_set_ui(g_factor[i], 1);

  base::Prime primes(1 << level);
  for (int64 p = primes.GetNextPrime(); p > 0; p = primes.GetNextPrime()) {
    if (p == 2 || p == 3)
      continue;
    int64 width = 2;
    for (int64 i = 0; i < kTournamentLevel; ++i, width *= 2) {
      if ((width / p) % 2 == 1)
        mpz_mul_ui(g_factor[i], g_factor[i], p);
    }
  }
}

void Chudnovsky::PerfectTournament(int64 low, int64 width, int64 level,
                                   mpz_t a0, mpz_t b0, mpz_t c0) {
  DCHECK_EQ(1 << level, width);

  if (width == 1) {
    SetValues(low + 1, a0, b0, c0);
    return;
  }

  mpz_t a1, b1, c1;
  mpz_inits(a1, b1, c1, NULL);

  int64 half = width / 2;
  PerfectTournament(low,        half, level - 1, a0, b0, c0);
  PerfectTournament(low + half, half, level - 1, a1, b1, c1);

  mpz_mul(b0, b0, a1);
  mpz_mul(b1, b1, c0);
  mpz_add(b0, b0, b1);
  mpz_mul(a0, a0, a1);
  mpz_mul(c0, c0, c1);

  DCHECK_GE(kTournamentLevel, level);
  DCHECK_LT(0, level);
  mpz_divexact(a0, a0, g_factor[level - 1]);
  mpz_divexact(c0, c0, g_factor[level - 1]);

  mpz_clears(a1, b1, c1, NULL);
}

void Chudnovsky::SetValues(int64 k, mpz_t a, mpz_t b, mpz_t c) {
  // a[k] = k^3 * C^3 / 24
  int64 base = kConstC * k;
  mpz_set_ui(a, base / 24);
  mpz_mul_ui(a, a, base);
  mpz_mul_ui(a, a, base);

  // b[k] = A * k + B;
  mpz_set_ui(b, kConstA * k);
  mpz_add_ui(b, b, kConstB);

  // c[k] = -(6k+1)(2k+1)(6k+5)
  mpz_set_si(c, -(6 * k + 1));
  mpz_mul_ui(c, c, 2 * k + 1);
  mpz_mul_ui(c, c, 6 * k + 5);
}

}  // namespace pi

