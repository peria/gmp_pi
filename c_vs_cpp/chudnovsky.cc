/* Pi computation using Chudnovsky's algortithm.

 * Copyright 2002, 2005 Hanhong Xue (macroxue at yahoo dot com)

 * Slightly modified 2005 by Torbjorn Granlund to allow more than 2G
   digits to be computed.

 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 * 1. Redistributions of source code must retain the above copyright notice,
 * this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 * this list of conditions and the following disclaimer in the documentation
 * and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE AUTHORS ``AS IS'' AND ANY EXPRESS OR
 * IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
 * MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO
 * EVENT SHALL THE AUTHORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
 * OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
 * WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
 * OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
 * ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <vector>

#include "gmp.h"
#include "gmpxx.h"

constexpr double BITS_PER_DIGIT = 3.32192809488736234787;

enum OutputFlag {
  kNone = 0,
  kResult = 1,
  kDebug = 2,
};

class Chudnovsky {
  struct Sieve {
    int fac = 0;
    int pow = 0;
    int nxt = 0;
  };
  struct Factorized {
    static constexpr int64_t INIT_FACS = 32;

    Factorized() : fac(INIT_FACS), pow(INIT_FACS) {}

    // g = GCD(a, b); a /= g; b /= g;
    static void removeGcd(mpz_class& a,
                          Factorized& fa,
                          mpz_class& b,
                          Factorized& fb);

    // remove factors of power 0
    void compact();

    void toMpzClass(mpz_class& r) const;
    void resize(int64_t s);

    int64_t max_facs = 0;
    int64_t num_facs = 0;
    std::vector<int64_t> fac;
    std::vector<int64_t> pow;

   private:
    void toMpzClassBs(mpz_class& r, int64_t a, int64_t b) const;
  };

 public:
  Chudnovsky(int64_t digits, int output_flags);
  void compute();

  void outputResult();

  bool doesOutputResult() { return output_flags_ & OutputFlag::kResult; }
  bool doesOutputDebug() { return output_flags_ & OutputFlag::kDebug; }

 private:
  void buildSieve(int64_t n);
  void bs(uint64_t a, uint64_t b, bool gflag, int64_t level);
  void bsSet(uint64_t b);

  // r = base^pow
  void setFactor(Factorized& r, int64_t base, int64_t pow) const;
  // r *= base^pow
  void mulFactor(Factorized& r, int64_t base, int64_t pow) const;
  // r *= f
  void mulFactor(Factorized& r, Factorized& f) const;
  // r = f * g
  void mulFactor(Factorized& r, const Factorized& f, const Factorized& g) const;

  static constexpr int64_t A = 13591409;
  static constexpr int64_t B = 545140134;
  static constexpr int64_t C = 640320;
  static constexpr int64_t D = 12;

  static constexpr double DIGITS_PER_ITER = 14.1816474627254776555;

  const int64_t digits_;
  const int output_flags_;

  const int64_t terms_;
  const double progress_percent_;
  double progress_ = 0;

  int64_t depth_ = 1;
  std::vector<Sieve> sieve_;

  // Pseudo stacks
  std::vector<mpz_class> pstack_;
  std::vector<mpz_class> qstack_;
  std::vector<mpz_class> gstack_;
  std::vector<Factorized> fpstack_;
  std::vector<Factorized> fgstack_;
  int64_t stack_top_ = 0;

  mpf_class qi_;
};

using Clock = std::chrono::system_clock;
template <typename T>
double timeDiff(const T& from, const T& to) {
  using MS = std::chrono::milliseconds;
  return std::chrono::duration_cast<MS>(to - from).count() * 1e-3;
}

void Chudnovsky::Factorized::removeGcd(mpz_class& p,
                                       Factorized& fp,
                                       mpz_class& g,
                                       Factorized& fg) {
  Factorized fmul;
  fmul.resize(std::min(fp.num_facs, fg.num_facs));

  int64_t k = 0;
  for (int64_t i = 0, j = 0; i < fp.num_facs && j < fg.num_facs;) {
    if (fp.fac[i] == fg.fac[j]) {
      int64_t c = std::min(fp.pow[i], fg.pow[j]);
      fp.pow[i] -= c;
      fg.pow[j] -= c;
      fmul.fac[k] = fp.fac[i];
      fmul.pow[k] = c;
      ++i;
      ++j;
      ++k;
    } else if (fp.fac[i] < fg.fac[j]) {
      ++i;
    } else {
      ++j;
    }
  }
  fmul.num_facs = k;

  if (fmul.num_facs) {
    mpz_class gcd;
    fmul.toMpzClass(gcd);
    mpz_divexact(p.get_mpz_t(), p.get_mpz_t(), gcd.get_mpz_t());
    mpz_divexact(g.get_mpz_t(), g.get_mpz_t(), gcd.get_mpz_t());
    fp.compact();
    fg.compact();
  }
}

void Chudnovsky::Factorized::compact() {
  int64_t j = 0;
  for (int64_t i = 0; i < num_facs; ++i) {
    if (pow[i] == 0)
      continue;
    if (j < i) {
      fac[j] = fac[i];
      pow[j] = pow[i];
    }
    ++j;
  }
  num_facs = j;
}

void Chudnovsky::Factorized::toMpzClass(mpz_class& r) const {
  toMpzClassBs(r, 0, num_facs);
}

void Chudnovsky::Factorized::toMpzClassBs(mpz_class& r,
                                          int64_t a,
                                          int64_t b) const {
  if (b - a <= 32) {
    r = 1;
    for (int64_t i = a; i < b; ++i)
      for (int64_t j = 0; j < pow[i]; ++j)
        r *= fac[i];
  } else {
    mpz_class r2;
    const int64_t m = (a + b) / 2;
    toMpzClassBs(r2, a, m);
    toMpzClassBs(r, m, b);
    r *= r2;
  }
}

void Chudnovsky::Factorized::resize(int64_t s) {
  if (max_facs < s) {
    if (s < INIT_FACS)
      s = INIT_FACS;

    max_facs = s;
    num_facs = 0;
    fac.resize(s);
    pow.resize(s);
  }
}

Chudnovsky::Chudnovsky(int64_t digits, int output_flags)
    : digits_(digits),
      output_flags_(output_flags),
      terms_(digits / DIGITS_PER_ITER),
      progress_percent_(terms_ * 0.01) {
  for (depth_ = 1; (1LL << depth_) < terms_; ++depth_) {
  }
  ++depth_;

  printf("#terms=%ld, depth=%ld\n", terms_, depth_);
  pstack_.resize(depth_);
  qstack_.resize(depth_);
  gstack_.resize(depth_);
  fpstack_.resize(depth_);
  fgstack_.resize(depth_);
}

void Chudnovsky::compute() {
  auto begin = Clock::now();
  printf("sieve   ");
  fflush(stdout);
  int64_t sieve_size = std::max<int64_t>(3 * 5 * 23 * 29 + 1, terms_ * 6);
  buildSieve(sieve_size);
  auto mid0 = Clock::now();
  printf("time = %6.3f\n", timeDiff(begin, mid0));
  bs(0, terms_, false, 0);
  auto mid1 = Clock::now();
  printf("\nbs      time = %6.3f\n", timeDiff(mid0, mid1));

  pstack_.resize(1);
  qstack_.resize(1);
  gstack_.clear();
  fpstack_.clear();
  fgstack_.clear();

  //       p*(C/D)*sqrt(C)
  // pi = -----------------
  //           (q+A*p)

  // prepare to convert integers to floats
  int64_t default_prec = static_cast<int64_t>(digits_ * BITS_PER_DIGIT + 16);
  mpf_set_default_prec(default_prec);

  mpz_class& p1 = pstack_.front();
  mpz_class& q1 = qstack_.front();
  int64_t psize = mpz_sizeinbase(p1.get_mpz_t(), 10);
  int64_t qsize = mpz_sizeinbase(q1.get_mpz_t(), 10);

  q1 += p1 * A;
  p1 *= C / D;

  mpf_class pi(p1);
  pstack_.clear();

  mpf_class& qi = qi_;
  qi.set_prec(default_prec);
  qi = q1;
  qstack_.clear();

  auto mid2 = Clock::now();

  // final step
  printf("div     ");
  fflush(stdout);
  qi = pi / qi;
  auto mid3 = Clock::now();
  printf("time = %6.3f\n", timeDiff(mid2, mid3));

  printf("sqrt    ");
  fflush(stdout);
  mpf_sqrt_ui(pi.get_mpf_t(), C);
  auto mid4 = Clock::now();
  printf("time = %6.3f\n", timeDiff(mid3, mid4));

  printf("mul     ");
  fflush(stdout);
  qi *= pi;
  auto end = Clock::now();
  printf("time = %6.3f\n", timeDiff(mid4, end));

  printf("total   time = %6.3f\n", timeDiff(begin, end));
  fflush(stdout);

  printf(
      "   P size=%ld digits (%f)\n"
      "   Q size=%ld digits (%f)\n",
      psize, 1.0 * psize / digits_, qsize, 1.0 * qsize / digits_);
}

void Chudnovsky::outputResult() {
  printf("pi(0,%ld)=\n", terms_);
  mpf_out_str(stdout, 10, digits_ + 2, qi_.get_mpf_t());
  printf("\n");
}

void Chudnovsky::buildSieve(int64_t n) {
  const int64_t m = static_cast<int64_t>(std::sqrt(n));
  sieve_.resize(n);

  sieve_[1 / 2].fac = 1;
  sieve_[1 / 2].pow = 1;

  for (int64_t i = 3; i <= n; i += 2) {
    Sieve& si = sieve_[i / 2];
    if (si.fac)
      continue;
    si.fac = i;
    si.pow = 1;
    if (i > m)
      continue;
    for (int64_t j = i * i, k = i / 2; j <= n; j += 2 * i, ++k) {
      Sieve& sj = sieve_[j / 2];
      if (sj.fac)
        continue;
      sj.fac = i;
      Sieve& sk = sieve_[k];
      if (sk.fac == i) {
        sj.pow = sk.pow + 1;
        sj.nxt = sk.nxt;
      } else {
        sj.pow = 1;
        sj.nxt = k;
      }
    }
  }
}

void Chudnovsky::bs(uint64_t a, uint64_t b, bool gflag, int64_t level) {
  if (b - a == 1) {
    bsSet(b);
  } else {
    // p(a,b) = p(a,m) * p(m,b)
    // g(a,b) = g(a,m) * g(m,b)
    // q(a,b) = q(a,m) * p(m,b) + q(m,b) * g(a,m)
    int64_t mid = a + (b - a) * 0.5224;  // tuning parameter
    bs(a, mid, true, level + 1);

    ++stack_top_;
    bs(mid, b, gflag, level + 1);
    --stack_top_;

    mpz_class& p1 = pstack_[stack_top_];
    mpz_class& q1 = qstack_[stack_top_];
    mpz_class& g1 = gstack_[stack_top_];
    Factorized& fp1 = fpstack_[stack_top_];
    Factorized& fg1 = fgstack_[stack_top_];
    mpz_class& p2 = pstack_[stack_top_ + 1];
    mpz_class& q2 = qstack_[stack_top_ + 1];
    mpz_class& g2 = gstack_[stack_top_ + 1];
    Factorized& fp2 = fpstack_[stack_top_ + 1];
    Factorized& fg2 = fgstack_[stack_top_ + 1];

    if (level == 0)
      puts("");
    if (level >= 4)  // tuning parameter
      Factorized::removeGcd(p2, fp2, g1, fg1);

    mpz_mul(p1.get_mpz_t(), p1.get_mpz_t(), p2.get_mpz_t());
    mpz_mul(q1.get_mpz_t(), q1.get_mpz_t(), p2.get_mpz_t());
    mpz_mul(q2.get_mpz_t(), q2.get_mpz_t(), g1.get_mpz_t());
    mpz_add(q1.get_mpz_t(), q1.get_mpz_t(), q2.get_mpz_t());
    mulFactor(fp1, fp2);

    if (gflag) {
      mpz_mul(g1.get_mpz_t(), g1.get_mpz_t(), g2.get_mpz_t());
      mulFactor(fg1, fg2);
    }
  }
}

void Chudnovsky::bsSet(uint64_t b) {
  // g(b-1,b) = (6b-5)(2b-1)(6b-1)
  // p(b-1,b) = b^3 * C^3 / 24
  // q(b-1,b) = (-1)^b*g(b-1,b)*(A+Bb).
  mpz_class& p1 = pstack_[stack_top_];
  mpz_class& q1 = qstack_[stack_top_];
  mpz_class& g1 = gstack_[stack_top_];
  Factorized& fp1 = fpstack_[stack_top_];
  Factorized& fg1 = fgstack_[stack_top_];

  p1 = b;
  p1 *= b;
  p1 *= b;
  p1 *= (C / 24) * (C / 24);
  p1 *= C * 24;

  g1 = 2 * b - 1;
  g1 *= 6 * b - 1;
  g1 *= 6 * b - 5;

  q1 = b;
  q1 *= B;
  q1 += A;
  q1 *= g1;
  if (b % 2)
    q1 = -q1;

  int64_t i = b;
  while ((i & 1) == 0)
    i >>= 1;
  setFactor(fp1, i, 3);
  mulFactor(fp1, 3 * 5 * 23 * 29, 3);
  fp1.pow[0]--;
  setFactor(fg1, 2 * b - 1, 1);  // 2b-1
  mulFactor(fg1, 6 * b - 1, 1);  // 6b-1
  mulFactor(fg1, 6 * b - 5, 1);  // 6b-5

  if (b > progress_) {
    printf(".");
    fflush(stdout);
    progress_ += progress_percent_ * 2;
  }
}

void Chudnovsky::setFactor(Factorized& f, int64_t base, int64_t pow) const {
  int64_t i;
  for (i = 0, base /= 2; base > 0; ++i, base = sieve_[base].nxt) {
    f.fac[i] = sieve_[base].fac;
    f.pow[i] = sieve_[base].pow * pow;
  }
  f.num_facs = i;
}

void Chudnovsky::mulFactor(Factorized& f, int64_t base, int64_t pow) const {
  Factorized tmp;
  setFactor(tmp, base, pow);
  mulFactor(f, tmp);
}

void Chudnovsky::mulFactor(Factorized& f, Factorized& g) const {
  Factorized tmp;
  tmp.resize(f.num_facs + g.num_facs);
  mulFactor(tmp, f, g);
  std::swap(f, tmp);
}

void Chudnovsky::mulFactor(Factorized& r,
                           const Factorized& f,
                           const Factorized& g) const {
  int64_t i = 0;
  int64_t j = 0;
  int64_t k = 0;

  for (i = j = k = 0; i < f.num_facs && j < g.num_facs; ++k) {
    if (f.fac[i] == g.fac[j]) {
      r.fac[k] = f.fac[i];
      r.pow[k] = f.pow[i] + g.pow[j];
      ++i;
      ++j;
    } else if (f.fac[i] < g.fac[j]) {
      r.fac[k] = f.fac[i];
      r.pow[k] = f.pow[i];
      ++i;
    } else {
      r.fac[k] = g.fac[j];
      r.pow[k] = g.pow[j];
      ++j;
    }
  }
  for (; i < f.num_facs; ++i, ++k) {
    r.fac[k] = f.fac[i];
    r.pow[k] = f.pow[i];
  }
  for (; j < g.num_facs; ++j, ++k) {
    r.fac[k] = g.fac[j];
    r.pow[k] = g.pow[j];
  }
  r.num_facs = k;
}

int main(int argc, char* argv[]) {
  const int64_t digits = (argc > 1) ? std::strtoll(argv[1], nullptr, 10) : 100;
  const int output_flags = (argc > 2) ? std::atoi(argv[2]) : 0;

  Chudnovsky pi(digits, output_flags);

  pi.compute();

  if (pi.doesOutputResult()) {
    pi.outputResult();
  }

  return 0;
}
