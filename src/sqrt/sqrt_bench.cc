#include <chrono>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include <gmp.h>
#include <gmpxx.h>

using int64 = std::int64_t;

namespace {

const int64 kMinBit = 1 << 20;
const int64 kMaxBit = 1 << 25;

const int64 d = 10005;
mpf_class answer;

}  // namespace

struct Timer {
  using Clock = std::chrono::system_clock;
  using Ms = std::chrono::milliseconds;
  Timer() : start(Clock::now()) {}
  double getTime() const {
    return std::chrono::duration_cast<Ms>(Clock::now() - start).count() * 1e-3;
  }
private:
  std::chrono::time_point<Clock> start;
};

struct Bench {
  virtual double compute(const int64 n) = 0;
  void getDiff(const mpf_class& r) {
    mpf_class rdiff = abs(answer - r);
    diff = rdiff.get_str(dexp, 10, 5);
  }
  std::string diff;
  mp_exp_t dexp;
};

struct GmpMpf : public Bench {
  double compute(const int64 n) override {
    mpf_class r;
    r.set_prec(n);
    Timer timer;
    mpf_sqrt_ui(r.get_mpf_t(), d);
    double t = timer.getTime();
    getDiff(r);
    return t;
  }
};

struct MpfInvNewton : public Bench {
  double compute(const int64 n) override {
    mpf_class r(1 / std::sqrt(d), 60);
    Timer timer;
    for (int64 k = 53; k < n; k *= 2) {
      r.set_prec(2 * k);
      r = (3 - d * r * r) * r / 2;
    }
    r.set_prec(n);
    r *= d;
    double t = timer.getTime();
    getDiff(r);
    return t;
  }
};

struct IntNewton : public Bench {
  double compute(const int64 n) override {
    mpf_class r;
    r.set_prec(n);
    Timer timer;
    mpz_class u = std::sqrt(d), v = 1;
    while (true) {
      long e;
      mpz_get_d_2exp(&e, u.get_mpz_t());
      if (e >= n)
        break;
      mpz_class u1 = u * u + d * v * v;
      mpz_class v1 = 2 * u * v;
      u = u1;
      v = v1;
    }
    mpf_class ru(u, n), rv(v, n);
    r = ru / rv;
    double t = timer.getTime();
    getDiff(r);
    return t;
  }
};

int main() {
  answer.set_prec(kMaxBit * 3 / 2);
  mpf_sqrt_ui(answer.get_mpf_t(), d);

  std::vector<std::unique_ptr<Bench>> benches;
  benches.emplace_back(new GmpMpf);
  benches.emplace_back(new MpfInvNewton);
  benches.emplace_back(new IntNewton);

  for (int64 n = kMinBit; n <= kMaxBit; n *= 2) {
    std::cout << std::setw(10) << n;
    for (auto& b : benches) {
      std::cout << " " << b->compute(n);
    }
    std::cout << "\n";
  }

  for (auto& b : benches) {
    std::cout << " " << b->diff << "e" << b->dexp;
  }
  std::cout << "\n";

  return 0;
}
