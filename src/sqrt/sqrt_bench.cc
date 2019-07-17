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
    r.set_prec(n + n / 10);
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
    r.set_prec(n + n / 10);
    r *= d;
    double t = timer.getTime();
    getDiff(r);
    return t;
  }
};

struct IntNewton : public Bench {
  double compute(const int64 n) override {
    mpf_class r;
    r.set_prec(n + n / 10);
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
    mpf_class ru(u, n + n / 10), rv(v, n + n / 10);
    r = ru / rv;
    double t = timer.getTime();
    getDiff(r);
    return t;
  }
};

namespace {

using Matrix = std::array<std::array<mpz_class, 2>, 2>;

Matrix operator*(const Matrix& a, const Matrix& b) {
  Matrix c;
  c[0][0] = a[0][0] * b[0][0] + a[0][1] * b[1][0];
  c[1][0] = a[1][0] * b[0][0] + a[1][1] * b[1][0];
  c[0][1] = a[0][0] * b[0][1] + a[0][1] * b[1][1];
  c[1][1] = a[1][0] * b[0][1] + a[1][1] * b[1][1];
  return c;
}

Matrix Power(Matrix m, int64 e) {
  Matrix r;
  r[0][0] = 1;
  r[0][1] = 0;
  r[1][0] = 0;
  r[1][1] = 1;
  for (int64 b = 1 << 30; b; b >>= 1) {
    r = r * r;
    if (e & b)
      r = r * m;
  }
  return r;
}

}  // namespace

struct Eigen : public Bench {
  double compute(const int64 n) override {
    static const int64 a = std::ceil(std::sqrt(d));
    const double alpha = std::log2(a + std::sqrt(d));
    const double beta = std::log2(a - std::sqrt(d));
    int64 k = (n + 1 + std::log2(d) / 2) / (alpha - beta);
    Matrix m;
    m[0][0] = a;
    m[0][1] = d;
    m[1][0] = 1;
    m[1][1] = a;
    Timer timer;
    Matrix rs = Power(m, k);
    mpz_class zr(a * rs[0][0] + rs[0][1]);
    mpz_class zs(a * rs[1][0] + rs[1][1]);
    mpf_class r(zr, n + n / 10);
    mpf_class s(zs, n + n / 10);
    r /= s;
    double t = timer.getTime();
    getDiff(r);
    return t;
  }
};

struct CFrac : public Bench {
  CFrac() {
    {
      const int64 sd = std::sqrt(d);
      int64 s0 = 0;
      int64 t0 = 1;
      int64 q0 = sd;
      q.push_back(q0 * 2);
      while (true) {
        int64 s1 = q0 * t0 - s0;
        int64 t1 = (d - s1 * s1) / t0;
        int64 q1 = (s1 + sd) / t1;
        s0 = s1;
        t0 = t1;
        q0 = q1;
        if (q0 == q[0])
          break;
        q.push_back(q0);
      }
    }
    m = makeMatrix(0, q.size());
    {
      double p0 = m[0][0].get_d();
      double q0 = m[0][1].get_d();
      double p1 = m[1][0].get_d();
      double q1 = m[1][1].get_d();
      alpha = (p0 + q1) / 2 + q0 * std::sqrt(d);
      beta = (p0 + q1) / 2 - q0 * std::sqrt(d);
    }
  }

  Matrix makeMatrix(int64 low, int64 high) {
    if (low + 1 == high) {
      Matrix d;
      d[0][0] = q[low];
      d[0][1] = 1;
      d[1][0] = 1;
      d[1][1] = 0;
      return d;
    }
    int64 mid = (low + high) / 2;
    return makeMatrix(mid, high) * makeMatrix(low, mid);
  }

  double compute(const int64 n) override {
    int64 k = (n + 1 - std::log2(d) / 2) / (alpha - beta) * 310 + 1;
    Timer timer;
    Matrix am = Power(m, k);
    mpf_class r(am[0][0], n + n / 10);
    mpf_class s(am[0][1], n + n / 10);
    r /= s;
    r -= static_cast<int>(std::sqrt(d));
    double t = timer.getTime();
    getDiff(r);
    return t;
  }

  std::vector<int64> q;
  Matrix m;
  double alpha;
  double beta;
};

int main() {
  answer.set_prec(kMaxBit * 3 / 2);
  mpf_sqrt_ui(answer.get_mpf_t(), d);

  std::vector<std::unique_ptr<Bench>> benches;
  benches.emplace_back(new GmpMpf);
  benches.emplace_back(new MpfInvNewton);
  benches.emplace_back(new IntNewton);
  benches.emplace_back(new Eigen);
  benches.emplace_back(new CFrac);

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
