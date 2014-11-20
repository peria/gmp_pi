#pragma once

#include "base/base.h"
#include <gmp.h>

namespace pi {

class Chudnovsky {
 public:
  static void Compute(int64 digits, mpf_t pi);

 protected:
  static void BinarySplit(int64 low, int64 up,
			  mpz_t a0, mpz_t b0, mpz_t c0, int cpu);

  // Set values of A[k], B[k], and C[k].
  static void SetValues(int64 k, mpz_t a, mpz_t b, mpz_t c);
};

}  // namespace pi
