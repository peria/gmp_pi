#include <gflags/gflags.h>
#include <glog/logging.h>
#include <gmp.h>
#include <iostream>
#include <cstdio>

#include "drm/chudnovsky.h"
#include "drm/ramanujan.h"

DEFINE_int64(digits, 100, "The number of digits to compute.");
DEFINE_bool(verbose, false, "If set true, outputs all digits.");
DEFINE_int32(formula, 0, "Choose a formula to compute.\n"
             "\t0: Chudnovsky (default)\n"
             "\t1: Ramanujan");
DEFINE_int32(check, 1, "Check if output pi has a right value.\n"
             "\t0: No check,\n"
             "\t1: Compare few digits,\n"
             "\t2: Compare full digits.");

namespace {
enum Formula {
  kChudnovsky = 0,
  kRamanujan = 1
};
}

int main(int argc, char* argv[]) {
  google::ParseCommandLineFlags(&argc, &argv, true);
  google::InitGoogleLogging(argv[0]);

  mpf_t pi;
  mpf_init(pi);
  if (FLAGS_formula == kRamanujan) {
    pi::Ramanujan::Compute(FLAGS_digits, pi);
  } else {
    pi::Chudnovsky::Compute(FLAGS_digits, pi);
  }

  if (FLAGS_check) {
    mpf_t answer;
    mpf_init2(answer, 300);
    mpf_set_str(answer, "3.141592653589793238462643383279502884197169399", 10);
    mpf_sub(answer, answer, pi);
    int64 e;
    double d = mpf_get_d_2exp(&e, answer);
    if (e > -280) {
      std::cout << "** Wrong value **\n";
    }
    mpf_clear(answer);
  }

  if (FLAGS_verbose || FLAGS_digits < 100)
    mpf_out_str(stdout, 10, FLAGS_digits, pi);
  else
    mpf_out_str(stdout, 10, 100, pi);
  mpf_clear(pi);
  
  return 0;
}
