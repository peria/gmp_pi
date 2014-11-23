#include <gflags/gflags.h>
#include <glog/logging.h>
#include <gmp.h>
#include <cstdio>

#include "drm/chudnovsky.h"
#include "drm/ramanujan.h"

DEFINE_int64(digits, 100, "The number of digits to compute.");
DEFINE_bool(verbose, false, "If set true, outputs all digits.");
DEFINE_int32(formula, 0,
             "Choose a formula to compute.\n"
             "\t0: Chudnovsky (default)\n"
             "\t1: Ramanujan");

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
  if (FLAGS_verbose || FLAGS_digits < 100)
    mpf_out_str(stdout, 10, FLAGS_digits, pi);
  else
    mpf_out_str(stdout, 10, 100, pi);
  mpf_clear(pi);
  
  return 0;
}
