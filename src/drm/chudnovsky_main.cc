#include "drm/chudnovsky.h"

#include <gflags/gflags.h>
#include <glog/logging.h>
#include <gmp.h>
#include <cstdio>

DEFINE_int64(digits, 100, "The number of digits to compute.");
DEFINE_bool(verbose, false, "If set true, outputs all digits.");

int main(int argc, char* argv[]) {
  google::ParseCommandLineFlags(&argc, &argv, true);
  google::InitGoogleLogging(argv[0]);

  pi::Chudnovsky chudnovsky;
  chudnovsky.Init(FLAGS_digits);

  mpf_t pi;
  mpf_init(pi);
  chudnovsky.Compute(pi);
  if (FLAGS_verbose || FLAGS_digits < 100)
    mpf_out_str(stdout, 10, FLAGS_digits, pi);
  else
    mpf_out_str(stdout, 10, 100, pi);
  mpf_clear(pi);
  
  return 0;
}
