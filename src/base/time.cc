#include "base/time.h"

#include <sys/time.h>
#include "base/base.h"

namespace pi {
namespace base {

double GetTime() {
  timeval t;
  gettimeofday(&t, nullptr);
  return t.tv_sec + t.tv_usec * 1e-6;
}

}  // namespace base
}  // namespace pi
