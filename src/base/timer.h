#pragma once

#include <chrono>

namespace pi {
namespace base {

class Timer {
  using Clock = std::chrono::steady_clock;
  using Micro = std::chrono::microseconds;
 public:
  Timer() : start_(Clock::now()) {}
  void Stop() { end_ = Clock::now(); }

  double TimeInSec() {
    return std::chrono::duration_cast<Micro>(end_ - start_).count() * 1e-6;
  }

 private:
  Clock::time_point start_;
  Clock::time_point end_;
};

}  // namespace base
}  // namespace pi
