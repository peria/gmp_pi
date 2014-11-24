#include "base/prime.h"

#include <glog/logging.h>

#include <algorithm>
#include <cmath>
#include <vector>

#include "base/base.h"

namespace pi {
namespace base {

namespace {

const size_t kOffset[] = {1, 7, 11, 13, 17, 19, 23, 29};
const uint8 kBitMask[][8] = {
    {0xFE, 0xFD, 0xFB, 0xF7, 0xEF, 0xDF, 0xBF, 0x7F},
    {0xFD, 0xDF, 0xEF, 0xFE, 0x7F, 0xF7, 0xFB, 0xBF},
    {0xFB, 0xEF, 0xFE, 0xBF, 0xFD, 0x7F, 0xF7, 0xDF},
    {0xF7, 0xFE, 0xBF, 0xDF, 0xFB, 0xFD, 0x7F, 0xEF},
    {0xEF, 0x7F, 0xFD, 0xFB, 0xDF, 0xBF, 0xFE, 0xF7},
    {0xDF, 0xF7, 0x7F, 0xFD, 0xBF, 0xFE, 0xEF, 0xFB},
    {0xBF, 0xFB, 0xF7, 0x7F, 0xFE, 0xEF, 0xDF, 0xFD},
    {0x7F, 0xBF, 0xDF, 0xEF, 0xF7, 0xFB, 0xFD, 0xFE}};
const int64  kNextIndex[][8] = {
    {0, 0, 0, 0, 0, 0, 0, 1},
    {1, 1, 1, 0, 1, 1, 1, 1},
    {2, 2, 0, 2, 0, 2, 2, 1},
    {3, 1, 1, 2, 1, 1, 3, 1},
    {3, 3, 1, 2, 1, 3, 3, 1},
    {4, 2, 2, 2, 2, 2, 4, 1},
    {5, 3, 1, 4, 1, 3, 5, 1},
    {6, 4, 2, 4, 2, 4, 6, 1}};
const int64 kCoef[] = {6, 4, 2, 4, 2, 4, 6, 2};

}  // namespace

Prime::Prime(const int n) : n_(n), index_(-1), bit_(0) {
  bits_.clear();
  bits_.resize((n + 29) / 30, 0xff);
  bits_[0] = 0xfe;
  Sieve();
}

int Prime::GetNextPrime() {
  if (index_ < 0) {
    switch (bit_) {
    case 0:
      ++bit_;
      return 2;
    case 1:
      ++bit_;
      return 3;
    case 2:
      index_ = 0;
      bit_ = 0;
      return 5;
    }
  }

  // Search next prime
  do {
    if (++bit_ > 7) {
      bit_ = 0;
      if (++index_ >= static_cast<int>(bits_.size()))
        break;
    }
  } while ((bits_[index_] & (1 << bit_)) == 0);

  if (index_ >= static_cast<int>(bits_.size()))
    return -1;
  int prime = index_ * 30 + kOffset[bit_];

  return (prime <= n_) ? prime : -1;
}

namespace {
inline int Bit2Index(uint8 bit) {
  switch (bit) {
  case (1 << 0): return 0;
  case (1 << 1): return 1;
  case (1 << 2): return 2;
  case (1 << 3): return 3;
  case (1 << 4): return 4;
  case (1 << 5): return 5;
  case (1 << 6): return 6;
  case (1 << 7): return 7;
  }
  return -1;
}
}  // namespace

void Prime::Sieve() {
  const size_t imax = (std::sqrt(n_) + 30) / 30;
  for (size_t i = 0; i < imax; ++i) {
    for (uint8 bits = bits_[i]; bits; bits &= bits - 1) {
      uint8 bit = bits & -bits;
      int bid = Bit2Index(bit);
      const int64 offset = kOffset[bid];
      const int64 p = 30 * i + offset;
      size_t k = (p + offset) * i + offset * offset / 30;
      int l = bid;
      while (k < bits_.size()) {
        bits_[k] &= kBitMask[bid][l];
        k += kCoef[l] * i + kNextIndex[bid][l];
        l = (l + 1) % 8;
      }
    }
  }
}

int Prime::GetPrimes(int n, std::vector<int>* primes) {
  if (n > n_)
    return -1;

  primes->clear();
  if (2 <= n)
    primes->push_back(2);
  if (3 <= n)
    primes->push_back(3);
  if (5 <= n)
    primes->push_back(5);
  for (int i = 0; i < n / 30; ++i) {
    for (uint8 bits = bits_[i]; bits; bits &= bits - 1) {
      int bid = Bit2Index(bits & -bits);
      int p = 30 * i + kOffset[bid];
      primes->push_back(p);
    }
  }

  size_t i = n / 30;
  for (uint8 bits = bits_[i]; bits; bits &= bits - 1) {
    int bid = Bit2Index(bits & -bits);
    int p = 30 * i + kOffset[bid];
    if (p > n)
      break;
    primes->push_back(p);
  }

  return primes->size();
}

namespace {
int BitCnt(uint8 bits) {
  bits = (bits & 0x55) + ((bits >> 1) & 0x55);
  bits = (bits & 0x33) + ((bits >> 2) & 0x33);
  return (bits + (bits >> 4)) & 0xf;
}
}  // namespace

int Prime::CountPrimes(int n) {
  if (n > n_)
    return -1;

  int ret = 0;
  if (2 <= n)
    ++ret;
  if (3 <= n)
    ++ret;
  if (5 <= n)
    ++ret;
  for (int i = 0; i < n / 30; ++i)
    ret += BitCnt(bits_[i]);

  size_t i = n / 30;
  for (uint8 bits = bits_[i]; bits; bits &= bits - 1) {
    int bid = Bit2Index(bits & -bits);
    int p = 30 * i + kOffset[bid];
    if (p > n)
      break;
    ++ret;
  }

  return ret;
}

bool Prime::IsPrime(int64 n) {
  if (n < 2)
    return false;
  if (n >= n_)
    return IsPrimeSt(n);
  if (n == 2 || n == 3 || n == 5)
    return true;
  if (n % 2 == 0 || n % 3 == 0 || n % 5 == 0)
    return false;

  uint8 bit = bits_[n / 30];
  switch (n % 30) {
  case  1: return (bit & 0x01) != 0;
  case  7: return (bit & 0x02) != 0;
  case 11: return (bit & 0x04) != 0;
  case 13: return (bit & 0x08) != 0;
  case 17: return (bit & 0x10) != 0;
  case 19: return (bit & 0x20) != 0;
  case 23: return (bit & 0x40) != 0;
  case 29: return (bit & 0x80) != 0;
  }

  return false;
}

void Prime::ResetIndex() {
  index_ = -1;
  bit_ = 0;
}

// static method
bool Prime::IsPrimeSt(int64 n) {
  if (n < 2)
    return false;
  if (n == 2)
    return true;
  if (n % 2 == 0)
    return false;
  for (int64 p = 3; p * p <= n; p += 2)
    if (n % p == 0)
      return false;
  return true;
}

}  // namespace base
}  // namespace pi
