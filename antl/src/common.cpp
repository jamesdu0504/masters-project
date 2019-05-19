/**
 * @file common.cpp
 * @author Michael Jacobson
 * @brief implementation of non-templated functions defined in common.hpp
 */

#include <ANTL/common.hpp>

namespace ANTL {

  void DivRem(long & q, long & r, long a, long b)
  {
    q = (long)std::floor((double)a/(double)b);
    r = a - q*b;
  }

} // ANTL
