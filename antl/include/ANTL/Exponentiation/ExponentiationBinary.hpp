/**
 * @file ExponentiationBinary.hpp
 * @author Michael Jacobson
 * @brief class for left-to-right binary exponentiation
 */

#ifndef EXPONENTIATION_BINARY_H
#define EXPONENTIATION_BINARY_H

#include "Exponentiation.hpp"

namespace ANTL
{

  /**
   * @brief class for left-to-right binary exponentiation
   * @remarks This concrete class defines a method for standard left-to-right binary
   * exponentiation.  The base type is templated, and must have the following
   * functions defined:
   *   - assign(T &C, const T &A): sets C = A
   *   - mul(T &C, const T &A, const T &B): computes C = AB
   *   - sqr(T &C, const T &A): computes C = AA
   */
  template < class T >
  class ExponentiationBinary : public Exponentiation<T> 
  {

  public:
    ExponentiationBinary() {};
    ~ExponentiationBinary() {};

    /**
     * @brief Computes A^n using left-to-right binary exponentiation
     * @param[out] C result of computing A^n
     * @param[in] A base for exponentiation
     * @param[in] n exponent
     *
     * @pre n >= 0
     */
    void power (T &C, const T &A, const ZZ & n);
  };

} // ANTL

// Unspecialized template definitions.
#include "../../../src/Exponentiation/ExponentiationBinary_impl.hpp"

#endif // EXPONENTIATION_BINARY_H

