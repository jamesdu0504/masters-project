/**
 * @file ExponentiationSB3.hpp
 * @author Michael Jacobson/Reginald Lybbert
 * @brief class for left-to-right SB3 exponentiation
 */

#ifndef EXPONENTIATION_SB3_H
#define EXPONENTIATION_SB3_H

#include <vector>
#include "Exponentiation.hpp"

namespace ANTL
{

  /**
   * @brief class for left-to-right SB3 exponentiation
   * @remarks This concrete class defines a method for left-to-right signed base 3
   * exponentiation.  The method first computes the signed ternary representation of the
   * exponent using the basic right-to-left method, and then applies it to the base
   * left-to-right. The base type is templated, and must have the following
   * functions defined:
   *   - assign(T &C, const T &A): sets C = A
   *   - mul(T &C, const T &A, const T &B): computes C = AB
   *   - inv(T &C, const T &A): computes C = inverse of A
   *   - cub(T &C, const T &A): computes C = AAA
   */
  template < class T >
  class ExponentiationSB3 : public Exponentiation<T> 
  {

  protected:
    vector<short> e;  /**< vector containing SB3 expansion of exponent */
    T Ainv;           /**< inverse of the base element */

  public:
    ExponentiationSB3() {};
    ~ExponentiationSB3() {};

    /**
     * @brief Setup the class to compute A^n.  Compute inverse of A and allocates a digit vector
     * of size 2(NumBits(n))/3 + 2
     * @param[in] A base for exponentiation
     * @param[in] n exponent
     *
     * @pre n >= 0
     */
    void initialize(const T &A, const ZZ &n);

    /**
     * @brief Computes A^n.
     * @param[out] C result of computing A^n using left-to-right SB3 exponentiation
     * @param[in] A base for exponentiation
     * @param[in] n exponent
     *
     * @pre n >= 0
     */
    void power (T &C, const T &A, const ZZ &n);

  };

} // ANTL

// Unspecialized template definitions.
#include "../../../src/Exponentiation/ExponentiationSB3_impl.hpp"

#endif // EXPONENTIATION_NAF_H

