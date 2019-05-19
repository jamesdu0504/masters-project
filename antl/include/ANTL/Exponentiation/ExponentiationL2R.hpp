/**
 * @file ExponentiationL2R.hpp
 * @author Michael Jacobson/Reginald Lybbert
 * @brief class for left-to-right NAF exponentiation with NAF representation also computed left-to-right
 */

#ifndef EXPONENTIATION_L2R_H
#define EXPONENTIATION_L2R_H

#include <vector>
#include "Exponentiation.hpp"

namespace ANTL
{

  /**
   * @brief class for left-to-right NAF exponentiation with NAF representation also computed left-to-right
   * @remarks This concrete class defines a method for left-to-right non-adjacent
   * form exponentiation.  The method computes the NAF representation of the
   * exponent using a left-to-right method, meanwhile applying it to the base
   * left-to-right. The base type is templated, and must have the following
   * functions defined:
   *   - assign(T &C, const T &A): sets C = A
   *   - mul(T &C, const T &A, const T &B): computes C = AB
   *   - sqr(T &C, const T &A): computes C = AA
   *   - inv(T &C, const T &A): computes C = inverse of A
   */
  template < class T >
  class ExponentiationL2R : public Exponentiation<T> 
  {

  protected:
    T Ainv;           /**< inverse of the base element */

  public:
    ExponentiationL2R() {};
    ~ExponentiationL2R() {};

    /**
     * @brief Setup the class to compute A^n.  Compute inverse of A.
     * @param[in] A base for exponentiation
     */
    void initialize(const T &A);

    /**
     * @brief Computes A^n.
     * @param[out] C result of computing A^n using left-to-right NAF exponentiation
     * @param[in] A base for exponentiation
     * @param[in] n exponent
     *
     * @pre n >= 0
     */
    void power (T &C, const T &A, const ZZ &n);

  };

} // ANTL

// Unspecialized template definitions.
#include "../../../src/Exponentiation/ExponentiationL2R_impl.hpp"

#endif // EXPONENTIATION_L2R_H

