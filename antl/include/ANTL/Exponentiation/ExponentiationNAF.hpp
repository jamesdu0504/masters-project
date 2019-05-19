/**
 * @file ExponentiationNAF.hpp
 * @author Michael Jacobson
 * @brief class for left-to-right NAF exponentiation
 */

#ifndef EXPONENTIATION_NAF_H
#define EXPONENTIATION_NAF_H

#include <vector>
#include "Exponentiation.hpp"

namespace ANTL
{

  /**
   * @brief class for left-to-right NAF exponentiation
   * @remarks This concrete class defines a method for left-to-right non-adjacent
   * form exponentiation.  The method first computes the NAF representation of the
   * exponent using the basic right-to-left method, and then applies it to the base
   * left-to-right. The base type is templated, and must have the following
   * functions defined:
   *   - assign(T &C, const T &A): sets C = A
   *   - mul(T &C, const T &A, const T &B): computes C = AB
   *   - sqr(T &C, const T &A): computes C = AA
   *   - inv(T &C, const T &A): computes C = inverse of A
   */
  template < class T >
  class ExponentiationNAF : public Exponentiation<T> 
  {

  protected:
    vector<short> e;  /**< vector containing NAF expansion of exponent */
    T Ainv;           /**< inverse of the base element */

  public:
    ExponentiationNAF() {};
    ~ExponentiationNAF() {};

    /**
     * @brief Setup the class to compute A^n.  Compute inverse of A and allocates a digit vector
     * of size NumBits(n)+1
     * @param[in] A base for exponentiation
     * @param[in] n exponent
     *
     * @pre n >= 0
     */
    void initialize(const T &A, const ZZ &n);

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
#include "../../../src/Exponentiation/ExponentiationNAF_impl.hpp"

#endif // EXPONENTIATION_NAF_H

