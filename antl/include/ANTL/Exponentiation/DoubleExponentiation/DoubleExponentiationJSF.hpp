/**
 * @file DoubleExponentiationJSF.hpp
 * @author Michael Jacobson/Reginald Lybbert
 * @brief class for left-to-right JSF double exponentiation
 */

#ifndef DOUBLEEXPONENTIATION_JSF_H
#define DOUBLEEXPONENTIATION_JSF_H

#include <vector>
#include "DoubleExponentiation.hpp"

namespace ANTL
{

  /**
   * @brief class for left-to-right JSF double exponentiation
   * @remarks This concrete class defines a method for left-to-right joint sparse form
   * double exponentiation.  The method first computes the joint sparse form of the
   * exponents using the basic right-to-left method, and then applies it to the bases
   * left-to-right. The base type is templated, and must have the following
   * functions defined:
   *   - assign(T &C, const T &A): sets C = A
   *   - mul(T &C, const T &A, const T &B): computes C = AB
   *   - sqr(T &C, const T &A): computes C = AA
   *   - inv(T &C, const T &A): computes C = inverse of A
   */
  template < class T >
  class DoubleExponentiationJSF : public DoubleExponentiation<T> 
  {

  protected:
    vector<short> em;  /**< vector containing a NAF expansion of m */
    vector<short> en;  /**< vector containing a NAF expansion of n */
    T AInv;	       /**< templated class containing the inverse of A */
    T BInv;	       /**< templated class containing the inverse of B */
    T AB; 	       /**< templated class containing the product of A and B */
    T ABInv;           /**< templated class containing the product of A and BInv */
    T AInvB;           /**< templated class containing the product of AInv and B */
    T AInvBInv;        /**< templated class containing the product of AInv and BInv */

  public:
    DoubleExponentiationJSF() {};
    ~DoubleExponentiationJSF() {};

    /**
     * @brief Setup the class to compute A^mB^n.  Computes necessary precomputed values of A and B and allocates a digit vector
     * of size NumBits(n)+1
     * @param[in] A first base for exponentiation
     * @param[in] B second base for exponentiation
     * @param[in] m first exponent
     * @param[in] n second exponent
     *
     * @pre m,n >= 0
     */
    void initialize(const T &A, const T &B, const ZZ &m, const ZZ &n);

    /**
     * @brief Computes A^mB^n
     * @param[out] C result of computing A^m^n using left-to-right JSF double exponentiation
     * @param[in] A first base for exponentiation
     * @param[in] B second base for exponentiation
     * @param[in] m first exponent
     * @param[in] n second exponent
     *
     * @pre n >= 0
     */
    void power (T &C, const T &A, const T &B, const ZZ &m, const ZZ &n);

  };

} // ANTL

// Unspecialized template definitions.
#include "../../../../src/DoubleExponentiation/DoubleExponentiationJSF_impl.hpp"

#endif // DOUBLEEXPONENTIATION_JSF_H

