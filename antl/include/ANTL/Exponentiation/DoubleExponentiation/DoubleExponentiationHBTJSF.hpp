/**
 * @file DoubleExponentiationHBTJSF.hpp
 * @author Michael Jacobson/Reginald Lybbert
 * @brief class for hybrid binary-ternary joint sparse form double exponentiation
 * @remarks This algorithm is adapted from the paper "Hybrid Binary-Ternary Joint
 * Sparse Form and its Application in Elliptic Curve Cryptography", Jithra Adikari,
 * Vassil Dimitrov, Laurent Imbert; 2008; https://eprint.iacr.org/2008/285.pdf
 *
 */

#ifndef DOUBLEEXPONENTIATION_HBTJSF_H
#define DOUBLEEXPONENTIATION_HBTJSF_H

#include <vector>
#include "DoubleExponentiation.hpp"

namespace ANTL
{
    

  /**
   * @brief class for hybrid binary-ternat joint sparse form double exponentiation
   * @remarks This concrete class defines a method for hybrid binart-ternary joint sparse form
   * double exponentiation.  The method first computes the hybrid binary-ternary joint sparse form 
   * of the exponents using the basic right-to-left method, and then applies it to the bases
   * left-to-right. The base type is templated, and must have the following
   * functions defined:
   *   - assign(T &C, const T &A): sets C = A
   *   - mul(T &C, const T &A, const T &B): computes C = AB
   *   - sqr(T &C, const T &A): computes C = AA
   *   - cub(T &C, const T &A): computes C = AAA
   *   - inv(T &C, const T &A): computes C = inverse of A
   *   - equals(const T &C, const T &A): returns true if C == A 
   */
  template < class T >
  class DoubleExponentiationHBTJSF : public DoubleExponentiation<T> 
  {

  protected:
    vector<short> em;   /**< vector containing a hybrid binary-ternary expansion of m */
    vector<short> en;   /**< vector containing a hybrid binary-ternary expansion of n */
    vector<short> base; /**< vector containing the bases used for each digit of em and en*/
    ZZ initM;           /**< value currently held in em */
    ZZ initN;           /**< value currently held in en */
    T precomps[16];     /**< array of templated class containing the precomputations for this algorithm */
    
    
    
    void readPrecomps(T &C, short mDigit, short nDigit);  /**<helper function used to read precomps array>*/


  public:
    DoubleExponentiationHBTJSF() {};
    ~DoubleExponentiationHBTJSF() {};

    /**
     * @brief Setup the class to compute A^mB^n.  Computes necessary precomputed values of A and B
     * @param[in] A first base for exponentiation
     * @param[in] B second base for exponentiation
     */
    void initializeBase(const T &A, const T &B);

    /**
     * @brief Setup the class to compute A^mB^n. Computes the HBTJSF representation of m and n
     * @param[in] m first exponent
     * @param[in] n second exponent
     *
     * @pre m,n >= 0
     */
     void initializeExponent(const ZZ &m, const ZZ &n);
   
    /**
     * @brief Computes A^mB^n
     * @param[out] C result of computing A^mB^n using HBTJSF double exponentiation
     * @param[in] A first base for exponentiation
     * @param[in] B second base for exponentiation
     * @param[in] m first exponent
     * @param[in] n second exponent
     *
     * @pre m,n >= 0
     */
    void power (T &C, const T &A, const T &B, const ZZ &m, const ZZ &n);

  };

} // ANTL

// Unspecialized template definitions.
#include "../../../../src/DoubleExponentiation/DoubleExponentiationHBTJSF_impl.hpp"

#endif // DOUBLEEXPONENTIATION_HBTJSF_H

