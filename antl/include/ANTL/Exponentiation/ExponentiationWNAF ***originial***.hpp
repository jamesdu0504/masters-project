/**
 * @file ExponentiationWNAF.hpp
 * @author Michael Jacobson/Reginald Lybbert
 * @brief class for left-to-right WNAF exponentiation
 */

#ifndef EXPONENTIATION_WNAF_H
#define EXPONENTIATION_WNAF_H

#include <vector>
#include "Exponentiation.hpp"

namespace ANTL
{

  /**
   * @brief class for left-to-right WNAF exponentiation
   * @remarks This concrete class defines a method for left-to-right windowed non-adjacent
   * form exponentiation.  The method first computes the WNAF representation of the
   * exponent using the basic right-to-left method, and then applies it to the base
   * left-to-right. The base type is templated, and must have the following
   * functions defined:
   *   - assign(T &C, const T &A): sets C = A
   *   - mul(T &C, const T &A, const T &B): computes C = AB
   *   - sqr(T &C, const T &A): computes C = AA
   *   - inv(T &C, const T &A): computes C = inverse of A
   */
  template < class T >
  class ExponentiationWNAF : public Exponentiation<T> 
  {

  protected:
    vector<short> e;  /**< vector containing WNAF expansion of exponent */
    T * precomp;      /**< array containing precomputed values used for WNAF */
    T * precompInv;   /**< array containing the inverses of the values in precomp */ 
    short w;	      /**< width of WNAF representation */

  public:
    ExponentiationWNAF() {};
    ~ExponentiationWNAF()
    {
        
    };

    /**
     * @brief Setup the class to compute A^n.  Computes necessary precomputed values of A and allocates a digit vector
     * of size NumBits(n)+1
     * @param[in] A base for exponentiation
     * @param[in] n exponent
     * @param[in] w windowing width
     *
     * @pre n >= 0
     * @pre w >= 2
     */
    void initialize(const T &A, const ZZ &n, const short w);

    /**
     * @brief Computes A^n.
     * @param[out] C result of computing A^n using left-to-right WNAF exponentiation
     * @param[in] A base for exponentiation
     * @param[in] n exponent
     *
     * @pre n >= 0
     */
    void power (T &C, const T &A, const ZZ &n);

  };

} // ANTL

// Unspecialized template definitions.
#include "../../../src/Exponentiation/ExponentiationWNAF_impl.hpp"

#endif // EXPONENTIATION_WNAF_H

