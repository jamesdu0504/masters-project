/**
 * @file ExponentiationDoubleBaseStrictChain.hpp
 * @author Reginald Lybbert
 * @brief class for double base exponentiation using the strict chain (2,3) representation
 */

#ifndef EXPONENTIATION_DoubleBaseStrictChain_H
#define EXPONENTIATION_DoubleBaseStrictChain_H

#include <vector>
#include "Exponentiation.hpp"

namespace ANTL
{

  /**
   * @brief class for double base exponentiation using the strict chain (2,3) representation
   * @remarks This concrete class defines a method for double base exponentiation using the 
   * strict chain (2,3) representation  The method first computes the double base representation 
   * of the exponent using the basic right-to-left method, and then applies it to the base
   * left-to-right. The base type is templated, and must have the following
   * functions defined:
   *   - assign(T &C, const T &A): sets C = A
   *   - mul(T &C, const T &A, const T &B): computes C = AB
   *   - sqr(T &C, const T &A): computes C = AA
   *   - cub(T &C, const T &A): computes C = AAA
   *   - inv(T &C, const T &A): computes C = inverse of A
   */
  template < class T >
  class ExponentiationDoubleBaseStrictChain : public Exponentiation<T> 
  {

  protected:
    
    T Ainv;           /**< inverse of the base element */

    struct chainElement{
          short sign;
          short powerOfTwo;
          short powerOfThree;
    };

  public:
    ExponentiationDoubleBaseStrictChain() {};
    ~ExponentiationDoubleBaseStrictChain() {};

    /**
     * @brief Setup the class to compute A^n.  Compute inverse of A
     * @param[in] A base for exponentiation

     */
    void initialize(const T &A);

    /**
     * @brief Computes A^n.
     * @param[out] C result of computing A^n using double base strict chain
     * @param[in] A base for exponentiation
     * @param[in] n exponent
     *
     * @pre n >= 0
     */
    void power (T &C, const T &A, const ZZ &n);

  };

} // ANTL

// Unspecialized template definitions.
#include "../../../src/Exponentiation/ExponentiationDoubleBaseStrictChain_impl.hpp"

#endif // EXPONENTIATION_DoubleBaseStrictChain_H

