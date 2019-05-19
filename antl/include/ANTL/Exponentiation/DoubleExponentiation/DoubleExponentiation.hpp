/**
 * @file DoubleExponentiation.hpp
 * @author Michael Jacobson/Reginald Lybbert
 * @brief virtual superclass for generic double exponentiation
 */

#ifndef DOUBLEEXPONENTIATION_H
#define DOUBLEEXPONENTIATION_H

#include "../../common.hpp"

namespace ANTL 
{

  /**
   * @brief Virtual superclass for generic double exponentiation
   * @remarks This virtual class defines a method for generic  double exponentiation.
   * Particular double exponentiation algorithms are defined and implemented in concrete
   * subclasses.  The base type is templated, and must have as a minimum the following
   * functions defined:
   *   - assign(T &C, const T &A): sets C = A
   *   - mul(T &C, const T &A, const T &B): computes C = AB
   *   - sqr(T &C, const T &A): computes C = AA
   * Concrete subclasses may have further requirements (eg. a cube function) that are
   * described in the documentation of each concrete class. 
   */
  template < class T >
  class DoubleExponentiation 
  {

  public:
    DoubleExponentiation() {};
    virtual ~DoubleExponentiation() {};


    /**
     * @brief Computes A^mB^n.  Virtual method, to be overridden by concrete classes corresponding
     * to different exponentiation methods.
     * @param[out] C result of computing A^mB^n
     * @param[in] A first base for exponentiation
     * @param[in] B second base for exponentiation
     * @param[in] m first exponent
     * @param[in] n second exponent
     *
     * @pre m,n >= 0
     */
    virtual void power (T &C, const T &A, const T &B, const ZZ & m, const ZZ & n) = 0;
  };

} // ANTL

#endif // DOUBLEEXPONENTIATION_H

