/**
 * @file MultiExponentiation.hpp
 * @author Michael Jacobson/Reginald Lybbert
 * @brief virtual superclass for generic multi exponentiation
 */

#ifndef MULTIEXPONENTIATION_H
#define MULTIEXPONENTIATION_H

#include <ANTL/common.hpp>
#include <vector>

namespace ANTL 
{

  /**
   * @brief Virtual superclass for generic multi-exponentiation
   * @remarks This virtual class defines a method for generic multi-exponentiation.
   * Particular multi-exponentiation algorithms are defined and implemented in concrete
   * subclasses.  The base type is templated, and must have as a minimum the following
   * functions defined:
   *   - assign(T &C, const T &A): sets C = A
   *   - mul(T &C, const T &A, const T &B): computes C = AB
   *   - sqr(T &C, const T &A): computes C = AA
   * Concrete subclasses may have further requirements (eg. a cube function) that are
   * described in the documentation of each concrete class. 
   */
  template < class T >
  class MultiExponentiation 
  {

  public:
    MultiExponentiation() {};
    virtual ~MultiExponentiation() {};


    /**
     * @brief Computes A[1]^n[j][1]*A[2]^n[j][2]*...*A[i]^n[j][i].  Virtual method, to be overridden 
     * by concrete classes corresponding to different exponentiation methods.
     * @param[out] C vector of results of computing A[1]^n[j][1]*A[2]^n[j][2]*...*A[i]^n[j][i].
     * @param[in] A vector of bases for exponentiation
     * @param[in] n vector of vectors of exponents
     *
     * @pre len(C) = len(n) (number of equations)
     * @pre len(A) = len(n[j]) for all j  (number of variables)
     * @pre n[j][i] > 0 for all i
     */
    virtual void power (vector<T> &C, const vector<T> &A, const vector<vector<ZZ> > &n) = 0;
  };

} // ANTL

#endif // MULTIEXPONENTIATION_H

