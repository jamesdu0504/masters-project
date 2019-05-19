/**
 * @file MultiExponentiationPippenger.hpp
 * @author Michael Jacobson/Reginald Lybbert
 * @brief class for Pippenger's Multi-exponentiation Algorithm, as described in
 * "Pippenger's Multiproduct and Multiexponentiation Algorithms"; Ryan Henry; 2010
 * 
 *   cacr.uwaterloo.ca/techreports/2010/cacr2010-26.pdf
 */

#ifndef MULTIEXPONENTIATION_PIPPENGER_H
#define MULTIEXPONENTIATION_PIPPENGER_H

#include <ANTL/MultiExponentiation/MultiExponentiation.hpp>

namespace ANTL
{

  /**
   * @brief class for Pippenger's multi-exponentiation algorithm
   * @remarks This concrete class defines a method for Pippenger's multi-exponentiation
   * algorithm.
   *   - assign(T &C, const T &A): sets C = A
   *   - mul(T &C, const T &A, const T &B): computes C = AB
   *   - sqr(T &C, const T &A): computes C = AA
   *   - id(T &C): sets C = 1  (or whatever the identity is in the group T represents)
   */
  template < class T >
  class MultiExponentiationPippenger : public MultiExponentiation<T> 
  {

  protected:
      /* @param[out] x_prime: the set of new inputs
       * @param[out] y_prime: the set of new outputs
       * @param[in] x: the set of inputs
       * @param[in] y: the set of desired outputs
       * @param[in] r: the radix
       * @param[in] b: the word-length
       *
       * @pre r*b >= lg(y[i][j] + 1) for all i,j 
       *
       */
    void decompose( vector<T> &x_prime, vector<vector<long> > &y_prime,
                    const vector<T> &x, const vector<vector<ZZ> > &y, const long r, const long b );

      /* @param[out] C: the set of outputs
       * @param[in] x: the set of inputs
       * @param[in] y: the set of desired outputs, expressed in terms of indices of x
       */
    void combine( vector<T> &C, const vector<T> &x, const long &r);

      /* @param[out] C: the set of outputs
       * @param[in]  x: the set of inputs
       * @param[in] r:  the radix
       */

    void multiprod( vector<T> &C, const vector<T> &x, const vector<vector<long> > &y); 

      /* @param[out] ell: the iteration limit
       * @param[out] c: the initial clumping factor
       * @param[out] alpha: the clumping factors
       * @param[out] beta: the grouping factors
       * @param[in] x: the set of inputs
       * @param[in] y: the set of desired outputs
       *
       */
    void getParams( long &ell, long &c, vector<long> &alpha, vector<long> &beta, const vector<T> &x, const vector<vector<long> > &y);

      /* @param[out] C: the set of outputs
       * @param[in] i: the current iteration index
       * @param[in] ell: the iteration limit
       * @param[in] c: the initial clumping factor
       * @param[in] alpha: the clumping factors
       * @param[in] beta: the grouping factors
       * @param[in] x: the set of inputs
       * @param[in] y: the set of desired outputs
       *
       */
    void computeMultiProd( vector<T> &C, const long i, const long &ell, const long &c, const vector<long> &alpha, 
                            const vector<long> &beta, const vector<T> &x, const vector<vector<long> > &y);

      /* @param[out] x_prime: the set of new inputs
       * @param[out] y_prime: the set of new outputs
       * @param[in] x: the set of inputs
       * @param[in] y: the set of desired outputs
       * @param[in] c: the initial clumping factor
       *
       */
    void inputPartition(vector<T> &x_prime, vector<vector<long> > &y_prime, const vector<T> &x,
                          const vector<vector<long> > &y, const long &c);
 
      /* @param[out] y_prime: the set of new outputs
       * @param[out] y_doubleprime: the set of original outputs, expressed in terms of both new and old outputs
       * @param[in] x: the set of inputs
       * @param[in] y: the set of desired outputs
       * @param[in] beta: the grouping factor
       * @param[in] alpha: the clumping factor
       *
       */
    void outputClump(vector<vector<long> > &y_prime, vector<vector<long> > &y_doubleprime, 
                        const vector<T> &x, const vector<vector<long> > &y, const long alpha, const long beta);
  
      /* @param[out] x_prime: the set of new inputs
       * @param[out] y_prime: the set of new outputs
       * @param[in] x: the set of inputs
       * @param[in] y: the set of desired outputs
       * @param[in] beta: the grouping factor
       * @param[in] alpha: the clumping factor
       *
       */
    void inputClump(vector<T> &x_prime, vector<vector<long> > &y_prime, const vector<T> &x,
                        const vector<vector<long> > &y, const long alpha, const long beta);
    

      /* @param[out] C: the set of outputs
       * @param[in] x: the set of inputs
       * @param[in] y: the set of desired outputs, expressed in terms of indices of x
       */
    void naiveMultiply(vector<T> &C, const vector<T> &x, const vector<vector<long> > &y);

  public:
    MultiExponentiationPippenger() {};
    ~MultiExponentiationPippenger() {};


    /* @brief Computes C[j] = A[1]^n[j][1]*A[2]^n[j][2]*...*A[i]^n[j][i] using Pippenger's algorithm
     * @param[out] C vector of results of computing A[1]^n[j][1]*A[2]^n[j][2]*...*A[i]^n[j][i].
     * @param[in] A vector of bases for exponentiation
     * @param[in] n vector of exponents
     *
     * @pre len(A) = len(n)
     * @pre n[i] > 0 for all i
     */
    void power(vector<T> &C, const vector<T> &A, const vector<vector<ZZ> > &n);

  };

} // ANTL

// Unspecialized template definitions.
#include "../../../src/MultiExponentiation/MultiExponentiationPippenger_impl.hpp"

#endif // MULTIEXPONENTIATION_PIPPENGER_H

