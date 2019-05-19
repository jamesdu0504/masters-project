/**
 * @file ExponentiationDoubleBaseUnsignedGreedy.hpp
 * @author Reginald Lybbert
 * @brief class for double base exponentiation using the unsigned greedy (2,3) representation
 */

#ifndef EXPONENTIATION_DoubleBaseUnsignedGreedy_H
#define EXPONENTIATION_DoubleBaseUnsignedGreedy_H

#include <vector>
#include <set>
#include <math.h>
#include "Exponentiation.hpp"

namespace ANTL
{

  /**
   * @brief class for double base exponentiation using the unsigned greedy (2,3) representation
   * @remarks This concrete class defines a method for double base exponentiation using the 
   * unsigned greedy (2,3) representation  The method first computes the double base representation 
   * of the exponent using a greedy algorithm, and then applies it to the base
   * The base type is templated, and must have the following
   * functions defined:
   *   - assign(T &C, const T &A): sets C = A
   *   - mul(T &C, const T &A, const T &B): computes C = AB
   *   - sqr(T &C, const T &A): computes C = AA
   *   - cub(T &C, const T &A): computes C = AAA
   */
  template < class T >
  class ExponentiationDoubleBaseUnsignedGreedy : public Exponentiation<T> 
  {

  protected:
    
    struct powerOfThreeElem {
	ZZ adjustedPower;
	long exponent;
    };
    
    struct comparePowers
    {
       bool operator()(const powerOfThreeElem& x, const powerOfThreeElem& y){ return x.adjustedPower < y.adjustedPower; }
    };

    struct repElement {
        long powerOfTwo;
	long powerOfThree;
        long uniqueID;
       
    };
    struct comparePowersOfTwo
    {
      bool operator()(const repElement& x, const repElement& y){
        if(x.powerOfTwo == y.powerOfTwo){ return x.uniqueID < y.uniqueID; }
        else{ return x.powerOfTwo < y.powerOfTwo; }
      }
    };


    std::set<powerOfThreeElem, comparePowers> powersOfThree;
    long powerOfTwoBound;
    long powerOfThreeBound;

  public:
    ExponentiationDoubleBaseUnsignedGreedy() {};
    ~ExponentiationDoubleBaseUnsignedGreedy() {};

    /**
     * @brief Setup the class to compute A^n.  Compute powerOfThree and set given bounds.
     * @param[in] twoBound - the upper bound on powers of two in the representation
     * @param[in] threeBound - the upper bound on powers of three in this representation

     * @pre twoBound > 0 
     * @pre threeBound > 0

     */
    void initialize(long twoBound, long threeBound);

    /**
     * @brief Computes A^n.
     * @param[out] C result of computing A^n using double base unsigned greedy
     * @param[in] A base for exponentiation
     * @param[in] n exponent
     *
     * @pre n >= 0
     */
    void power (T &C, const T &A, const ZZ &n);

  };

} // ANTL

// Unspecialized template definitions.
#include "../../../src/Exponentiation/ExponentiationDoubleBaseUnsignedGreedy_impl.hpp"

#endif // EXPONENTIATION_DoubleBaseUnsignedGreedy_H

