/**
 * @file ExponentiationExtendedDoubleBase.hpp
 * @author Reginald Lybbert
 * @brief class for double base exponentiation using the Signed greedy (2,3) representation, with an extended digit set.
 */

#ifndef EXPONENTIATION_ExtendedDoubleBase_H
#define EXPONENTIATION_ExtendedDoubleBase_H

#include <vector>
#include <set>
#include <math.h>
#include "Exponentiation.hpp"

namespace ANTL
{

  /**
   * @brief class for double base exponentiation using the extended signed greedy (2,3) representation
   * @remarks This concrete class defines a method for double base exponentiation using the extended
   * signed greedy (2,3) representation  The method first computes the double base representation 
   * of the exponent using a greedy algorithm, and then applies it to the base
   * The base type is templated, and must have the following
   * functions defined:
   *   - assign(T &C, const T &A): sets C = A
   *   - mul(T &C, const T &A, const T &B): computes C = AB
   *   - sqr(T &C, const T &A): computes C = AA
   *   - cub(T &C, const T &A): computes C = AAA
   *   - inv(T &C, const T &A): computes C = inverse of A
   *   - equals(const T &C, const T &A): returns true of C == A, false otherwise. 
   */
  template < class T >
  class ExponentiationExtendedDoubleBase : public Exponentiation<T> 
  {

  protected:
    
    long digitSetSize;   //How many elements in the sequence 1,5,7,11,13,17,19,23,25,29,... to include as digits (non-multiples of 2 or 3).

    struct powerOfThreeElem {
	ZZ adjustedPower;
	long exponent;
        long digitIndex;
    };
    
    struct comparePowers
    {
       bool operator()(const powerOfThreeElem& x, const powerOfThreeElem& y){ return x.adjustedPower < y.adjustedPower; }
    };

    struct repElement {
        long powerOfTwo;
	long powerOfThree;
        long digitIndex;
        long uniqueID;
        short sign;  //should only take values 1,-1
       
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
    ZZ maxTerm;
    bool boundsSet;
 
    std::set<repElement, comparePowersOfTwo> representation;
    ZZ initializedExponent;  //the exponent whose representation is currently stored in representation

    std::vector<std::vector<T> > cubes;
    std::vector<std::vector<T> > inverseCubes;

  public:
    ExponentiationExtendedDoubleBase() { boundsSet = false; digitSetSize = 4; };
    ~ExponentiationExtendedDoubleBase() {};


    /**
    * @brief Optional constructor, sets the upper bounds for powers of two and three allowed in this representation
    * @param[in] digitSize - the number of digits to include in the extended digit set.
    * @param[in] twoBound - the upper bound on powers of two in the representation
    * @param[in] threeBound - the upper bound on powers of three in the representation
    *
    * @pre digitSetSize > 0;
    * @pre twoBound > 0 
    * @pre threeBound > 0
    *
    */
    ExponentiationExtendedDoubleBase(long digitSize, long twoBound, long threeBound) {
          digitSetSize = digitSize;
          powerOfTwoBound = twoBound; 
          powerOfThreeBound = threeBound;
          boundsSet = true;
    }

    /**
    * @brief Sets the upper bounds for powers of two and three allowed in this representation
    * @param[in] twoBound - the upper bound on powers of two in the representation
    * @param[in] threeBound - the upper bound on powers of three in the representation
    *
    * @pre twoBound > 0 
    * @pre threeBound > 0
    *
    */
    void setBounds(long twoBound, long threeBound) {
          powerOfTwoBound = twoBound; 
          powerOfThreeBound = threeBound; 
    }


    /**
     * @brief Setup the class to compute A^n.  Compute cubes of A up to powerOfThreeBound.
     * @param[in] A - base for exponentiation.
     *
     */
    void initializeBase(const T &A);

   /**
     * @brief Setup the class to compute A^n.  Compute representation of n.
     * @param[in] n - exponent.
     *
     */
    void initializeExponent(const ZZ &n);

    /**
     * @brief Computes A^n.
     * @param[out] C result of computing A^n using extended double base
     * @param[in] A base for exponentiation
     * @param[in] n exponent
     *
     * @pre n >= 0
     */
    void power (T &C, const T &A, const ZZ &n);

  };

} // ANTL

// Unspecialized template definitions.
#include "../../../src/Exponentiation/ExponentiationExtendedDoubleBase_impl.hpp"

#endif // EXPONENTIATION_ExtendedDoubleBase_H

