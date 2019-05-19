/**
 * @file ExponentiationDoubleBaseUnsignedGreedy_impl.hpp
 * @author Reginald Lybbert
 * @brief generic implementation of templated methods from ExponentiationUnsignedGreedy class
 */

using namespace ANTL;



//
// initialize class for computing A^n using double base strict chain
//
template < class T >
void 
ExponentiationDoubleBaseUnsignedGreedy<T>::initialize(long twoBound, long threeBound)
{
    powerOfTwoBound = twoBound;
    powerOfThreeBound = threeBound;
    ZZ currPowerOfThree = ZZ(1);
    ZZ adjustedPowerOfThree = ZZ(1);
    long currExponent = 0;
    powerOfThreeElem currElem;
    long highestNumBits = ceil(1.585*threeBound); // This should be the number of bits of 3^(threeBound).  
						 // Every element in powersOfThree will have this number of bits
    while(currExponent <= threeBound){
       adjustedPowerOfThree = currPowerOfThree << (highestNumBits - NumBits(currPowerOfThree));
       currElem.adjustedPower = adjustedPowerOfThree;
       currElem.exponent = currExponent;
       powersOfThree.insert(currElem);
       currExponent++;
       currPowerOfThree *= 3;
    }

}

//
// compute A^n using double base unsigned greedy method
//
template < class T > 
void 
ExponentiationDoubleBaseUnsignedGreedy<T>::power (T &C, const T &A, const ZZ & n)
{

  ZZ ex = abs(n);
  ZZ bestApprox;
  std::set<repElement, comparePowersOfTwo> representation;
  repElement currRepElem;
  powerOfThreeElem testElem;
  long twoOffset;
  long highestNumBits = ceil(1.585*powerOfThreeBound);
  bool foundFlag = false;
  long twosInRep = 0;
  long threesInRep = 0;
  int i;
  long j = 0;


    //Compute unsigned greedy (2,3) representation of n
  while(ex > 0){
      twoOffset = highestNumBits - NumBits(ex);
      testElem.adjustedPower = ex << twoOffset;

        //use upper_bound to find the location in the set for ex
      typename std::set<powerOfThreeElem>::iterator upper;
      upper = powersOfThree.upper_bound(testElem);

      //use this to find the best approximation to ex
      i = 0;
      foundFlag = false;
      while(!foundFlag){
        while(upper != powersOfThree.begin() && !foundFlag){
          upper--;
          twosInRep = NumTwos((*upper).adjustedPower) - twoOffset - i; //3^(upper.exponent)*2^(twosInRep) is first term in greedy rep of ex
          if(twosInRep >= 0 && twosInRep <= powerOfTwoBound){ 
               foundFlag = true; 
               threesInRep = (*upper).exponent;
               bestApprox = (*upper).adjustedPower >> (twoOffset + i); 
          }
        }
        upper = powersOfThree.end();
        i++;
      }
      
        //ex = ex - bestApprox, repeat until ex is 0
      currRepElem.powerOfTwo = twosInRep;
      currRepElem.powerOfThree = threesInRep;
      currRepElem.uniqueID = j++;
      representation.insert(currRepElem);
      ex = ex - bestApprox;
   }

   T cubes[powerOfThreeBound + 1];
  //Compute C = A^n using that representation
        //Repeatedly cube A as many times as necessary
   assign(cubes[0], A);
   for(long i = 0; i < powerOfThreeBound; i++){ 
       cub(cubes[i+1],cubes[i]);
   }
        //From largest power of two down, multiply in the powers of three, squaring in between 


   typename std::set<repElement>::iterator twoIter = representation.end();
   twoIter--;

   long x = (*twoIter).powerOfTwo;
   assign(C,cubes[(*twoIter).powerOfThree]);
   twoIter--;
   while(x > 0){
     if(x == (*twoIter).powerOfTwo){
          mul(C,C,cubes[(*twoIter).powerOfThree]);
          --twoIter;
     }else{
          sqr(C,C);
          x--;
     }
   }
   while(twoIter != representation.end()){
      mul(C,C,cubes[(*twoIter).powerOfThree]);
      twoIter--;
   }

}
