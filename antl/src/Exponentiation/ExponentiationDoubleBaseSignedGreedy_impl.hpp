/**
 * @file ExponentiationDoubleBaseSignedGreedy_impl.hpp
 * @author Reginald Lybbert
 * @brief generic implementation of templated methods from ExponentiationSignedGreedy class
 */

using namespace ANTL;



//
// initialize class for computing A^n using double base strict chain
//
template < class T >
void 
ExponentiationDoubleBaseSignedGreedy<T>::initializeExponent(const ZZ &n)
{

    initializedExponent = n;
    //check if bounds have been explicitly set
    if(!boundsSet){
        powerOfThreeBound = ceil(0.5*NumBits(n));
        powerOfTwoBound = NumBits(n);
    }

    if(powersOfThree.size() != powerOfThreeBound + 1){  //If the powers of three are not yet, or not properly initialized
        powersOfThree.clear();
        ZZ currPowerOfThree = ZZ(1);
        ZZ adjustedPowerOfThree = ZZ(1);
        long currExponent = 0;
        powerOfThreeElem currElem;
        long highestNumBits = ceil(1.585*powerOfThreeBound); // This should be the number of bits of 3^(threeBound).  
    						 // Every element in powersOfThree will have this number of bits
        while(currExponent <= powerOfThreeBound){
           adjustedPowerOfThree = currPowerOfThree << (highestNumBits - NumBits(currPowerOfThree));
           currElem.adjustedPower = adjustedPowerOfThree;
           currElem.exponent = currExponent;
           powersOfThree.insert(currElem);
           currExponent++;
           currPowerOfThree *= 3;
        }
        maxTerm = (currPowerOfThree << powerOfTwoBound)/3;
        
    }

    ZZ ex = abs(n);
    ZZ bestUpperApprox;
    ZZ bestLowerApprox;
    repElement currRepElem;
    currRepElem.sign = 1;
    powerOfThreeElem testElem;
    long twoOffset;
    long highestNumBits = ceil(1.585*powerOfThreeBound);
    bool foundFlag = false;
    long twosInUpperRep = 0;
    long threesInUpperRep = 0;
    long twosInLowerRep = 0;
    long threesInLowerRep = 0;
    int i;
    long j = 0;


  //Compute Signed greedy (2,3) representation of n
    while(ex > 0){
      twoOffset = highestNumBits - NumBits(ex);
      testElem.adjustedPower = ex << twoOffset;

        //use upper_bound to find the location in the set for ex
      typename std::set<powerOfThreeElem>::iterator upper, lower;
      upper = powersOfThree.upper_bound(testElem);
      lower = powersOfThree.upper_bound(testElem);

      //use this to find the best approximation to ex less than ex.
      i = 0;
      foundFlag = false;
      while(!foundFlag){
        while(lower != powersOfThree.begin() && !foundFlag){
          lower--;
          
          twosInLowerRep = NumTwos((*lower).adjustedPower) - twoOffset - i; //3^(upper.exponent)*2^(twosInRep) is first term in greedy rep of ex
        if(twosInLowerRep >= 0 && twosInLowerRep <= powerOfTwoBound){ 
               foundFlag = true; 
               threesInLowerRep = (*lower).exponent;
               bestLowerApprox = (*lower).adjustedPower >> (twoOffset + i); 
          }
        }
        lower = powersOfThree.end();
        i++;
      }

      if(ex < maxTerm){
         //use this to find the best approximation to ex greater than ex.
         i = 0;
         foundFlag = false;
         while(!foundFlag){
           while(upper != powersOfThree.end() && !foundFlag){
             twosInUpperRep = NumTwos((*upper).adjustedPower) - twoOffset + i; //3^(upper.exponent)*2^(twosInRep) is first term in greedy rep of ex
             if(twosInUpperRep >= 0 && twosInUpperRep <= powerOfTwoBound){ 
               foundFlag = true; 
               threesInUpperRep = (*upper).exponent;
               bestUpperApprox = (*upper).adjustedPower >> (twoOffset - i); 
             }
            upper++;
          }
          upper = powersOfThree.begin();
          i++;
        }
      }else{
             bestUpperApprox = ex + bestLowerApprox + 1;
      }
      
      //ex = ex - bestApprox, repeat until ex is 0

      if(ex - bestLowerApprox <= bestUpperApprox - ex){
          currRepElem.powerOfTwo = twosInLowerRep;
          currRepElem.powerOfThree = threesInLowerRep;
          currRepElem.uniqueID = j++;
          representation.insert(currRepElem);
          ex = ex - bestLowerApprox;
      }else{
          currRepElem.powerOfTwo = twosInUpperRep;
          currRepElem.powerOfThree = threesInUpperRep;
          currRepElem.uniqueID = j++;
          representation.insert(currRepElem);
          currRepElem.sign = -currRepElem.sign;
          ex = bestUpperApprox - ex;
   }
  }


}
template < class T >
void
ExponentiationDoubleBaseSignedGreedy<T>::initializeBase(const T &A){



   cubes.clear();
   inverseCubes.clear();
   cubes.resize(powerOfThreeBound + 1);
   inverseCubes.resize(powerOfThreeBound + 1);
   //Repeatedly cube A as many times as necessary
   assign(cubes[0], A);
   inv(inverseCubes[0],A);
   for(long i = 0; i < powerOfThreeBound; i++){ 
       cub(cubes[i+1],cubes[i]);
       inv(inverseCubes[i+1],cubes[i+1]);
   }

}

//
// compute A^n using double base Signed greedy method
//
template < class T > 
void 
ExponentiationDoubleBaseSignedGreedy<T>::power (T &C, const T &A, const ZZ & n)
{

   //if statements checking for initialization
   if(n != initializedExponent){
        initializeExponent(n);
   }
   if(cubes.empty()){
        initializeBase(A);
   }else if(!equals(A, cubes[0])){
        initializeBase(A);
   }


   typename std::set<repElement>::iterator twoIter = representation.end();

   twoIter--;


   long x = (*twoIter).powerOfTwo;
   assign(C,cubes[(*twoIter).powerOfThree]);
   if((*twoIter).sign == -1){
      inv(C,C);
   }
   twoIter--;
   while(x > 0){
     if(x == (*twoIter).powerOfTwo){
          if((*twoIter).sign == 1){
             mul(C,C,cubes[(*twoIter).powerOfThree]);
          }else{
             mul(C,C,inverseCubes[(*twoIter).powerOfThree]);
          }
          --twoIter;
     }else{
          sqr(C,C);
          x--;
     }
   }
   while(twoIter != representation.end()){
      if((*twoIter).sign == 1){
             mul(C,C,cubes[(*twoIter).powerOfThree]);
          }else{
             mul(C,C,inverseCubes[(*twoIter).powerOfThree]);
          }
      twoIter--;
   }

}
