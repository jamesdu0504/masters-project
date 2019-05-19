/**
 * @file ExponentiationExtendedDoubleBase_impl.hpp
 * @author Reginald Lybbert
 * @brief generic implementation of templated methods from ExponentiationExtendedDoubleBase class
 */

using namespace ANTL;

//
// initialize class for computing A^n using extended double base representation
//
template < class T >
void 
ExponentiationExtendedDoubleBase<T>::initializeExponent(const ZZ &n)
{

    initializedExponent = n;
    //check if bounds have been explicitly set
    if(!boundsSet){
        powerOfThreeBound = ceil(0.5*NumBits(n));
        powerOfTwoBound = NumBits(n);
    }

    if(powersOfThree.size() != digitSetSize*(powerOfThreeBound + 1)){  //If the powers of three are not yet, or not properly initialized
        powersOfThree.clear();
        ZZ currPowerOfThree = ZZ(1);
        ZZ adjustedPowerOfThree = ZZ(1);
        long currExponent = 0;
        powerOfThreeElem currElem;
        long highestNumBits = ceil(1.585*powerOfThreeBound); // This should be the number of bits of 3^(threeBound).  
    						 // Every element in powersOfThree will have this number of bits
        while(currExponent <= powerOfThreeBound){
           for(int i = 0; i < digitSetSize; i++){
               adjustedPowerOfThree = ((i*3 + 1 + i%2)*currPowerOfThree) << (highestNumBits - NumBits((i*3+1+i%2)*currPowerOfThree));
               currElem.adjustedPower = adjustedPowerOfThree;
               currElem.exponent = currExponent;
               currElem.digitIndex = i;
               powersOfThree.insert(currElem);
           }
           currExponent++;
           currPowerOfThree *= 3;
        }
        maxTerm = (((digitSetSize*3 + 1 - digitSetSize%2)*currPowerOfThree) << powerOfTwoBound)/3;
        
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
    long digitInUpperRep = 0;
    long twosInLowerRep = 0;
    long threesInLowerRep = 0;
    long digitInLowerRep = 0;
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
               digitInLowerRep = (*lower).digitIndex;
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
               digitInUpperRep = (*upper).digitIndex;
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
          currRepElem.digitIndex = digitInLowerRep;
          currRepElem.uniqueID = j++;
          representation.insert(currRepElem);
          ex = ex - bestLowerApprox;
      }else{
          currRepElem.powerOfTwo = twosInUpperRep;
          currRepElem.powerOfThree = threesInUpperRep;
          currRepElem.digitIndex = digitInUpperRep;
          currRepElem.uniqueID = j++;
          representation.insert(currRepElem);
          currRepElem.sign = -currRepElem.sign;
          ex = bestUpperApprox - ex;
   }
  }

/*  //Prints the representation
typename std::set<repElement, comparePowersOfTwo>::iterator printer = representation.begin();
while(printer != representation.end()){
   cout << "(" << (*printer).sign << "," << (*printer).digitIndex << "," << (*printer).powerOfTwo << "," << (*printer).powerOfThree << ")" << endl;
   printer++;
}
*/

}
template < class T >
void
ExponentiationExtendedDoubleBase<T>::initializeBase(const T &A){


//Fix this to take varying digitSetSize, including those greater than 8
   std::vector<T> digitPreComputations;
   digitPreComputations.resize(8);
   assign(digitPreComputations[0], A);  //A^1
   sqr(digitPreComputations[6],A);      //A^2  (temp)
   sqr(digitPreComputations[7],digitPreComputations[6]);       //A^4  (temp)
   mul(digitPreComputations[1],digitPreComputations[7],A);     //A^5
   mul(digitPreComputations[2],digitPreComputations[6],digitPreComputations[1]);  //A^7
   mul(digitPreComputations[3],digitPreComputations[7],digitPreComputations[2]);  //A^11
   mul(digitPreComputations[4],digitPreComputations[6],digitPreComputations[3]);  //A^13
   mul(digitPreComputations[5],digitPreComputations[7],digitPreComputations[4]);  //A^17
   mul(digitPreComputations[6],digitPreComputations[6],digitPreComputations[5]);  //A^19   (overwrite A^2)
   mul(digitPreComputations[7],digitPreComputations[7],digitPreComputations[6]);  //A^23   (overwrite A^4)
   
/*   
  //Print the precomputations
  for(int x = 0; x < digitSetSize; x++){
      cout << "A^" << (3*x+1+x%2) << " = " << digitPreComputations[x] << endl;
  }
*/


   cubes.clear();
   inverseCubes.clear();
   cubes.resize(powerOfThreeBound + 1);
   inverseCubes.resize(powerOfThreeBound + 1);
   cubes[0].resize(digitSetSize);
   inverseCubes[0].resize(digitSetSize);
   //Repeatedly cube A as many times as necessary, and multiply by relevant digits.
   for(long j = 0; j < digitSetSize; j++){
       assign(cubes[0][j], digitPreComputations[j]);
       inv(inverseCubes[0][j],cubes[0][j]);
   }
   for(long i = 0; i < powerOfThreeBound; i++){ 
       cubes[i+1].resize(digitSetSize);
       inverseCubes[i+1].resize(digitSetSize);
       for(long j = 0; j < digitSetSize; j++){
           cub(cubes[i+1][j],cubes[i][j]);
           inv(inverseCubes[i+1][j],cubes[i+1][j]);
       }
   }

}

//
// compute A^n using double base Signed greedy method
//
template < class T > 
void 
ExponentiationExtendedDoubleBase<T>::power (T &C, const T &A, const ZZ & n)
{

   //if statements checking for initialization
   if(n != initializedExponent){
        initializeExponent(n);
   }
   if(cubes.empty()){
        initializeBase(A);
   }else if(!equals(A, cubes[0][0])){
        initializeBase(A);
   }


   typename std::set<repElement>::iterator twoIter = representation.end();

   twoIter--;


   long x = (*twoIter).powerOfTwo;
   assign(C,cubes[(*twoIter).powerOfThree][(*twoIter).digitIndex]);
   if((*twoIter).sign == -1){
      inv(C,C);
   }
   twoIter--;
   while(x > 0){
     if(x == (*twoIter).powerOfTwo){
          if((*twoIter).sign == 1){
             mul(C,C,cubes[(*twoIter).powerOfThree][(*twoIter).digitIndex]);
          }else{
             mul(C,C,inverseCubes[(*twoIter).powerOfThree][(*twoIter).digitIndex]);
          }
          --twoIter;
     }else{
          sqr(C,C);
          x--;
     }
   }
   while(twoIter != representation.end()){
      if((*twoIter).sign == 1){
             mul(C,C,cubes[(*twoIter).powerOfThree][(*twoIter).digitIndex]);
          }else{
             mul(C,C,inverseCubes[(*twoIter).powerOfThree][(*twoIter).digitIndex]);
          }
      twoIter--;
   }

}
