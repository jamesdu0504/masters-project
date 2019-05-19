/**
 * @file ExponentiationDoubleBaseStrictChain_impl.hpp
 * @author Reginald Lybbert
 * @brief generic implementation of templated methods from ExponentiationDoubleBaseStrictChain class
 */

using namespace ANTL;

//
// initialize class for computing A^n using double base strict chain
//
template < class T >
void 
ExponentiationDoubleBaseStrictChain<T>::initialize(const T &A)
{
  // compute A^-1
  inv(Ainv,A);

}


//
// compute A^n using double base strict chain method
//
template < class T > 
void 
ExponentiationDoubleBaseStrictChain<T>::power (T &C, const T &A, const ZZ & n)
{
  // compute double base expansion of n (right-to-left)
  ZZ ex = abs (n);
  std::vector<chainElement> e;
    short sign = 1;
    short powerOfTwo;
    short powerOfThree;
    chainElement elem;
    while(ex > 0){
    	powerOfTwo = 0;
    	powerOfThree = 0;
    	while(ex%3 == 0){
    		powerOfThree++;
    		ex = ex/3;
    	}
    	while(ex%2 == 0){
    		powerOfTwo++;
    		ex = ex/2;
    	}
    	if(ex%6 == 5){
    		sign = -1;
    		ex = ex + 1;
    	}
    	else{
    		sign = 1;
    		ex = ex - 1;
    	}
    	elem.sign = sign;
    	elem.powerOfTwo = powerOfTwo;
    	elem.powerOfThree = powerOfThree;
    	e.push_back(elem);
    }

  // compute C = A^n from left-to-right using double base chain elements in e
    assign(C,A);
    elem = e.back();
    e.pop_back();
    while(elem.powerOfThree > 0){
    		cub(C,C);
    		elem.powerOfThree = elem.powerOfThree - 1;
    }
    while(elem.powerOfTwo > 0){
    		sqr(C,C);
    		elem.powerOfTwo = elem.powerOfTwo - 1;
    }
    while(!e.empty()){
    	elem = e.back();
    	e.pop_back();
        if(elem.sign == 1){
           mul(C,C,A);
        }
        else
        {
           mul(C,C,Ainv);
        } 
    	while(elem.powerOfThree > 0){
    		cub(C,C);
    		elem.powerOfThree = elem.powerOfThree - 1;
    	}
    	while(elem.powerOfTwo > 0){
    		sqr(C,C);
    		elem.powerOfTwo = elem.powerOfTwo - 1;
    	}
    }
}
