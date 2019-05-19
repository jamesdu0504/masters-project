/**
 * @file zz_p_expo_test.cpp
 * @author Michael Jacobson
 * @brief Test program for exponentiation routines using zz_p as base type.
 */

#include <NTL/lzz_p.h>

inline void cub(NTL::zz_p& x, NTL::zz_p a)
{
     sqr(x,a);
     mul(x,x,a);
}

inline bool equals(const NTL::zz_p& x,const NTL::zz_p& a)
{
    return x == a;

}  

#include "../../include/ANTL/Exponentiation/ExponentiationBinary.hpp"
#include "../../include/ANTL/Exponentiation/ExponentiationNAF.hpp"
#include "../../include/ANTL/Exponentiation/ExponentiationL2R.hpp"
#include "../../include/ANTL/Exponentiation/ExponentiationWNAF.hpp"
//#include "../../include/ANTL/Exponentiation/ExponentiationSB3.hpp"
#include "../../include/ANTL/Exponentiation/ExponentiationYao.hpp"
//#include "../../include/ANTL/Exponentiation/ExponentiationDoubleBaseStrictChain.hpp"
//#include "../../include/ANTL/Exponentiation/ExponentiationDoubleBaseUnsignedGreedy.hpp"
//#include "../../include/ANTL/Exponentiation/ExponentiationDoubleBaseSignedGreedy.hpp"
//#include "../../include/ANTL/Exponentiation/ExponentiationExtendedDoubleBase.hpp"



NTL_CLIENT
using namespace ANTL;

	
int main (int argc, char **argv)
{
  zz_p a, b, c, b_bin,b_naf, b_l2r, b_wnaf, b_sb3, b_yao, b_dbsc, b_dbug, b_dbsg, b_exdb;
  ZZ n,m;

  // use GF(1073741827) for these tests
  zz_p::init(1073741827);

  // set base to be a random value in GF(1073741827)
  do {
    random(a);
  } while (IsOne(a));

    b = a;
    m = n;

  // generate random exponent of size 512 bits
  RandomLen (n, 512);

  cout << "Using:" << endl;
  cout << " p = " << zz_p::modulus() << endl;
  cout << " a = " << a << endl;
  cout << " n = " << n << endl;


  // initialize exponentiation classes
  ExponentiationBinary<zz_p> ebin;
  ExponentiationNAF<zz_p> enaf;
  
  // compute a^n with available methods
  ebin.power(b_bin,a,n);

  enaf.initialize(a,n);
  enaf.power(b_naf,a,n);
    
    enaf.initialize(b,m);
    enaf.power(c,b,m);

  

  // check and output results
  cout << "a^n (binary) = " << b_bin << endl;
  cout << "a^n (naf)    = " << b_naf << endl;
    cout << "b^m (naf)    = " << c << endl;
  
  if ((b_bin == b_naf))
    cout << "RESULTS MATCH!" << endl;
  else
  {
    cout << "ERROR:  RESULTS DO NOT MATCH!" << endl;
    exit(1);
  }
}
