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
#include "../../include/ANTL/Exponentiation/ExponentiationSB3.hpp"
#include "../../include/ANTL/Exponentiation/ExponentiationYao.hpp"
//#include "../../include/ANTL/Exponentiation/ExponentiationDoubleBaseStrictChain.hpp"
//#include "../../include/ANTL/Exponentiation/ExponentiationDoubleBaseUnsignedGreedy.hpp"
//#include "../../include/ANTL/Exponentiation/ExponentiationDoubleBaseSignedGreedy.hpp"
//#include "../../include/ANTL/Exponentiation/ExponentiationExtendedDoubleBase.hpp"



NTL_CLIENT
using namespace ANTL;

	
int main (int argc, char **argv)
{
  zz_p a,b_bin,b_naf, b_l2r, b_wnaf, b_sb3, b_yao, b_dbsc, b_dbug, b_dbsg, b_exdb;
  ZZ n;

  // use GF(1073741827) for these tests
  zz_p::init(1073741827);

  // set base to be a random value in GF(1073741827)
  do {
    random(a);
  } while (IsOne(a));


  // generate random exponent of size 512 bits
  RandomLen (n, 512);

  cout << "Using:" << endl;
  cout << " p = " << zz_p::modulus() << endl;
  cout << " a = " << a << endl;
  cout << " n = " << n << endl;


  // initialize exponentiation classes
  ExponentiationBinary<zz_p> ebin;
  ExponentiationNAF<zz_p> enaf;
  ExponentiationL2R<zz_p> el2r;
  ExponentiationWNAF<zz_p> ewnaf;
  ExponentiationSB3<zz_p> esb3;
  ExponentiationYao<zz_p> eyao;
  ExponentiationDoubleBaseStrictChain<zz_p> edbsc;
  ExponentiationDoubleBaseUnsignedGreedy<zz_p> edbug;
  ExponentiationDoubleBaseSignedGreedy<zz_p> edbsg; 
  ExponentiationExtendedDoubleBase<zz_p> eexdb;

  // compute a^n with available methods
  ebin.power(b_bin,a,n);

  enaf.initialize(a,n);
  enaf.power(b_naf,a,n);

  el2r.initialize(a);
  el2r.power(b_l2r,a,n);

  ewnaf.initialize(a,n,5);
  ewnaf.power(b_wnaf,a,n);

  esb3.initialize(a,n);
  esb3.power(b_sb3,a,n);

  eyao.power(b_yao,a,n);

  edbsc.initialize(a);
  edbsc.power(b_dbsc,a,n);

  edbug.initialize(500,300);
  edbug.power(b_dbug,a,n);

  edbsg.setBounds(500,300);
  edbsg.power(b_dbsg,a,n);  

  eexdb.setBounds(500,300);
  eexdb.power(b_exdb,a,n);

  // check and output results
  cout << "a^n (binary) = " << b_bin << endl;
  cout << "a^n (naf)    = " << b_naf << endl;
  cout << "a^n (l2r)    = " << b_l2r << endl;
  cout << "a^n (wnaf)   = " << b_wnaf << endl;
  cout << "a^n (sb3)    = " << b_sb3 << endl;
  cout << "a^n (byao)   = " << b_yao << endl;
  cout << "a^n (dbsc)   = " << b_dbsc << endl; 
  cout << "a^n (dbug)   = " << b_dbug << endl;
  cout << "a^n (dbsg)   = " << b_dbsg << endl;
  cout << "a^n (exdb)   = " << b_exdb << endl;
 
  if ((b_bin == b_naf) && (b_naf == b_l2r) && (b_wnaf == b_bin) && (b_bin == b_sb3) && (b_bin == b_yao) && (b_dbsc == b_bin) && (b_dbug == b_bin) && (b_bin == b_dbsg) && (b_bin == b_exdb))
    cout << "RESULTS MATCH!" << endl;
  else
  {
    cout << "ERROR:  RESULTS DO NOT MATCH!" << endl;
    exit(1);
  }
}
