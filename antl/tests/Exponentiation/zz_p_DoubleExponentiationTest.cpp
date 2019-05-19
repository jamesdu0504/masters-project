/**
 * @file zz_p_DoubleExponentiationTest.cpp
 * @author Michael Jacobson/Reginald Lybbert
 * @brief Test program for double exponentiation routines using zz_p as base type.
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

#include <ANTL/Exponentiation/ExponentiationBinary.hpp>
#include <ANTL/DoubleExponentiation/DoubleExponentiationIL.hpp>
#include <ANTL/DoubleExponentiation/DoubleExponentiationJSF.hpp>
#include <ANTL/DoubleExponentiation/DoubleExponentiationHBTJSF.hpp>

NTL_CLIENT
using namespace ANTL;
	
int main (int argc, char **argv)
{
  zz_p a,b,ca_bin,cb_bin,c_bin, c_il, c_jsf, c_hbtjsf;
  ZZ m,n;

  // use GF(1073741827) for these tests
  zz_p::init(1073741827);

  // set base to be a random value in GF(5)
  do {
    random(a);
  } while (IsOne(a));

  do {
    random(b);
  } while (IsOne(b));

  // generate random exponent of size 128 bits
  RandomLen (m, 512);
  RandomLen (n, 512);

  cout << "Using:" << endl;
  cout << " p = " << zz_p::modulus() << endl;
  cout << " a = " << a << endl;
  cout << " b = " << b << endl;
  cout << " m = " << m << endl;
  cout << " n = " << n << endl;

  // initialize exponentiation classes
  ExponentiationBinary<zz_p> ebin;
  DoubleExponentiationIL<zz_p> deil;
  DoubleExponentiationJSF<zz_p> djsf;
  DoubleExponentiationHBTJSF<zz_p> dhbt;

  // compute a^mb^n with available methods
  ebin.power(ca_bin,a,m);
  ebin.power(cb_bin,b,n);
  mul(c_bin,ca_bin,cb_bin);

  deil.initialize(a,b,m,n,5,5);
  deil.power(c_il,a,b,m,n);

  djsf.initialize(a,b,m,n);
  djsf.power(c_jsf,a,b,m,n);

  dhbt.power(c_hbtjsf,a,b,m,n);

  // check and output results
  cout << "a^mb^n (naive) = " << c_bin << endl;
  cout << "a^mb^n (interleaving)    = " << c_il << endl;
  cout << "a^mb^n (joint sparse form) = " << c_jsf << endl;
  cout << "a^mb^n (hybrid binary-ternary joint sparse form) = " << c_hbtjsf << endl;
  
  if ((c_bin == c_il) and (c_il == c_jsf) and (c_jsf == c_hbtjsf))
    cout << "RESULTS MATCH!" << endl;
  else
    cout << "ERROR:  RESULTS DO NOT MATCH!" << endl;
}
