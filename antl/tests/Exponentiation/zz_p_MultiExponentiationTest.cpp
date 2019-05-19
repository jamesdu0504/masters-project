/**
 * @file zz_p_MultiExponentiationTest.cpp
 * @author Michael Jacobson/Reginald Lybbert
 * @brief Test program for multiexponentiation routines using zz_p as base type.
 */

#include <NTL/lzz_p.h>

inline void id(NTL::zz_p& A){
    A = 1;
}

#include <ANTL/Exponentiation/ExponentiationBinary.hpp>
#include <ANTL/MultiExponentiation/MultiExponentiationPippenger.hpp>

NTL_CLIENT
using namespace ANTL;
	
int main (int argc, char **argv)
{
  zz_p a,b,c,ca1_bin,cb1_bin,cc1_bin,ca2_bin,cb2_bin,cc2_bin,ca3_bin,cb3_bin,cc3_bin,c1_bin,c2_bin,c3_bin;
  ZZ m1,n1,r1,m2,n2,r2,m3,n3,r3;
  vector<zz_p> A,C;
  vector<vector<ZZ> > n;

  // use GF(1073741827) for these tests
  zz_p::init(1073741827);

  // set base to be a random value in GF(5)
  do {
    random(a);
  } while (IsOne(a));

  do {
    random(b);
  } while (IsOne(b));

  do {
    random(c);
  } while (IsOne(c));

  // generate random exponents of size 128 bits
  RandomLen (m1, 512);
  RandomLen (n1, 512);
  RandomLen (r1, 512);
  RandomLen (m2, 512);
  RandomLen (n2, 512);
  RandomLen (r2, 512);
  RandomLen (m3, 512);
  RandomLen (n3, 512);
  RandomLen (r3, 512);

  cout << "Using:" << endl;
  cout << " p = " << zz_p::modulus() << endl;
  cout << " a = " << a << endl;
  cout << " b = " << b << endl;
  cout << " c = " << c << endl;
  cout << " m1 = " << m1 << endl;
  cout << " n1 = " << n1 << endl;
  cout << " r1 = " << r1 << endl;
  cout << " m2 = " << m2 << endl;
  cout << " n2 = " << n2 << endl;
  cout << " r2 = " << r2 << endl;
  cout << " m3 = " << m3 << endl;
  cout << " n3 = " << n3 << endl;
  cout << " r3 = " << r3 << endl;

  // initialize exponentiation classes
  ExponentiationBinary<zz_p> ebin;
  MultiExponentiationPippenger<zz_p> pipp;

  // compute a^{mi}b^{ni}c^{ri} with available methods
  ebin.power(ca1_bin,a,m1);
  ebin.power(cb1_bin,b,n1);
  ebin.power(cc1_bin,c,r1);
  mul(c1_bin,ca1_bin,cb1_bin);
  mul(c1_bin,c1_bin,cc1_bin);

  ebin.power(ca2_bin,a,m2);
  ebin.power(cb2_bin,b,n2);
  ebin.power(cc2_bin,c,r2);
  mul(c2_bin,ca2_bin,cb2_bin);
  mul(c2_bin,c2_bin,cc2_bin);

  ebin.power(ca3_bin,a,m3);
  ebin.power(cb3_bin,b,n3);
  ebin.power(cc3_bin,c,r3);
  mul(c3_bin,ca3_bin,cb3_bin);
  mul(c3_bin,c3_bin,cc3_bin);


  A.push_back(a);
  A.push_back(b);
  A.push_back(c);

  n.resize(3);
  n[0].push_back(m1);
  n[0].push_back(n1);
  n[0].push_back(r1); 
  n[1].push_back(m2);
  n[1].push_back(n2);
  n[1].push_back(r2);
  n[2].push_back(m3);
  n[2].push_back(n3);
  n[2].push_back(r3); 

  pipp.power(C,A,n);

  // check and output results
  cout << endl;
  cout << "a^{m1}b^{n1}c^{r1} (binary) = " << c1_bin << endl;
  cout << "a^{m2}b^{n2}c^{r2} (binary) = " << c2_bin << endl;
  cout << "a^{m3}b^{n3}c^{r3} (binary) = " << c3_bin << endl << endl;

  cout << "a^{m1}b^{n1}c^{r1} (pippenger) = " << C[0] << endl;
  cout << "a^{m2}b^{n2}c^{r2} (pippenger) = " << C[1] << endl;
  cout << "a^{m3}b^{n3}c^{r3} (pippenger) = " << C[2] << endl << endl;

  if ((c1_bin == C[0]) and (c2_bin == C[1]) and (c3_bin == C[2]))
    cout << "RESULTS MATCH!" << endl;
  else
    cout << "ERROR:  RESULTS DO NOT MATCH!" << endl;

}


