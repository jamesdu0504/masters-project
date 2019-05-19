/**
 * @file ExponentiationSB3_impl.hpp
 * @author Michael Jacobson/Reginald Lybbert
 * @brief generic implementation of templated methods from ExponentiationSB3 class
 */

using namespace ANTL;

//
// initialize class for computing A^n using left-to-right SB3
//
template < class T >
void 
ExponentiationSB3<T>::initialize(const T &A, const ZZ & n)
{
  // compute A^-1
  inv(Ainv,A);

  // initialize digit vector to size (2*NumBits(n))/3+2
  e.reserve((2*NumBits(n))/3+2);
}


//
// compute A^n using standard left-to-right SB3 method
//
template < class T > 
void 
ExponentiationSB3<T>::power (T &C, const T &A, const ZZ & n)
{
  // compute SB3 expansion of n (right-to-left)
  ZZ ex = abs (n);
  while (ex > 0)
    {
      short ei = (short) rem(ex,3);
      if (ei == 2) {
	ei = -1;
      }
      e.push_back(ei);
      sub(ex,ex,ei);
      //cout << ei;
     div(ex,ex,3);
    }
  //cout << endl;
  // compute C = A^n from left-to-right using SB3 digits in e
  assign(C,A);
  for (register long j = e.size()-2; j >= 0; --j)
    {
      cub(C, C);
      if (e.at(j) == 1)
	mul(C, C, A);
      else if (e.at(j) == -1)
	mul(C,C,Ainv);
    }
}
