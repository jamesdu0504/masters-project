/**
 * @file ExponentiationBinary_impl.hpp
 * @author Michael Jacobson
 * @brief generic implementation of templated method from ExponentiationBinary class
 */

using namespace ANTL;

//
// compute A^n using standard left-to-right binary method
//
template < class T >
void ExponentiationBinary<T>::power (T &C, const T &A, const ZZ & n)
{
  assign(C,A);
  for (long i = NumBits(n)-2 ; i >= 0 ; i--)
    {
      sqr(C, C);
      if (bit(n, i) == 1)
          mul(C, C, A);
    }
}
