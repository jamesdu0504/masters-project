/**
 * @file ExponentiationNAF_impl.hpp
 * @author Michael Jacobson
 * @brief generic implementation of templated methods from ExponentiationNAF class
 */

using namespace ANTL;

//
// initialize class for computing A^n using left-to-right NAF
//
template < class T >
void ExponentiationNAF<T>::initialize (const T &A, const ZZ & n)
{
    // compute A^-1
    inv(Ainv,A);
    
    // initialize digit vector to size NumBits(n)+1
    e.reserve(NumBits(n)+1);
}


//
// compute A^n using standard left-to-right NAF method
//
template <class T>
void ExponentiationNAF<T>::power (T &C, const T &A, const ZZ & n)
{
    //--- Added during Parth's project
    e.clear();
    //cout << "Ainv = " << Ainv;
    //---
    
    // compute NAF expansion of n (right-to-left)
    ZZ ex = abs (n);
    while (ex > 0)
    {
        if (IsOdd (ex))
        {
            short ei = 2 - (short) rem(ex,4);
            e.push_back(ei);
            sub(ex,ex,ei);
        }
        else
            e.push_back(0);
        RightShift(ex,ex,1);
    }
    
    // compute C = A^n from left-to-right using NAF digits in e
    assign(C, A);
    //cout << "After assign = " << C;
    
    for (long j = e.size()-2; j >= 0; --j)
    {
        sqr(C, C);
        if (e.at(j) == 1)
        {
            mul(C, C, A);
        }
            
        else if (e.at(j) == -1)
        {
            mul(C, C, Ainv);
        }
    }
}
