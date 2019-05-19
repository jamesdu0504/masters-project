/**
 * @file ExponentiationWNAF_impl.hpp
 * @author Michael Jacobson/Reginald Lybbert
 * @brief generic implementation of templated methods from ExponentiationWNAF class
 */

using namespace ANTL;

//
// initialize class for computing A^n using left-to-right WNAF
//
template < class T >
void ExponentiationWNAF<T>::initialize(const T &A, const ZZ & n, const short width)
{
    w = width;
    
    //create precomputed array, precomp[i] corresponds to A^{2i+1}
    int len = 1 << (w - 2);
    
    precomp = new T[len];
    precompInv = new T[len];
    
    T square;
    sqr(square, A);
    assign(precomp[0], A);
    inv(precompInv[0],precomp[0]);
    
    for(int i = 1; i < len; i++)
    {
        mul(precomp[i],precomp[i-1],square);
        inv(precompInv[i],precomp[i]);
    }
    
    // initialize digit vector to size NumBits(n)+1
    e.reserve(NumBits(n)+1);
}


//
// compute A^n using standard left-to-right NAF method
//
template < class T > 
void ExponentiationWNAF<T>::power (T &C, const T &A, const ZZ & n)
{
    //--- Added during Parth's project
    e.clear();
    //---
    
    // compute NAF expansion of n (right-to-left)
    ZZ ex = abs (n);
    while (ex > 0)
    {
        if (IsOdd (ex))
        {
            short ei = (short) rem(ex,(1 << w));
            if (ei > (1 << (w-1)))
                ei = ei - (1 << w);
            e.push_back(ei);
            sub(ex,ex,ei);
        }
        
        else
            e.push_back(0);
        RightShift(ex,ex,1);
    }
    
    int x;
    
    // compute C = A^n from left-to-right using NAF digits in e
    x = (e.at(e.size() - 1) - 1)/2;
    assign(C,precomp[x]);
    for ( long j = e.size()-2; j >= 0; --j)
    {
        sqr(C, C);
        if (e.at(j) > 0)
        {
            x = (e.at(j) - 1)/2;
            mul(C, C, precomp[x]);
        }
        
        else if (e.at(j) < 0)
        {
            x = ((-1)*e.at(j) - 1)/2;
            mul(C,C,precompInv[x]);
        }
    }
}
