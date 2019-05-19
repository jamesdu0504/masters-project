/**
 * @file DoubleExponentiationIL_impl.hpp
 * @author Michael Jacobson/Reginald Lybbert
 * @brief generic implementation of templated methods from DoubleExponentiationIL class
 */

using namespace ANTL;

//
// initialize class for computing A^mB^n using left-to-right interleaving
//
template < class T >
void DoubleExponentiationIL<T>::initialize(const T &A, const T&B, const ZZ &m, const ZZ & n, const short w0, const short w1)
{
    wm = w0;
    wn = w1;
    
    //create precomputed array, precomp[i] corresponds to A^{2i+1}
    int lenm = 1 << (wm - 2);
    int lenn = 1 << (wn - 2);
    
    precompA.resize(lenm);
    precompAInv.resize(lenm);
    precompB.resize(lenn);
    precompBInv.resize(lenn);

    T squareA;
    sqr(squareA, A);
    assign(precompA[0], A);
    inv(precompAInv[0],precompA[0]);

    for(int i = 1; i < lenm; i++)
    {
        mul(precompA[i],precompA[i-1],squareA);
        inv(precompAInv[i],precompA[i]);
    }
    
    T squareB;
    sqr(squareB, B);
    assign(precompB[0], B);
    inv(precompBInv[0],precompB[0]);

    for(int i = 1; i < lenn; i++)
    {
        mul(precompB[i],precompB[i-1],squareB);
        inv(precompBInv[i],precompB[i]);
    }
    
    // initialize digit vector to size NumBits(n)+1
    em.reserve(NumBits(max(m,n))+1);
    en.reserve(NumBits(max(m,n))+1);
}


//
// compute A^mB^n using standard left-to-right NAF method
//
template < class T > 
void DoubleExponentiationIL<T>::power (T& C, const T& A, const T& B, const ZZ& m, const ZZ& n)
{
    // Clear before starting - added during Parth's project
    em.clear();
    en.clear();
    
    // compute WNAF expansion of m and n (right-to-left)
    ZZ ex = abs (m);
    while (ex > 0) 
    {
        if (IsOdd (ex)) 
        {
            short ei = (short) rem(ex,(1 << wm));
            if (ei > (1 << (wm-1)))
                ei = ei - (1 << wm);
            em.push_back(ei);
            sub(ex,ex,ei);
        }
        
        else
            em.push_back(0);
        RightShift(ex,ex,1);
    }
    
    ex = abs (n);
    while (ex > 0)
    {
        if (IsOdd (ex))
        {
            short ei = (short) rem(ex,(1 << wn));
            if (ei > (1 << (wn-1)))
                ei = ei - (1 << wn);
            en.push_back(ei);
            sub(ex,ex,ei);
        }
        
        else
            en.push_back(0);
        RightShift(ex,ex,1);
    }
    
    while(em.size() < en.size())
    {
        em.push_back(0);
    }
    
    while(en.size() < em.size())
    {
        en.push_back(0);
    }
    
    int x;
    // compute C = A^mB^n from left-to-right using IL digits in em and en
    if (em.at(em.size() - 1) != 0)
    {
        x = (em.at(em.size() - 1) - 1)/2;
        assign(C,precompA[x]);
        if(en.at(en.size() - 1) != 0)
        {
            x = (en.at(en.size() - 1) - 1)/2;
            mul(C,C,precompB[x]);
        }
    }
    
    else
    {
        x = (en.at(en.size() - 1) - 1)/2;
        assign(C, precompB[x]);
    }
    
    for (long j = em.size()-2; j >= 0; --j)
    {
        sqr(C, C);
        if (em.at(j) > 0)
        {
            x = (em.at(j) - 1)/2;
            mul(C, C, precompA[x]);
        }
        
        else if (em.at(j) < 0)
        {
            x = ((-1)*em.at(j) - 1)/2;
            mul(C,C,precompAInv[x]);
        }
        
        if (en.at(j) > 0)
        {
            x = (en.at(j) - 1)/2;
            mul(C, C, precompB[x]);
        }
        
        else if (en.at(j) < 0)
        {
            x = ((-1)*en.at(j) - 1)/2;
            mul(C,C,precompBInv[x]);
        }
    }
}
