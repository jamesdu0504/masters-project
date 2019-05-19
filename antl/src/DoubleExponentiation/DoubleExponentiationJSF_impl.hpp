/**
 * @file DoubleExponentiationJSF_impl.hpp
 * @author Michael Jacobson/Reginald Lybbert
 * @brief generic implementation of templated methods from DoubleExponentiationJSF class
 */

using namespace ANTL;

//
// initialize class for computing A^mB^n using left-to-right joint sparse form
//
template < class T >
void DoubleExponentiationJSF<T>::initialize(const T &A, const T&B, const ZZ &m, const ZZ & n)
{
    //create precomputed values
    inv(AInv,A);
    inv(BInv,B);
    mul(AB, A, B);
    mul(ABInv, A, BInv);
    inv(AInvB, ABInv);
    inv(AInvBInv, AB);
    
    // initialize digit vector to size NumBits(max(m,n))+1
    em.reserve(NumBits(max(m,n))+1);
    en.reserve(NumBits(max(m,n))+1);
}


//
// compute A^mB^n using JSF method
//
template < class T > 
void DoubleExponentiationJSF<T>::power (T &C, const T &A, const T &B, const ZZ & m, const ZZ & n)
{
    // Clear before starting - added during Parth's project
    em.clear();
    en.clear();
    
    // compute JSF expansion of m and n (right-to-left)
    ZZ mt = m;
    ZZ nt = n;
    
    short d1 = 0;
    short d2 = 0;
    ZZ l1, l2;
    int u;

    while((mt + d1 > 0) or (nt + d2 > 0))
    {
        add(l1, mt, d1);
        add(l2, nt, d2);
        if(IsOdd(l1))
        {
            u =  l1 % 4;
            if(u == 3)
                u = -1;
            if(((l1 % 8 == 3) or (l1 % 8 == 5)) and (l2 % 4 == 2))
                u = -u;
        }
        else
            u = 0;
        em.push_back(u);
        if(2 * d1 == 1 + u)
            d1 = 1 - d1;
        RightShift(mt, mt, 1);
        if(IsOdd(l2))
        {
            u =  l2 % 4;
            if(u == 3)
                u = -1;
            if(((l2 % 8 == 3) or (l2 % 8 == 5)) and (l1 % 4 == 2))
                u = -u;
        }
        else
            u = 0;
        en.push_back(u);
        if(2 * d2 == 1 + u)
            d2 = 1 - d2;
        RightShift(nt,nt,1);
    }
    
    // compute C = A^mB^n from left-to-right using JSF digits in em and en
    if(em.at(em.size() - 1) == 1)
    {
        if(en.at(em.size() - 1) == 1)
            assign(C,AB);
        else
            assign(C,A);
    }
    else
        assign(C,B);
    for (long j = em.size()-2; j >= 0; --j)
    {
        sqr(C, C);
        if(em.at(j) == 1)
        {
            if(en.at(j) == -1)
                mul(C, C, ABInv);
            else if(en.at(j) == 0)
                mul(C, C, A);
            else
                mul(C, C, AB);
        }
        
        else if(em.at(j) == 0)
        {
            if(en.at(j) == -1)
                mul(C, C, BInv);
            else if(en.at(j) == 1)
                mul(C,C,B);
        }
        
        else if(em.at(j) == -1)
        {
            if(en.at(j) == -1)
                mul(C,C,AInvBInv);
            else if(en.at(j) == 0)
                mul(C,C,AInv);
            else if(en.at(j) == 1)
                mul(C,C,AInvB);
        }
    }
}
