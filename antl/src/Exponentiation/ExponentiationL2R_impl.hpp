/**
 * @file ExponentiationL2R_impl.hpp
 * @author Michael Jacobson/Reginald Lybbert
 * @brief generic implementation of templated methods from ExponentiationL2R class
 */

using namespace ANTL;

//
// initialize class for computing A^n using left-to-right NAF
//
template < class T >
void ExponentiationL2R<T>::initialize(const T &A)
{
    // compute A^-1
    inv(Ainv,A);
}


//
// compute A^n using left-to-right NAF method without storing the full NAF representation
//
template < class T > 
void ExponentiationL2R<T>::power (T &C, const T &A, const ZZ & n)
{
    //--- Added during Parth's project
    //---
    
    long m, i, j;
    short b;
  
    //Initialize lookup tables
    short Btable [16] = {0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 1, 1};
    short Rtable [16] = {0, 0, 0, 1, 1, 1, 0, 0, 0, 0, -1, -1, -1, 0, 0, 0};

    assign(C,A);

    m = NumBits(n) - 2;
  
    //Initialize after first digit of NAF representation
    if (bit(n,m) == 0)
        b = 0;
    else
    {
        b = 1;
        m = m + 1;
    }
    
    //Find index for next digit in the lookup table
    for(j = m; j >= 0; j--)
    {
        i = (b << 3) + (bit(n,j) << 2) + (bit(n,j-1) << 1) + bit(n,j-2);
        //if j-1 or j-2 is less than zero, bit returns 0
        
        b = Btable[i];
        
        //perform appropriate operations on the base
        sqr(C,C);
        if (Rtable[i] == 1)
            mul(C, C, A);
        else if(Rtable[i] == -1)
            mul(C,C,Ainv);
    }
}
