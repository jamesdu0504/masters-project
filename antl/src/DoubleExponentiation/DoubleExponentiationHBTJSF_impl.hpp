/**
 * @file DoubleExponentiationHBTJSF_impl.hpp
 * @author Michael Jacobson/Reginald Lybbert
 * @brief generic implementation of templated methods from DoubleExponentiationHBTJSF class
 */

using namespace ANTL;


//look up appropriate precomputed value (inverting if necessary)
template < class T >
void DoubleExponentiationHBTJSF<T>::readPrecomps(T &C, short mDigit, short nDigit)
{
    int switcher = 6*mDigit + nDigit;
    switch(switcher)
    {
        case -13:  //(-2,-1)
            assign(C,precomps[4]);
            inv(C,C);
            break;
        case -11:  //(-2,1)
           assign(C,precomps[5]);
           inv(C,C);
           break;
        case -9:   // (-2,3)
           assign(C,precomps[15]);
           inv(C,C);
           break;
        case -8:   // (-1,-2)
           assign(C,precomps[8]);
           inv(C,C);
           break;
        case -7:   // (-1,-1)
           assign(C,precomps[2]);
           inv(C,C);
           break;
        case -6:   // (-1,0)
           assign(C,precomps[0]);
           inv(C,C);
           break;
        case -5:   // (-1,1)
           assign(C,precomps[3]);
           inv(C,C);
           break;
        case -4:   // (-1,2)
           assign(C,precomps[9]);
           inv(C,C);
           break;
        case -3:   // (-1,3)
           assign(C,precomps[13]);
           inv(C,C);
           break;
        case -1:   // (0,-1)
           assign(C,precomps[1]);
           inv(C,C);
           break;
        case 1:   // (0,1)
           assign(C,precomps[1]);
           break;
        case 4:   // (1,-2)
           assign(C,precomps[9]);
           break;
        case 5:   // (1,-1)
           assign(C,precomps[3]);
           break;
        case 6:   // (1,0)
           assign(C,precomps[0]);
           break;
        case 7:   // (1,1)
           assign(C,precomps[2]);
           break;
        case 8:   // (1,2)
           assign(C,precomps[8]);
           break;
        case 9:   // (1,3)
           assign(C,precomps[12]);
           break;
        case 11:   // (2,-1)
           assign(C,precomps[5]);
           break;
        case 13:   // (2,1)
           assign(C,precomps[4]);
           break;
        case 15:   // (2,3)
           assign(C,precomps[14]);
           break;
        case 16:   // (3,-2)
           assign(C,precomps[11]);
           break;
        case 17:   // (3,-1)
           assign(C,precomps[7]);
           break;
        case 19:   // (3,1)
           assign(C,precomps[6]);
           break;
        case 20:   // (3,2)
           assign(C,precomps[10]);
           break;
        default:   // (3,2)
            assign(C,precomps[0]);
            cout << "Error in parsing representation; got digits: " << mDigit << "," << nDigit << endl;
            cout << "All digits should be in {-2,-1,0,1,2,3}, and should be coprime." << endl;
    }
}

//
// initialize class for computing A^mB^n using hybrid binary-ternary joint sparse form
//
template < class T >
void DoubleExponentiationHBTJSF<T>::initializeBase(const T &A, const T&B)
{
    //create precomputed values
    T Binv;
    
    inv(Binv,B);
    assign(precomps[0],A);
    assign(precomps[1],B);
    mul(precomps[2],precomps[0],precomps[1]);     //AB
    mul(precomps[3],precomps[0],Binv);            //AB^-1
    mul(precomps[4],precomps[2],precomps[0]);     //A^2B
    mul(precomps[5],precomps[3],precomps[0]);     //A^2B^-1
    mul(precomps[6],precomps[4],precomps[0]);     //A^3B
    mul(precomps[7],precomps[5],precomps[0]);     //A^3B^-1
    mul(precomps[8],precomps[2],precomps[1]);     //AB^2
    mul(precomps[9],precomps[3],Binv);            //AB^-2
    mul(precomps[10],precomps[6],precomps[1]);    //A^3B^2
    mul(precomps[11],precomps[7],Binv);           //A^3B^-2
    mul(precomps[12],precomps[8],precomps[1]);    //AB^3
    mul(precomps[13],precomps[9],Binv);           //AB^-3
    mul(precomps[14],precomps[12],precomps[0]);   //A^2B^3
    mul(precomps[15],precomps[13],precomps[0]);   //A^2B^-3
} 


//Produces the HBTJSF representation of m and n, note that the most significant digit is at the back of the vectors.
template < class T >
void DoubleExponentiationHBTJSF<T>::initializeExponent(const ZZ &m, const ZZ &n)
{
    initM = m;
    initN = n;
    
    ZZ k1 = m;
    ZZ k2 = n;
    
    while(k1 > 0 || k2 > 0)
    {
        if(k1 % 3 == 0 && k2 % 3 == 0)
        {
            base.push_back(3);
            em.push_back(0);
            en.push_back(0);
            k1 = k1/3;
            k2 = k2/3;
        }
        
        else if(k1 % 2 == 0 && k2 % 2 == 0)
        {
            base.push_back(2);
            em.push_back(0);
            en.push_back(0);
            
            k1 = k1/2;
            k2 = k2/2;
        }
        
        else
        {
            base.push_back(2);
            em.push_back(k1%6>3?(k1%6-6):k1%6);
            en.push_back(k2%6>3?(k2%6-6):k2%6);
            k1 = (k1 - em.back())/2;
            k2 = (k2 - en.back())/2;
        }
    }
}


//
// compute A^mB^n using hybrid binary-ternary joint sparse form (HBTJSF) method
//
template < class T > 
void DoubleExponentiationHBTJSF<T>::power (T &C, const T &A, const T &B, const ZZ & m, const ZZ & n)
{
    // Clear before starting - added during Parth's project
    initM = 0;
    initN = 0;
    
    //Check if bases have been initialized
    if(!equals(precomps[0],A) || !equals(precomps[1],B))
    {
        initializeBase(A,B);
    }
    
    //Check if exponents have been initialized
    if(initM != m || initN != n)
    {
        initializeExponent(m,n);
    }
    
    //Do the exponentiation
    T temp;
    readPrecomps(C, em.back(), en.back());
    base.pop_back();
    em.pop_back();
    en.pop_back();
    while(!base.empty())
    {
        if(base.back() == 2)
        {
            sqr(C,C);
        }
        
        else if(base.back() == 3)
        {
            cub(C,C);
        }
        
        if(em.back() != 0 || en.back() != 0)
        {
            readPrecomps(temp, em.back(), en.back());
            mul(C,C,temp);
        }
        
        base.pop_back();
        em.pop_back();
        en.pop_back();
    }
}	
