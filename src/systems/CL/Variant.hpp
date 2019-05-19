#ifndef Variant_HPP
#define Variant_HPP

#include "CL.hpp"


class Variant : public CL
{
    public:
    
    Variant() {};
    Variant(const int WIN, const int MS, const int DS, const int ES) : CL (WIN, MS, DS, ES) {};
    ~Variant() {};
        
    void          keygen    (CLPublickey& clpk,  CLSecretkey& clsk);
    
    // CL
    CLCiphertext& encrypt   (CLPlaintext& clpt,  CLPublickey& clpk);
    CLPlaintext&  decrypt   (CLCiphertext& clct, CLPublickey& clpk, CLSecretkey& clsk);
    CLPlaintext&  mdecrypt   (CLCiphertext& clct, CLPublickey& clpk, CLSecretkey& clsk);
    
    // DJS
    CLCiphertext& djsencrypt   (CLPlaintext& clpt,  CLPublickey& clpk);
    CLPlaintext&  djsdecrypt   (CLCiphertext& clct, CLPublickey& clpk, CLSecretkey& clsk);
    
    // CRT versions of CL and DJS by DJS
    CLPlaintext&  decryptCRT  (CLCiphertext& clct, CLPublickey& clpk, CLSecretkey& clsk);
    CLPlaintext&  mdecryptCRT (CLCiphertext& clct, CLPublickey& clpk, CLSecretkey& clsk);
    CLCiphertext& encrypt2    (CLPlaintext& clpt,  CLPublickey& clpk);
    CLPlaintext&  decrypt2    (CLCiphertext& clct, CLPublickey& clpk, CLSecretkey& clsk);
    CLCiphertext& encrypt3    (CLPlaintext& clpt,  CLPublickey& clpk);
    CLPlaintext&  decrypt3    (CLCiphertext& clct, CLPublickey& clpk, CLSecretkey& clsk);

    
    // Homomorphic functions
    CLCiphertext& evalsum   (CLCiphertext& clct1, CLCiphertext& clct2);
    CLCiphertext& evalscal  (CLCiphertext& clct, const mpz_t alpha);
    CLCiphertext& evalsum3  (CLCiphertext& clct1, CLCiphertext& clct2);
    CLCiphertext& evalscal3 (CLCiphertext& clct, const mpz_t alpha);

}; // end of class

inline void Variant::keygen(CLPublickey& clpk,  CLSecretkey& clsk)
{
    // Set n and t parameters of public key
    clpk.set(n,t);
    
    // Get upper bound on odd part; sets CL variable B
    getB(DK, mu);
    
    // Set upper bound on size of G
    mpz_mul(UB, B, con);
    ZZ_limbs_set(zUB, UB[0]._mp_d, UB[0]._mp_size);
    
    // Set DK to -DK
    mpz_neg(DK, DK);
    
    // Compute Dcon
    mpz_mul(Dcon, con, con);
    mpz_mul(Dcon, Dcon, DK);
    
    // Compute Dconsqr
    mpz_mul(Dconsqr, con, con);
    mpz_mul(Dconsqr, Dconsqr, Dcon);
        
    // Set F to Castagnos ideal (con^2, con) of discriminant Dcon
    mpz_mul(temp1, con, con);
    F.assign(temp1, con, Dcon);
        
    // Compute ideal R
    mpz_set_ui(r, 2);
    mpz_gcd(temp1, r, con);
    while ( !( (mpz_legendre(DK, r) == 1) && (mpz_cmp_ui(temp1, 1) == 0) ) )
    {
        mpz_nextprime(r, r);
        mpz_gcd(temp1, r, con);
    }
    prime_r = mpz_get_ui(r);
    prime_index = prime_index_ge(prime_r);
    prime_index = R.next_prime_QF(prime_index, DK);
        
    // Set G = R^2 of discriminant DK
    sqr(G, R);
        
    // Generate secret key x in [0, Bcon-1] = [0,UB)
    if (ebits == 0)
        zx = RandomBnd(zUB);
    else
        zx = RandomLen_ZZ(ebits);
    
    // Compute H = G^x
    qfeG.initialize(G, zx);
    qfeG.power(H, G, zx);
    
    // Generate a dummy exponent for initialisation
    e1 = RandomBnd(zUB);

    // Set Fi's for N > 1
    for (int i = 0; i < n; i++)
    {
        mpz_mul(temp1, pntn[i], pntn[i]);
        Fn[i].assign(temp1, pntn[i], Dcon);
        
        qfeFn[i].initialize(Fn[i], e1);
    }

    // Set exponentiation objects for F and H bases for use in encrypt
    qfeF.initialize(F, e1); // F -- uses qfeF
    H.setD(DK);
    qfeH.initialize(H, e1); // H -- uses qfeH
    
    // Set
    clpk.set(UB, con, DK, Dcon, F, G, H);
    clsk.set(zx);
}



// ------------------------------------------------------------------------------
// ------------------------------------------------------------------------------
// ------------------------------------------------------------------------------
// ---------------------- Encryption and Decryption by CL -----------------------
// ------------------------------------------------------------------------------
// ------------------------------------------------------------------------------
// ------------------------------------------------------------------------------

/**
 *
 * Should be used when
 *  1. fbits <= fmaxbits
 *  2. fbits >  fmaxbits and fbits < DKbits
 *
 */
inline CLCiphertext& Variant::encrypt(CLPlaintext& clpt, CLPublickey& clpk)
{
    // Get pk, m
    clpk.get(UB, con, Dcon, G, H);
    clpk.get(n, t);
    clpt.get(m, zm);

    // Select r from [0, Bcon-1]
    if (ebits == 0)
        e1 = RandomBnd(zUB);
    else
        e1 = RandomLen_ZZ(ebits);

    // Compute C1 = G^r ; has D = DK
    G.setD(DK);
    qfeG.power(T1, G, e1);
    
    // Compute C2 = F^m * psi(H^r)
    qfeH.power(T4, H, e1);
    sendQFtosmallerord(T5, temp1, T4, DK, Dcon, con); // phi_inv(H^r)
    qfe.initialize(T5, zcon);
    qfe.power(T3, T5, zcon); // psi(H^r)
    qfe.CLpower(T2, con, Dcon, m); // F^m
    mul(T3, T3, T2);
    
    // Set ciphertext
    ctt.set(T1, T3);
    
    // Return ciphertext
    return ctt;
}



inline CLPlaintext& Variant::decrypt(CLCiphertext& clct, CLPublickey& clpk, CLSecretkey& clsk)
{
    // Get pk, sk, ciphertexts
    clpk.get(con);
    clsk.get(zx);
    clct.get(C1, C2);
    
    // Compute M = C2 * psi(C1^-x)
    C1.setD(DK);
    qfe.initialize(C1, zx);
    qfe.power(M, C1, zx);
    inv(M, M); // C1^-x
    sendQFtosmallerord(T1, temp1, M, DK, Dcon, con);
    qfe.initialize(T1, zcon);
    qfe.power(M, T1, zcon); // psi(C1^-x)
    mul(M, M, C2); // C2 * psi(C1^-x)
    
    
    // Get decrypted m from M
    CL::solve(dm, M, con);
    
    // Set decrypted result
    ptt.set(dm);
    
    // Return decrypted plaintext
    return ptt;
}




inline CLPlaintext& Variant::mdecrypt(CLCiphertext& clct, CLPublickey& clpk, CLSecretkey& clsk)
{
    // Get pk, sk, ciphertexts
    clpk.get(con);
    clsk.get(zx);
    clct.get(C1, C2);
    
    // Compute M = C2 * psi(C1^-x)
    C1.setD(DK);
    qfe.initialize(C1, zx);
    qfe.power(M, C1, zx);
    inv(M, M); // C1^-x
    sendQFtosmallerord(T1, temp1, M, DK, Dcon, con);
    qfe.initialize(T1, zcon);
    qfe.power(M, T1, zcon); // psi(C1^-x)
    mul(M, M, C2); // C2 * psi(C1^-x)
    
    // Lift M to smaller order so that (f^2,f) is reduced
    sendQFtosmallerord(T1, temp1, M, Dcon, Dconsqr, con);
    qfe.initialize(T1, zcon);
    qfe.power(M, T1, zcon);
    
    // Get decrypted m from M
    CL::solve(dm, M, con);
    
    // Set decrypted result
    ptt.set(dm);
    
    // Return decrypted plaintext
    return ptt;
}


// ------------------------------------------------------------------------------
// ------------------------------------------------------------------------------
// ------------------------------------------------------------------------------
// ---------------------- Encryption and Decryption by DJS ----------------------
// ------------------------------------------------------------------------------
// ------------------------------------------------------------------------------
// ------------------------------------------------------------------------------


inline CLCiphertext& Variant::djsencrypt(CLPlaintext& clpt, CLPublickey& clpk)
{
    // Get pk, m
    clpk.get(UB, con, Dcon, G, H);
    clpk.get(n, t);
    clpt.get(m, zm);
    
    // Select r from [0, Bcon-1]
    if (ebits == 0)
        e1 = RandomBnd(zUB);
    else
        e1 = RandomLen_ZZ(ebits);
    
    // Compute C1 = G^r ; has D = DK
    G.setD(DK);
    qfeG.power(T1, G, e1);
    
    // Compute C2 = F^m * psi(H^r)
    qfeH.power(T4, H, e1);
    sendQFtosmallerord(T5, temp1, T4, DK, Dcon, con); // phi_inv(H^r)
    qfed.initialize(F, T5, zm, zcon);
    qfed.power(T3, F, T5, zm, zcon);
    
    // Set ciphertext
    ctt.set(T1, T3);
    
    // Return ciphertext
    return ctt;
}




inline CLPlaintext& Variant::djsdecrypt(CLCiphertext& clct, CLPublickey& clpk, CLSecretkey& clsk)
{
    // Get pk, sk, ciphertexts
    clpk.get(con);
    clsk.get(zx);
    clct.get(C1, C2);
    
    // Compute M = C2 * psi(C1^-x)
    C1.setD(DK);
    qfe.initialize(C1, zx);
    qfe.power(M, C1, zx);
    inv(M, M); // C1^-x
    sendQFtosmallerord(T1, temp1, M, DK, Dcon, con);
    qfe.initialize(T1, zcon);
    qfe.power(M, T1, zcon); // psi(C1^-x)
    mul(M, M, C2); // C2 * psi(C1^-x)
    
    
    
    // Get decrypted m from M
    switch(n)
    {
            case (1):
            CL::solve(dm, F, M, con, p, t);  // Prime power PH
            break;
                
            default:
            CL::solve(dm, F, M, con, pntn, pn, n, t);  // PH
            break;
    }
    
    // Set decrypted result
    ptt.set(dm);
    
    // Return decrypted plaintext
    return ptt;
}


// ------------------------------------------------------------------------------
// ------------------------------------------------------------------------------
// ------------------------------------------------------------------------------
// ----- CRT variants of Enc. and Dec. for CL and DJS functions given by DJS ----
// ------------------------------------------------------------------------------
// ------------------------------------------------------------------------------
// ------------------------------------------------------------------------------

/**
 *
 * Should be used when N > 1 and and fbits < fmaxbits
 *
 */
inline CLPlaintext& Variant::decryptCRT(CLCiphertext& clct, CLPublickey& clpk, CLSecretkey& clsk)
{
    // Get pk, sk, ciphertexts
    clpk.get(con);
    clsk.get(zx);
    clct.get(C1, C2);
    
    // Compute M = C2 * psi(C1^-x)
    C1.setD(DK);
    qfe.initialize(C1, zx);
    qfe.power(M, C1, zx);
    inv(M, M); // C1^-x
    sendQFtosmallerord(T1, temp1, M, DK, Dcon, con);
    qfe.initialize(T1, zcon);
    qfe.power(M, T1, zcon); // psi(C1^-x)
    mul(M, M, C2);
    
    // Get decrypted m from M
    CL::solveCRT1(dm, temp1, temp2, M, con, pntn, n, var);
    
    // Set decrypted result
    ptt.set(dm);
    
    // Return decrypted plaintext
    return ptt;
}


/**
 *
 * Decryption using CRT by lifting to a smaller order
 * Should be used when N > 1 and fbits > fmaxbits and fbits <= DKbits
 *
 */
inline CLPlaintext& Variant::mdecryptCRT(CLCiphertext& clct, CLPublickey& clpk, CLSecretkey& clsk)
{
    // Get pk, sk, ciphertexts
    clpk.get(con);
    clsk.get(zx);
    clct.get(C1, C2);
    
    // Compute M = C2 * psi(C1^-x)
    C1.setD(DK);
    qfe.initialize(C1, zx);
    qfe.power(M, C1, zx);
    inv(M, M); // C1^-x
    sendQFtosmallerord(T1, temp1, M, DK, Dcon, con);
    qfe.initialize(T1, zcon);
    qfe.power(M, T1, zcon); // psi(C1^-x)
    mul(M, M, C2); // C2 * psi(C1^-x)
    
    // Lift M to smaller order so that (f^2,f) is reduced
    sendQFtosmallerord(T1, temp1, M, Dcon, Dconsqr, con);
    qfe.initialize(T1, zcon);
    qfe.power(M, T1, zcon);
    
    // Get decrypted m from M
    CL::solveCRT1(dm, temp1, temp2, M, con, pntn, n, var);
    
    // Set decrypted result
    ptt.set(dm);
    
    // Return decrypted plaintext
    return ptt;
}



/**
 *
 * Modifies C2 = F^m * H^r to C2 = F1^m1 * ... * FN^mN * H^r
 * Should be used when N > 1 and fbits < fmaxbits
 *
 */
inline CLCiphertext& Variant::encrypt2(CLPlaintext& clpt, CLPublickey& clpk)
{
    // Get pk, m
    clpk.get(UB, con, Dcon, G, H, pn, n);
    clpk.get(n,t);
    clpt.get(m);
    
    // Compute m1 = m mod p1^t and m2 = m mod p2^t
    for (int i = 0; i < n; i++)
        mpz_mod(mn[i], m, pntn[i]);
    
    // Select r from [0, Bcon-1]
    if (ebits == 0)
        e1 = RandomBnd(zUB);
    else
        e1 = RandomLen_ZZ(ebits);
    
    // Compute C1 = G^r
    G.setD(DK);
    qfeG.power(T1, G, e1);
    
    // Compute C2 = F1^m1 * ... * Fn^mn * psi(H^r)
    qfeH.power(T4, H, e1);
    sendQFtosmallerord(T5, temp1, T4, DK, Dcon, con);
    qfe.initialize(T5, zcon);
    qfe.power(T2, T5, zcon); // psi(H^r)

    for (int i = 0; i < n; i++)
    {
        qfe.CLpower(Tn[i], pntn[i], Dcon, mn[i]);
        mul(T2, T2, Tn[i]);
    }
    
    // Set ciphertext
    ctt.set(T1, T2);
    
    // Return ciphertext
    return ctt;
}


inline CLPlaintext& Variant::decrypt2(CLCiphertext& clct, CLPublickey& clpk, CLSecretkey& clsk)
{
    // Get pk, sk, ciphertexts
    clpk.get(con, DK, Dcon, pn, n);
    clsk.get(zx);
    clct.get(C1, C2);
    
    // Compute M = C2/psi(C1^x) = (N^2, X)
    C1.setD(DK);
    qfe.initialize(C1, zx);
    qfe.power(M, C1, zx);
    inv(M, M); // C1^-x
    sendQFtosmallerord(T1, temp1, M, DK, Dcon, con);
    qfe.initialize(T1, zcon);
    qfe.power(M, T1, zcon); // psi(C1^-x)
    mul(M, M, C2);
    
    // Compute m from M
    CL::solveCRT2(dm, temp1, temp2, M, con, pntn, n, var);
    ptt.set(dm);
    
    // Return decrypted plaintext
    return ptt;
}



inline CLCiphertext& Variant::encrypt3(CLPlaintext& clpt, CLPublickey& clpk)
{
    // Get pk, m
    clpk.get(UB, con, DK, Dcon, F, G, H);
    clpk.get(n,t);
    clpt.get(m);
    
    // Compute mn = m mod pn^tn
    for (int i = 0; i < n; i++)
        mpz_mod(mn[i], m, pntn[i]);
    
    // Select r from [0, Bcon-1]
    if (ebits == 0)
        e1 = RandomBnd(zUB);
    else
        e1 = RandomLen_ZZ(ebits);
    
    // Compute C1 = G^r
    G.setD(DK);
    qfeG.power(T1, G, e1);
    
    // Compute Cn = Fn^mn * psi(H^r)
    qfeH.power(T4, H, e1);
    sendQFtosmallerord(T5, temp1, T4, DK, Dcon, con);
    qfe.initialize(T5, zcon);
    qfe.power(T2, T5, zcon); // psi(H^r)
    
    
    if (fbits < DKbits)
    {
        for (int i = 0; i < n; i++)
        {
            qfe.CLpower(Cn[i], pntn[i], Dcon, mn[i]);
            mul(Cn[i], Cn[i], T2);
        }
    }
    
    else
    {
        for (int i = 0; i < n; i++)
        {
            ZZ_limbs_set(e2, mn[i]->_mp_d, mn[i]->_mp_size);
            qfeFn[i].power(Cn[i], Fn[i], e2);
            mul(Cn[i], Cn[i], T2);
        }
    }
    
    // Set ciphertext
    ctt.set(T1, Cn, n);
    
    // Return ciphertext
    return ctt;
}



inline CLPlaintext& Variant::decrypt3(CLCiphertext& clct, CLPublickey& clpk, CLSecretkey& clsk)
{
    // Get pk, sk, ciphertexts
    clpk.get(F, con, pn, n, t);
    clsk.get(zx);
    clct.get(C1, Cn, n);
    
    // Compute psi(C1^-x)
    C1.setD(DK);
    qfe.initialize(C1, zx);
    qfe.power(CIX, C1, zx);
    inv(CIX, CIX);
    sendQFtosmallerord(T1, temp1, CIX, DK, Dcon, con);
    qfe.initialize(T1, zcon);
    qfe.power(CIX, T1, zcon); // CIX = psi(C1^-x)
    
    // Compute mn from Mn = Cn/psi(C1^x)
    if (fbits < DKbits)
    {
        for (int i = 0; i < n; i++)
        {
            mul(M, Cn[i], CIX);
            CL::solve(var[i], M, pntn[i]);
        }
    }
    
    else
    {
        for (int i = 0; i < n; i++)
        {
            mul(M, Cn[i], CIX);
            CL::solve(var[i], Fn[i], M, pntn[i], pn[i], t);
        }
    }
    
    // CRT
    CL::CRT_list(dm, temp1, temp2, n, var, pntn);
    ptt.set(dm);
    
    
    // Return decrypted plaintext
    return ptt;
}




// ------------------------------------------------------------------------------
// ------------------------------------------------------------------------------
// ------------------------------------------------------------------------------
// ----------------- HOMOMORPHIC PROPERTY VALIDATION FUNCTIONS ------------------
// ------------------------------------------------------------------------------
// ------------------------------------------------------------------------------
// ------------------------------------------------------------------------------


inline CLCiphertext& Variant::evalsum(CLCiphertext& clct1, CLCiphertext& clct2)
{
    clct1.get(T1, T2);
    clct2.get(T3, T4);

    T6.setD(Dcon);
    mul(T6, T2, T4);
    
    T5.setD(DK);
    mul(T5, T1, T3);
    
    // Select r from [0, Bcon-1]
    if (ebits == 0)
        e1 = RandomBnd(zUB);
    else
        e1 = RandomLen_ZZ(ebits);
        
    // Compute G^r
    qfeG.power(T1, G, e1);

    mul(T5, T5, T1);
    
    // Compute psi(H^r)
    qfeH.power(T2, H, e1);
    sendQFtosmallerord(T4, temp1, T2, DK, Dcon, con);
    qfe.initialize(T4, zcon);
    qfe.power(T2, T4, zcon);
        
    mul(T6, T6, T2);
    
    // Set ciphertext
    ctt.set(T5, T6);
        
    // Return ciphertext
    return ctt;
}


inline CLCiphertext& Variant::evalscal(CLCiphertext& clct, const mpz_t alpha)
{
    clct.get(T1, T2);
    ZZ_limbs_set(e1, alpha[0]._mp_d, alpha[0]._mp_size);
    
    T1.setD(DK);
    qfe.initialize(T1, e1);
    qfe.power(T3, T1, e1);
    
    T2.setD(Dcon);
    qfe2.initialize(T2, e1);
    qfe2.power(T4, T2, e1);
    
    // Set ciphertext
    ctt.set(T3, T4);
    
    // Return ciphertext
    return ctt;
}



inline CLCiphertext& Variant::evalsum3(CLCiphertext& clct1, CLCiphertext& clct2)
{
    clct1.get(T1, Tn1, n);
    clct2.get(T3, Tn2, n);
    
    T5.setD(DK);
    mul(T5, T1, T3);

    // Select r from [0, Bcon-1]
    if (ebits == 0)
        e1 = RandomBnd(zUB);
    else
        e1 = RandomLen_ZZ(ebits);

    // Compute G^r
    qfeG.power(T1, G, e1);
    mul(T5, T5, T1);

    qfeH.power(T2, H, e1);
    
    // Compute psi(H^r)
    T4.setD(Dcon);
    sendQFtosmallerord(T4, temp1, T2, DK, Dcon, con);
    qfe.initialize(T4, zcon);
    qfe.power(T6, T4, zcon);
    for (int i = 0; i < n; i++)
    {
        mul(Tn3[i], Tn1[i], Tn2[i]);
    }
    
    for (int i = 0; i < n; i++)
    {
        mul(Tn3[i], Tn3[i], T6);
    }
    
    // Set ciphertext
    ctt.set(T5, Tn3, n);
    
    // Return ciphertext
    return ctt;
}


inline CLCiphertext& Variant::evalscal3(CLCiphertext& clct, const mpz_t alpha)
{
    clct.get(T1, Tn1, n);
    ZZ_limbs_set(e1, alpha[0]._mp_d, alpha[0]._mp_size);
    
    T3.setD(DK);
    qfe.initialize(T1, e1);
    qfe.power(T3, T1, e1);
    T3.setD(Dcon);
    for (int i = 0; i < n; i++)
    {
        qfeTn[i].initialize(Tn1[i], e1);
        qfeTn[i].power(Tn2[i], Tn1[i], e1); // Fn^mn
    }
    
    ctt.set(T3, Tn2, n);
    return ctt;
}

#endif // Variant_HPP

