#ifndef Basic_HPP
#define Basic_HPP

#include "CL.hpp"


class Basic : public CL
{
    public:
    
    Basic() {};
    Basic (const int WIN, const int MS, const int DS, const int ES) : CL (WIN, MS, DS, ES) {};
    ~Basic() {};

    void          keygen      (CLPublickey& clpk,  CLSecretkey& clsk);
    
    // CL
    CLCiphertext& encrypt     (CLPlaintext& clpt,  CLPublickey& clpk);
    CLPlaintext&  decrypt     (CLCiphertext& clct, CLPublickey& clpk, CLSecretkey& clsk);
    CLPlaintext&  mdecrypt    (CLCiphertext& clct, CLPublickey& clpk, CLSecretkey& clsk);
    
    // DJS
    CLCiphertext& djsencrypt  (CLPlaintext& clpt,  CLPublickey& clpk);
    CLPlaintext&  djsdecrypt  (CLCiphertext& clct, CLPublickey& clpk, CLSecretkey& clsk);
    
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

inline void Basic::keygen(CLPublickey& clpk, CLSecretkey& clsk)
{
    // Set n and t parameters of public key
    clpk.set(n, t);
    
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
    while ( !( (mpz_legendre(DK, r) == 1) && (mpz_cmp_ui(temp1, 1) == 0)))
    {
        mpz_nextprime(r, r);
        mpz_gcd(temp1, r, con);
    }
    prime_r = mpz_get_ui(r);
    prime_index = prime_index_ge(prime_r);
    prime_index = R.next_prime_QF(prime_index, DK);
 
    // Compute ideal R^2 of discriminant DK
    sqr(RR, R);
    
    // Generate a random k in (Z/conZ)^*
    mpz_urandomm (k, rands, con);
    mpz_gcd(temp1, k, con);
    while ((mpz_cmp_ui(k,0) == 0) && (mpz_cmp_ui(temp1, 1) != 0))
    {
        mpz_urandomm (k, rands, con);
        mpz_gcd(temp1, k, con);
    }
    
    // Compute phi_inv(R^2)
    sendQFtosmallerord(psi_RR, temp1, RR, DK, Dcon, con);
    
    // Set G = phi_inv(R^2)^f * (F^k)
    ZZ_limbs_set(e1, k[0]._mp_d, k[0]._mp_size);
    qfed.initialize(psi_RR, F, zcon, e1);
    qfed.power(G, psi_RR, F, zcon, e1);
    
    // Generate secret key x in [0, Bcon-1] = [0,UB)
    if (ebits == 0)
        zx = RandomBnd(zUB);
    else
        zx = RandomLen_ZZ(ebits);
    
    // Compute H = G^x
    qfeG.initialize(G, zx);
    qfeG.power(H, G, zx);
    
    // Set Fi's for N > 1
    for (int i = 0; i < n; i++)
    {
        mpz_mul(temp1, pntn[i], pntn[i]);
        Fn[i].assign(temp1, pntn[i], Dcon);
        
        qfeFn[i].initialize(Fn[i], e1);
    }
    
    // Set exponentiation objects for F and H bases for use in encrypt
    qfed.initialize(F, H, e1, e1); // Double exp.
    qfeF.initialize(F, e1); // F -- uses qfeF
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
inline CLCiphertext& Basic::encrypt(CLPlaintext& clpt, CLPublickey& clpk)
{
    // Get pk, m
    clpk.get(UB, con, Dcon, G, H);
    clpt.get(m, zm);
    
    // Select r from [0, Bcon-1]
    if (ebits == 0)
        e1 = RandomBnd(zUB);
    else
        e1 = RandomLen_ZZ(ebits);
    
    // Compute C1 = G^r
    qfeG.power(T1, G, e1);
    
    // Compute C2 = F^m * H^r
    qfe.CLpower(T2, con, Dcon, m); // F^m
    qfeH.power(T3, H, e1); // H^r
    mul(T3, T3, T2); // F^m * H^r
    
    // Set ciphertext C1, C2
    ctt.set(T1, T3);
    
    // Return ciphertext
    return ctt;
}


/**
 *
 * Decryption using CL method
 * Should be used when fbits <= fmaxbits
 *
 */
inline CLPlaintext& Basic::decrypt(CLCiphertext& clct, CLPublickey& clpk, CLSecretkey& clsk)
{
    // Get pk, sk, ciphertexts
    clpk.get(con);
    clsk.get(zx);
    clct.get(C1, C2);
 
    // Compute M = C2/C1^x
    qfe.initialize(C1, zx);
    qfe.power(CIX, C1, zx);
    inv(CIX, CIX);
    mul(M, C2, CIX);
    
    // Get decrypted m from M - Modular inverse decryption
    CL::solve(dm, M, con);
    
    // Set decrypted result
    ptt.set(dm);
     
    // Return decrypted plaintext
    return ptt;
}


/**
 *
 * Decryption by lifting to a smaller order
 * Should be used when fbits > fmaxbits and fbits < DKbits
 *
 */
inline CLPlaintext& Basic::mdecrypt(CLCiphertext& clct, CLPublickey& clpk, CLSecretkey& clsk)
{
    // Get pk, sk, ciphertexts
    clpk.get(con);
    clsk.get(zx);
    clct.get(C1, C2);
    
    // Compute M = C2/C1^x
    qfe.initialize(C1, zx);
    qfe.power(CIX, C1, zx);
    inv(CIX, CIX);
    mul(M, C2, CIX);
    
    // Lift M to smaller order so that (f^2,f) is reduced
    sendQFtosmallerord(T1, temp1, M, Dcon, Dconsqr, con);
    qfe.initialize(T1, zcon);
    qfe.power(M, T1, zcon);
    
    // Get decrypted m from M - Modular inverse decryption
    CL::solve(dm, M, con);
    T1.setD(Dcon);
    
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


inline CLCiphertext& Basic::djsencrypt(CLPlaintext& clpt, CLPublickey& clpk)
{
    // Get pk, m
    clpk.get(UB, con, Dcon, G, H);
    clpt.get(m, zm);
    
    // Select r from [0, Bcon-1]
    if (ebits == 0)
        e1 = RandomBnd(zUB);
    else
        e1 = RandomLen_ZZ(ebits);
    
    // Compute C1 = G^r
    qfeG.power(T1, G, e1);
    
    // Compute C2 = F^m * H^r
    qfed.power(T3, F, H, zm, e1);
    
    // Set ciphertext C1, C2
    ctt.set(T1, T3);
    
    // Return ciphertext
    return ctt;
}


inline CLPlaintext& Basic::djsdecrypt(CLCiphertext& clct, CLPublickey& clpk, CLSecretkey& clsk)
{
    // Get pk, sk, ciphertexts
    clpk.get(con);
    clsk.get(zx);
    clct.get(C1, C2);
    
    // Compute M = C2/C1^x
    qfe.initialize(C1, zx);
    qfe.power(CIX, C1, zx);
    inv(CIX, CIX);
    mul(M, C2, CIX);
    
    switch(n)
    {
            case (1): // Prime power PH
            CL::solve(dm, F, M, con, p, t);
            break;
                
            default: // Generic PH
            CL::solve(dm, F, M, con, pntn, pn, n, t);
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
inline CLPlaintext& Basic::decryptCRT(CLCiphertext& clct, CLPublickey& clpk, CLSecretkey& clsk)
{
    // Get pk, sk, ciphertexts
    clpk.get(con);
    clsk.get(zx);
    clct.get(C1, C2);
    
    // Compute M = C2/C1^x
    qfe.initialize(C1, zx);
    qfe.power(CIX, C1, zx);
    inv(CIX, CIX);
    mul(M, C2, CIX);
    
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
 * Should be used when N > 1 and fbits > fmaxbits and fbits < DKbits
 *
 */
inline CLPlaintext& Basic::mdecryptCRT(CLCiphertext& clct, CLPublickey& clpk, CLSecretkey& clsk)
{
    // Get pk, sk, ciphertexts
    clpk.get(con);
    clsk.get(zx);
    clct.get(C1, C2);
    
    // Compute M = C2/C1^x
    qfe.initialize(C1, zx);
    qfe.power(CIX, C1, zx);
    inv(CIX, CIX);
    mul(M, C2, CIX);
    
    // Lift M to smaller order so that (f^2,f) is reduced
    sendQFtosmallerord(T1, temp1, M, Dcon, Dconsqr, con);
    qfe.initialize(T1, zcon);
    qfe.power(M, T1, zcon);
    
    // Get decrypted m from M - Modular inverse decryption
    CL::solveCRT1(dm, temp1, temp2, M, con, pntn, n, var);
    T1.setD(Dcon);
    
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
inline CLCiphertext& Basic::encrypt2(CLPlaintext& clpt, CLPublickey& clpk)
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
    qfeG.power(T1, G, e1);

    // Compute C2 = F1^m1 * ... * Fn^mn * H^r
    qfeH.power(T4, H, e1);
    for (int i = 0; i < n; i++)
    {
        qfe.CLpower(Tn[i], pntn[i], Dcon, mn[i]);
        mul(T4, T4, Tn[i]);
    }
    
    // Set ciphertext
    ctt.set(T1, T4);

    // Return ciphertext
    return ctt;
}


inline CLPlaintext& Basic::decrypt2(CLCiphertext& clct, CLPublickey& clpk, CLSecretkey& clsk)
{
    // Get pk, sk, ciphertexts
    clpk.get(con, pn, n);
    clsk.get(zx);
    clct.get(C1, C2);
    
    // Compute M = C2/C1^x = (N^2, X)
    qfe.initialize(C1, zx);
    qfe.power(CIX, C1, zx);
    inv(CIX, CIX);
    mul(M, C2, CIX);
    
    // Compute m from M
    CL::solveCRT2(dm, temp1, temp2, M, con, pntn, n, var);
    ptt.set(dm);
    
    // Return decrypted plaintext
    return ptt;
}


/**
 *
 */
inline CLCiphertext& Basic::encrypt3(CLPlaintext& clpt, CLPublickey& clpk)
{
    // Get pk, m
    clpk.get(UB, con, F, G, H);
    clpk.get(n,t);
    clpt.get(m);
    
    // Compute mn = m mod pn^tn
    for (int i = 0; i < n; i++)
    {
        mpz_mod(mn[i], m, pntn[i]);
    }
    
    
    // Select r from [0, Bcon-1]
    if (ebits == 0)
        e1 = RandomBnd(zUB);
    else
        e1 = RandomLen_ZZ(ebits);
    
    
    // Compute C1 = G^r
    qfeG.power(T1, G, e1);
    
    
    // Compute Cn = Fn^mn * H^r
    qfeH.power(T2, H, e1);
    
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


inline CLPlaintext& Basic::decrypt3(CLCiphertext& clct, CLPublickey& clpk, CLSecretkey& clsk)
{
    // Get pk, sk, ciphertexts
    clpk.get(F, con, pn, n, t);
    clsk.get(zx);
    clct.get(C1, Cn, n);
    
    // Compute C1^-x
    qfe.initialize(C1, zx);
    qfe.power(CIX, C1, zx);
    inv(CIX, CIX);
    
    // Compute mn from Mn = Cn/C1^x
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

inline CLCiphertext& Basic::evalsum(CLCiphertext& clct1, CLCiphertext& clct2)
{
    clct1.get(T1, T2);
    clct2.get(T3, T4);
    
    //
    mul(T5, T1, T3);
    mul(T6, T2, T4);
    
    // Select r from [0, Bcon-1]
    if (ebits == 0)
        e1 = RandomBnd(zUB);
    else
        e1 = RandomLen_ZZ(ebits);
        
    // Compute G^r
    qfeG.power(T1, G, e1);
    
    // Compute H^r
    qfeH.power(T2, H, e1);
    
    //
    mul(T5, T5, T1);
    mul(T6, T6, T2);
    
    // Set ciphertext
    ctt.set(T5, T6);
        
    // Return ciphertext
    return ctt;
}


inline CLCiphertext& Basic::evalscal(CLCiphertext& clct, const mpz_t alpha)
{
    clct.get(T1, T2);
    ZZ_limbs_set(e1, alpha[0]._mp_d, alpha[0]._mp_size);
    
    qfe .initialize(T1, e1);
    qfe2.initialize(T2, e1);
    
    qfe .power(T3, T1, e1);
    qfe2.power(T4, T2, e1);
    
    // Set ciphertext
    ctt.set(T3, T4);
    
    // Return ciphertext
    return ctt;
}


inline CLCiphertext& Basic::evalsum3(CLCiphertext& clct1, CLCiphertext& clct2)
{
    
    clct1.get(T1, Tn1, n);
    clct2.get(T3, Tn2, n);
    
    mul(T5, T1, T3);
    for (int i = 0; i < n; i++)
    {
        mul(Tn3[i], Tn1[i], Tn2[i]);
    }
    
    // Select r from [0, Bcon-1]
    if (ebits == 0)
        e1 = RandomBnd(zUB);
    else
        e1 = RandomLen_ZZ(ebits);
    
    // Compute G^r
    qfeG.power(T1, G, e1);
    
    // Compute H^r
    qfeH.power(T2, H, e1);
    
    //
    mul(T5, T5, T1);
    
    for (int i = 0; i < n; i++)
    {
        mul(Tn3[i], Tn3[i], T2);
    }
    
    // Set ciphertext
    ctt.set(T5, Tn3, n);
    
    // Return ciphertext
    return ctt;
}


inline CLCiphertext& Basic::evalscal3(CLCiphertext& clct, const mpz_t alpha)
{
    clct.get(T1, Tn1, n);
    ZZ_limbs_set(e1, alpha[0]._mp_d, alpha[0]._mp_size);
    
    qfe.initialize(T1, e1);
    qfe.power(T3, T1, e1);
    for (int i = 0; i < n; i++)
    {
        qfeTn[i].initialize(Tn1[i], e1);
        qfeTn[i].power(Tn2[i], Tn1[i], e1); // Fn^mn
    }
    
    ctt.set(T3, Tn2, n);
    return ctt;
}

#endif // Basic_HPP
