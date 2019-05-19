/**
 *
 * AUTHOR - Parthasarathi Das
 * BRIEF  - A class for the Paillier (basic) cryptosystem
 */

#ifndef Pai0_HPP
#define Pai0_HPP

#include "Pai.hpp"

class Pai0 : public Pai
{
    public:
    Pai0(){};
    Pai0(const int MS, const int switcher) : Pai(MS, switcher) {};
    ~Pai0() {};
        
    void           keygen  (PaiPublickey& paipk, PaiSecretkey& paisk);
    PaiCiphertext& encrypt (PaiPlaintext& paipt, PaiPublickey& paipk);
    PaiPlaintext&  decrypt (PaiCiphertext& paict, PaiPublickey& paipk, PaiSecretkey& paisk);
        
}; // end of class


// Out of line class member definitions
inline void Pai0::keygen(PaiPublickey& paipk, PaiSecretkey& paisk)
{
    Pai::keygen(paipk, paisk);
    
    if (crt_switch == 1)
    {
        // Pre-compute hp
        mpz_sub_ui(temp1, p, 1);
        mpz_mul(temp2, p, p);
        mpz_powm(hp, g, temp1, temp2); // g^(p-1) % p^2
        mpz_sub_ui(hp, hp, 1);
        mpz_divexact(hp, hp, p);
        buffer = mpz_invert(hp, hp, p);
        
        // Pre-compute hq
        mpz_sub_ui(temp1, q, 1);
        mpz_mul(temp3, q, q);
        mpz_powm(hq, g, temp1, temp3); // g^(q-1) % q^2
        mpz_sub_ui(hq, hq, 1);
        mpz_divexact(hq, hq, q);
        buffer = mpz_invert(hq, hq, q);
        
        paisk.set(p, q, temp2, temp3);
    }

    else
    {
        // Precompute L(g^lam % n^2)^-1
        mpz_powm(temp1, g, lam, nn);
        mpz_sub_ui(temp1, temp1, 1);
        mpz_divexact(temp1, temp1, n);
        buffer = mpz_invert(precomp, temp1, n);
    }

    // Set keys
    paipk.set(n, nn, g);
    paisk.set(lam);
}

inline PaiCiphertext& Pai0::encrypt(PaiPlaintext& paipt, PaiPublickey& paipk)
{
    // Get pt and pk
    paipt.get(m);
    paipk.get(n, nn, g);
    
    // Choose r from [0, n-1]
    mpz_urandomm(r, rands, n);

    // Compute ciphertext c = (g^m * r^n) % n^2
    mpz_powm(c, g, m, nn);
    mpz_powm(temp1, r, n, nn);
    mpz_mul(c, c, temp1);
    mpz_mod(c, c, nn);
        
    // Set ciphertext
    ct.set(c);
    return ct;
}

inline PaiPlaintext& Pai0::decrypt(PaiCiphertext& paict, PaiPublickey& paipk, PaiSecretkey& paisk)
{
    // Get c, pk and sk
    paict.get(c);
    
    if (crt_switch == 1)
    {
        paisk.get(p, q, temp2, temp3);
        
        // Compute mp
        mpz_sub_ui(temp1, p, 1);
        mpz_powm(mp, c, temp1, temp2);
        mpz_sub_ui(mp, mp, 1);
        mpz_divexact(mp, mp, p);
        mpz_mul(mp, mp, hp);
        mpz_mod(mp, mp, p);
        
        // Compute mq
        mpz_sub_ui(temp1, q, 1);
        mpz_powm(mq, c, temp1, temp3);
        mpz_sub_ui(mq, mq, 1);
        mpz_divexact(mq, mq, q);
        mpz_mul(mq, mq, hq);
        mpz_mod(mq, mq, q);
        
        // Compute m from mp and mq
        Pai::CRT(d, mp, mq, p, q);
    }
    
    else
    {
        paipk.get(n, nn);
        paisk.get(lam);
        
        // Compute L(c^lam % n^2)
        mpz_powm(temp1, c, lam, nn);
        mpz_sub_ui(temp1, temp1, 1);
        mpz_divexact(temp1, temp1, n);
    
        // Compute d = (L(.) * precomp) % n
        mpz_mul(temp1, temp1, precomp);
        mpz_mod(d, temp1, n);
    }
    
    // Set plaintext
    pt.set(d);
    return pt;
}

#endif // Pai0_HPP
