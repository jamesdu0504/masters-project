/**
 *
 * AUTHOR - Parthasarathi Das
 * BRIEF  - A class for the BCP cryptosystem
 */

#ifndef BCP_HPP
#define BCP_HPP

#include "../../Variables/BCP/BCPPublickey.hpp"
#include "../../Variables/BCP/BCPSecretkey.hpp"
#include "../../Variables/BCP/BCPPlaintext.hpp"
#include "../../Variables/BCP/BCPCiphertext.hpp"

class BCP
{
    public:
    BCP(){};
    BCP(const int MS, const int ES);
    ~BCP();
    
    void           init     (const mpz_t modulus, const mpz_t prime1, const mpz_t prime2);
    void           keygen   (BCPPublickey& bcppk, BCPSecretkey& bcpsk);
    BCPCiphertext& encrypt  (BCPPlaintext& bcppt, BCPPublickey& bcppk);
    BCPPlaintext&  decrypt  (BCPCiphertext& bcpct, BCPPublickey& bcppk, BCPSecretkey& bcpsk);
    BCPCiphertext& evalsum  (BCPPublickey& bcppk, BCPCiphertext& bcpct1, BCPCiphertext& bcpct2);
    BCPCiphertext& evalscal (BCPPublickey& bcppk, BCPCiphertext& bcpct, const mpz_t alpha);
    
    protected:
    int              buffer, modbits, ebits, bits;
    mpz_t            a, A, A1, A2, B, B1, B2, d, g, h, m, N, p, q, r, Nsqr, ordG, pbar, qbar, alpha, seed, temp;
    gmp_randstate_t  rands;
    BCPPlaintext     pt;
    BCPCiphertext    ct;
            
            
}; // end of class


// Out of line class member definitions
inline BCP::BCP(const int MS, const int ES)
{
    modbits = MS;
    ebits = ES;
    mpz_inits(a, A, A1, A2, B, B1, B2, d, g, h, m, N, p, q, r, Nsqr, ordG, pbar, qbar, alpha, seed, temp, NULL);
    gmp_randinit_default(rands);
}

inline BCP::~BCP()
{
    mpz_clears(a, A, A1, A2, B, B1, B2, d, g, h, m, N, p, q, r, Nsqr, ordG, pbar, qbar, alpha, seed, temp, NULL);
    gmp_randclear(rands);
}

inline void BCP::init(const mpz_t modulus, const mpz_t prime1, const mpz_t prime2)
{
    mpz_set(N, modulus);
    mpz_set(p, prime1);
    mpz_set(q, prime2);
}

inline void BCP::keygen(BCPPublickey& bcppk, BCPSecretkey& bcpsk)
{
    // Compute pbar = (p - 1) / 2 and qbar = (q - 1) / 2
    mpz_sub_ui(pbar, p, 1);            // pbar = p - 1
    mpz_divexact_ui(pbar, pbar, 2);    // pbar = pbar / 2
    mpz_sub_ui(qbar, q, 1);            // qbar = q - 1
    mpz_divexact_ui(qbar, qbar, 2);    // qbar = qbar / 2
    
    // Compute ord(G) = N * pbar * qbar
    mpz_mul(ordG, pbar, qbar);
    mpz_mul(ordG, N, ordG);
    
    // Choose alpha
    mpz_mul(Nsqr, N, N);
    mpz_urandomm(alpha, rands, Nsqr);
    mpz_gcd (temp, alpha, Nsqr);
    while (mpz_cmp_ui(temp,1) != 0)
    {
        mpz_urandomm(alpha, rands, Nsqr);
        mpz_gcd (temp, alpha, Nsqr);
    }

    // Compute g
    mpz_powm_ui(g, alpha, 2, Nsqr);
    
    // Compute sk a
    mpz_add_ui(temp, ordG, 1);
    if (ebits == 0)
    {
        while (mpz_cmp_ui(a, 0) == 0)
            mpz_urandomm(a, rands, temp);
    }
    else
    {
        while (mpz_cmp_ui(a, 0)== 0)
            mpz_rrandomb(a, rands, ebits);
    }
    
    // Compute h
    mpz_powm(h, g, a, Nsqr);
    
    // Set keys
    bcppk.set(N, Nsqr, g, h);
    bcpsk.set(a);
}

inline BCPCiphertext& BCP::encrypt(BCPPlaintext& bcppt, BCPPublickey& bcppk)
{
    bcppt.get(m);
    bcppk.get(N, Nsqr, g, h);
    
    // Choose r
    if (ebits == 0)
    {
        while (mpz_cmp_ui(r, 0)== 0)
            mpz_urandomm(r, rands, temp);
    }
    else
    {
        while (mpz_cmp_ui(r, 0)== 0)
            mpz_rrandomb(r, rands, ebits);
    }
    
    // Compute  ciphertext A = g^r % N^2
    mpz_powm(A, g, r, Nsqr);
    
    // Compute ciphertext B = h^r * (1 + mn) % N^2 (and encrypt m)
    mpz_powm(B, h, r, Nsqr);
    mpz_mul(temp, m, N);        // temp = m * N
    mpz_add_ui(temp, temp, 1);  // temp = (m * N) + 1
    mpz_mul(B, B, temp);        // B = B * ((m * N) + 1)
    mpz_mod(B, B, Nsqr);        // B = (B * ((m * N) + 1)) % (N^2)
    
    ct.set(A, B);
    return ct;
}

inline BCPPlaintext& BCP::decrypt(BCPCiphertext& bcpct, BCPPublickey& bcppk, BCPSecretkey& bcpsk)
{
    bcpct.get(A, B);
    bcpsk.get(a);
    bcppk.get(N, Nsqr);
    
    mpz_powm(temp, A, a, Nsqr); // temp = A^a % N^2
    buffer = mpz_invert(temp, temp, Nsqr);
    mpz_mul(temp, B, temp);
    mpz_sub_ui(temp, temp, 1);
    mpz_mod(temp, temp, Nsqr);
    mpz_divexact(d, temp, N);
    
    pt.set(d);
    return pt;
}

inline BCPCiphertext& BCP::evalsum(BCPPublickey& bcppk, BCPCiphertext& bcpct1, BCPCiphertext& bcpct2)
{
    bcppk.get(N, Nsqr);
    bcpct1.get(A1, B1);
    bcpct2.get(A2, B2);
    
    mpz_mul(A, A1, A2);
    mpz_mul(B, B1, B2);
    mpz_mod(A, A, Nsqr);
    mpz_mod(B, B, Nsqr);
    
    ct.set(A, B);
    return ct;
}

inline BCPCiphertext& BCP::evalscal(BCPPublickey& bcppk, BCPCiphertext& bcpct, const mpz_t scalar)
{
    bcppk.get(N, Nsqr);
    bcpct.get(A, B);
    
    mpz_powm(A1, A, scalar, Nsqr);
    mpz_powm(B1, B, scalar, Nsqr);
    
    ct.set(A1, B1);
    return ct;
}

#endif // BCP_HPP
