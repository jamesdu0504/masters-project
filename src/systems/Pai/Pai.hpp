/**
 *
 * AUTHOR - Parthasarathi Das
 * BRIEF  - A class for the Paillier cryptosystem
 */

#ifndef Pai_HPP
#define Pai_HPP

#include "../../Variables/Pai/PaiPublickey.hpp"
#include "../../Variables/Pai/PaiSecretkey.hpp"
#include "../../Variables/Pai/PaiPlaintext.hpp"
#include "../../Variables/Pai/PaiCiphertext.hpp"

class Pai
{
    public:
    Pai(){};
    Pai(const int MS, const int switcher);
    ~Pai();
    
    void                    init     (const mpz_t modulus, const mpz_t prime1, const mpz_t prime2);
    void                    CRT      (mpz_t x, const mpz_t a, const mpz_t b, const mpz_t m, const mpz_t n);
    void                    keygen   (PaiPublickey& paipk, PaiSecretkey& paisk);
    virtual PaiCiphertext&  encrypt  (PaiPlaintext& paipt, PaiPublickey& paipk) = 0;
    virtual PaiPlaintext&   decrypt  (PaiCiphertext& paict, PaiPublickey& paipk, PaiSecretkey& paisk) = 0;
    PaiCiphertext&          evalsum  (PaiPublickey& paipk, PaiCiphertext& paict1, PaiCiphertext& paict2);
    PaiCiphertext&          evalscal (PaiPublickey& paipk, PaiCiphertext& paict, const mpz_t alpha);
        
    protected:
    int              buffer, modbits, bits, crt_switch;
    mpz_t            p, q, n, nn, lam, alpha, r, g, c, m, d, seed, precomp, hp, hq, mp, mq, crts, crtt, crtq, crtr, temp1, temp2, temp3;
    gmp_randstate_t  rands;
    PaiPlaintext     pt;
    PaiCiphertext    ct;
};



// Out of line class member definitions
inline Pai::Pai(const int MS, const int switcher)
{
    modbits = MS;
    crt_switch = switcher;
    
    mpz_inits(p, q, n, nn, lam, alpha, r, g, c, m, d, seed, precomp, hp, hq, mp, mq, temp1, temp2, temp3, NULL);
    mpz_inits(crts, crtt, crtq, crtr, NULL);
    gmp_randinit_default(rands);
}

inline Pai::~Pai()
{
    mpz_clears(p, q, n, nn, lam, alpha, r, g, c, m, d, seed, precomp, hp, hq, mp, mq, temp1, temp2, temp3, NULL);
    mpz_clears(crts, crtt, crtq, crtr, NULL);
    gmp_randclear(rands);
    
}

inline void Pai::init(const mpz_t modulus, const mpz_t prime1, const mpz_t prime2)
{
    mpz_set(n, modulus);
    mpz_set(p, prime1);
    mpz_set(q, prime2);
}

inline void Pai::CRT(mpz_t x, const mpz_t a, const mpz_t b, const mpz_t m, const mpz_t n)
{
    mpz_gcdext(x, crts, crtt, m, n);
    mpz_sub(crtt, b, a);
    mpz_fdiv_qr (crtq, crtr, crtt, x);
    mpz_lcm(crtt, m, n);
    mpz_mul(x, crtq, crts);
    mpz_mul(x, x, m);
    mpz_add(x, x, a);
    mpz_mod(x, x, crtt);
}

inline void Pai::keygen(PaiPublickey& paipk, PaiSecretkey& paisk)
{
    // Compute n^2
    mpz_mul(nn, n, n);
    
    // Compute lambda (Carmichael function of n)
    mpz_sub_ui(lam, p, 1);
    mpz_sub_ui(temp1, q, 1);
    mpz_lcm(lam, lam, temp1);
    
    // Set g = n + 1
    mpz_add_ui(g, n, 1);

    // Generate alpha from [1, lam]
    mpz_urandomm(alpha, rands, lam);
}

inline PaiCiphertext& Pai::evalsum(PaiPublickey& paipk, PaiCiphertext& paict1, PaiCiphertext& paict2)
{
    paipk.get(n, nn);
    paict1.get(temp1);
    paict2.get(temp2);
    
    mpz_mul(temp1, temp1, temp2);
    mpz_mod(temp1, temp1, nn);
    
    ct.set(temp1);
    return ct;
}

inline PaiCiphertext& Pai::evalscal(PaiPublickey& paipk, PaiCiphertext& paict, const mpz_t alpha)
{
    paipk.get(n, nn);
    paict.get(temp1);
    
    mpz_powm(temp1, temp1, alpha, nn);
    
    ct.set(temp1);
    return ct;
}

#endif // Pai_HPP
