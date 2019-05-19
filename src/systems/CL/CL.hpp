/**
 *
 * AUTHOR - Parthasarathi Das
 * BRIEF - A class for the CL cryptosystem
 *
 */

#ifndef CL_HPP
#define CL_HPP

#include "../../parameters/CL/CLPublickey.hpp"
#include "../../parameters/CL/CLSecretkey.hpp"
#include "../../parameters/CL/CLPlaintext.hpp"
#include "../../parameters/CL/CLCiphertext.hpp"
#include "../../QFExponentiation.hpp"


NTL_CLIENT

class CL
{
    public:
    CL  ();
    CL  (const int window, const int mbits, const int dkbits, const int expbits);
    ~CL ();

    void       initialise     (const int N, const int T, const mpz_t f, const mpz_t* pN,const mpz_t D);
    void       getB           (const mpz_t D, const int mu);
    
    void       CRT            (mpz_t x, mpz_t a, mpz_t b, mpz_t p1, mpz_t p2);
    void       CRT_list       (mpz_t x, mpz_t t1, mpz_t t2, const int n, mpz_t *v, mpz_t *moduli);
    
    void       solve          (mpz_t rop, QF& GE, const mpz_t con);
    void       solve          (mpz_t rop, QF& GG, const QF& GE, const mpz_t con, const mpz_t prime, const int exp);
    void       solve          (mpz_t rop, QF& F, QF& M, const mpz_t con, const mpz_t* pntn, const mpz_t* pn, const int n, const int t);
    
    void       solveCRT1      (mpz_t dm, mpz_t temp1, mpz_t temp2, QF& M, const mpz_t con, mpz_t *pn, const int n, mpz_t *var);
    void       solveCRT2      (mpz_t dm, mpz_t temp1, mpz_t temp2, QF& M, const mpz_t con, mpz_t *pn, const int n, mpz_t *var);
    
    static int check_legendre (mpz_t* pn, const mpz_t p, const int i, const int leg);
    
    
    
    protected:
    int                       prime_index;
    unsigned int              n, t, gamma, mu, fbits, DKbits, ebits, qbits, temp, prime_r, easy_max_flen;
    mpz_qform_t               form;
    mpz_qform_group_t         group;
    gmp_randstate_t           rands;
    mpz_t                     p, prod, q, r, k, B, UB, x, con, DK, Dcon, Dconsqr, seed, m, m1, m2, dm, dm1, dm2;
    mpz_t                     *pn, *pntn, *tn, *var, *mn;
    ZZ                        zp, zm, zx, zcon, zUB;
    QF                        F, F1, F2, G, H, psiH, R, RR, psi_RR, C1, C2, C3, M, M1;
    QFExponentiation          qfe, qfe2, qfe3, qfed, qfed2, qfeF, qfeF1, qfeF2, qfeG, qfeH, qfeM;
    vector<QF>                Fn, Tn, Tn1, Tn2, Tn3, Cn, Mn;
    vector<QFExponentiation>  qfeFn, qfeTn;
    CLPlaintext               ptt;
    CLCiphertext              ctt;
    int len;
    
    // Buffer variables
    mpz_t  temp1, temp2, temp3, temp4;
    QF     T1, T2, T3, T4, T5, T6, CIX;
    ZZ     e1, e2, e3, zrop;
            
    // Variables for setB()
    mpfr_t b, d, pi;

    // Variables for solve()
    mpz_t st, stemp1, stemp2;
    
    // Variables for Chinese Remainder-ing
    mpz_t crta, crtb, crts, crtt, crtq, crtr;

}; // end of class CL

inline CL::CL()
{
}

inline CL::CL(const int window, const int mbits, const int dkbits, const int expbits)
{
    mu     = mbits;
    fbits  = mbits;
    ebits  = expbits;
    DKbits = dkbits;
    len    = 8;

    gmp_randinit_default  (rands);
    mpz_qform_init        (&group, &form);
    mpz_qform_group_init  (&group);
    mpz_inits             (p, prod, q, r, k, B, UB, x, con, seed, DK, Dcon, Dconsqr, m, m1, m2, dm, dm1, dm2, NULL);
    mpz_inits             (temp1, temp2, temp3, temp4, NULL);
    mpz_inits             (crta, crtb, crts, crtt, crtq, crtr, NULL);
    

    // Initialize mpz arrays
    pn   = mpz_init_array (len);
    tn   = mpz_init_array (len);
    pntn = mpz_init_array (len);
    mn   = mpz_init_array (len);
    var  = mpz_init_array (len);
    
    
    //
    Fn    .resize (len);
    Tn    .resize (len);
    Tn1   .resize (len);
    Tn2   .resize (len);
    Tn3   .resize (len);
    Cn    .resize (len);
    Mn    .resize (len);
    qfeFn .resize (len);
    qfeTn .resize (len);

    
    // Initialize Ideal arrays
    for (int i = 0; i < len; i++)
    {
        Fn[i]  .init (group);
        Tn[i]  .init (group);
        Tn1[i] .init (group);
        Tn2[i] .init (group);
        Cn[i]  .init (group);
        Mn[i]  .init (group);
    }
    

    //
    R      .init   (group);
    RR     .init   (group);
    psiH   .init   (group);
    psi_RR .init   (group);
    C1     .init   (group);
    C2     .init   (group);
    C3     .init   (group);
    CIX    .init   (group);
    M      .init   (group);
    M1     .init   (group);
    F      .init   (group);
    F1     .init   (group);
    F2     .init   (group);
    G      .init   (group);
    H      .init   (group);
    T1     .init   (group);
    T2     .init   (group);
    T3     .init   (group);
    T4     .init   (group);
    T5     .init   (group);
    T6     .init   (group);
    ctt    .init   (group);
    
    
    //    
    qfe   .init  (window);
    qfe2  .init  (window);
    qfe3  .init  (window);
    qfed  .init  (window);
    qfed2 .init  (window);
    qfeF  .init  (window);
    qfeF1 .init  (window);
    qfeF2 .init  (window);
    qfeG  .init  (window);
    qfeH  .init  (window);
    qfeM  .init  (window);
    
    for (int i = 0; i < len; i++)
    {
        qfeFn[i].init (window);
        qfeTn[i].init (window);
    }
    
    
    // setB temps
    mpfr_init2 (b,  200);
    mpfr_init2 (d,  200);
    mpfr_init2 (pi, 200);

    
    // solve4
    mpz_inits (st, stemp1, stemp2, NULL);
}



inline CL::~CL()
{
    gmp_randclear          (rands);
    mpz_qform_clear        (&group, &form);
    mpz_qform_group_clear  (&group);
    mpz_clears             (p, prod, q, r, k, B, UB, x, con, seed, DK, Dcon, Dconsqr, m, m1, m2, dm, dm1, dm2, NULL);
    mpz_clears             (temp1, temp2, temp3, temp4, NULL);
    mpz_clears             (crta, crtb, crts, crtt, crtq, crtr, NULL);
    mpz_clear_array        (pn, len);
    mpz_clear_array        (tn, len);
    mpz_clear_array        (pntn, len);
    mpz_clear_array        (mn, len);
    mpz_clear_array        (var, len);
    
    mpfr_clear (b);
    mpfr_clear (d);
    mpfr_clear (pi);
    mpfr_free_cache ();

    mpz_clears (st, stemp1, stemp2, NULL);
}




inline void CL::initialise(const int N, const int T, const mpz_t f, const mpz_t* pN, const mpz_t D)
{
    // Set no of primes in conductor and the power in the conductor
    n = N;
    t = T;
    
    // Set conductor mpz_t and ZZ type
    mpz_set(con, f);
    ZZ_limbs_set(zcon, con[0]._mp_d, con[0]._mp_size);
    
    // Set fundamental discriminant
    mpz_set(DK, D);
    
    // Set primes and prime powers of conductor
    for (int i = 0; i < n; i++)
    {
        mpz_set(pn[i], pN[i]);
        mpz_pow_ui(pntn[i], pn[i], t);
    }

    // Set p if conductor is a single prime power
    if (n == 1)
        mpz_set(p, pn[0]);
}


/**
 *
 * Computes an upper bound B on odd part of the class number s
 * Source - Page 12, https://eprint.iacr.org/2015/047.pdf
 *
 */
inline void CL::getB(const mpz_t D, const int mu)
{
    int i;
    
    i = mpfr_const_pi(pi, MPFR_RNDU);       // Set pi value to pi
    i = mpfr_set_z(d, D, MPFR_RNDU);        // Get DK as mpfr_t type
    i = mpfr_log(b, d, MPFR_RNDU);          // 1. Compute log(|DK|)
    i = mpfr_sqrt(d, d, MPFR_RNDU);         // 2. Compute |DK|^1/2
    i = mpfr_mul(b, b, d, MPFR_RNDU);       // 3. Multiply b & d
    i = mpfr_div(b, b, pi, MPFR_RNDU);      // 4. Divide b by pi
    i = mpfr_div_ui(b, b, 4, MPFR_RNDU);    // 5. Divide b by 4
    i = mpfr_get_z(B, b, MPFR_RNDU);        // Set B = b
    
    mpz_mul_2exp(B, B, mu);                 // Set B = 2^mu * B
}


/**
 *
 * Checks if the Legendre symbols (pn[j]/p) == (p/pn[j]) == leg
 * Returns 1 if true, 0 otherwise
 *
 */
inline int CL::check_legendre(mpz_t* pn, const mpz_t p, const int i, const int leg)
{
    int j = 0;
    while (j < i)
    {
        if (!((mpz_legendre(pn[j], p) == leg) && (mpz_legendre(p, pn[j]) == leg)))
        {
            return 0;
        }

        j++;
    }
    return 1;
}


/**
 *
 * Solves the two simultaneous congruences, x = a (mod p1) and x = b (mod p2)
 * Source - http://doc.sagemath.org/html/en/reference/rings_standard/sage/arith/misc.html
 *
 */
inline void CL::CRT(mpz_t x, mpz_t a, mpz_t b, mpz_t p1, mpz_t p2)
{
    mpz_gcdext(x, crts, crtt, p1, p2);
    mpz_sub(crtt, b, a);
    mpz_fdiv_qr(crtq, crtr, crtt, x);
    mpz_lcm(crtt, p1, p2);
    mpz_mul(x, crtq, crts);
    mpz_mul(x, x, p1);
    mpz_add(x, x, a);
    mpz_mod(x, x, crtt);
}


/**
 *
 * Solves the simultaneous congruences, x = ai (mod pi)
 * Source - http://doc.sagemath.org/html/en/reference/rings_standard/sage/arith/misc.html
 *
 */
inline void CL::CRT_list(mpz_t x, mpz_t t1, mpz_t t2, const int n, mpz_t *v, mpz_t *moduli)
{
    mpz_set(t1, v[0]);
    mpz_set(t2, moduli[0]);
    for (int i = 1; i < n; i++)
    {
        CRT(x, t1, v[i], t2, moduli[i]);
        mpz_lcm(t2, t2, moduli[i]);
        
        mpz_set(t1, x);
    }
    mpz_mod(x, x, t2);
}


/**
 *
 * Computes message m from reduced ideal M = (f^2, L(m)f)
 * Source - https://eprint.iacr.org/2015/047.pdf
 *
 */
inline void CL::solve (mpz_t rop, QF& GE, const mpz_t con)
{
    // Get b part of M
    GE.getb(rop);
        
    // Get m from b by computing multiplicative inverse modulo con
    mpz_divexact(rop, rop, con);
    mpz_invert(rop, rop, con);
}


/**
 *
 * Computes message m from reduced ideal M = (f^2, L(m)f) using CRT1 method when f = (p1p2...pn)^t
 *
 */
inline void CL::solveCRT1 (mpz_t dm, mpz_t temp1, mpz_t temp2, QF& M, const mpz_t con, mpz_t *pntn, const int n, mpz_t *var)
{
    int temp;

    // Get L(m) part of M
    M.getb(temp1);
    mpz_divexact(temp1, temp1, con);
    
    // Compute L(m)^-1 mod pntn
    for (int i = 0; i < n; i++)
    {
        temp = mpz_invert(var[i], temp1, pntn[i]);
    }
    
    // Apply CRT to get m
    CRT_list(dm, temp1, temp2, n, var, pntn);
}


/**
 *
 * Computes message m from reduced ideal M = (f^2, L(m)f) using CRT2 method when f = (p1p2...pn)
 *
 */
inline void CL::solveCRT2 (mpz_t dm, mpz_t temp1, mpz_t temp2, QF& M, const mpz_t con, mpz_t *pn, const int n, mpz_t *var)
{
    int temp;

    // Get L(m) part of M
    M.getb(temp1);
    
    // Compute (X/pi)^-1 mod p1
    for (int i = 0; i < n; i++)
    {
        mpz_divexact(temp2, temp1, pn[i]);
        temp = mpz_invert(var[i], temp2, pn[i]);
    }
    
    // Apply CRT to get m
    CRT_list(dm, temp1, temp2, n, var, pn);
}


/**
 *
 * Computes message m when f = p^t.
 * Uses Pohlig-Hellman subroutine for groups of prime power order.
 * Pseudocode/Algorithm source - https://en.wikipedia.org/wiki/Pohlig%E2%80%93Hellman_algorithm
 *
 */
/*inline void CL::solve(mpz_t rop, QF& F, const QF& M, const mpz_t f, const mpz_t p, const int t)
{
    int var1;
    ZZ_limbs_set(zp, p[0]._mp_d, p[0]._mp_size);
    inv(T6, F);
    
    qfe2  .initialize (M, zp);
    qfed2 .initialize (T6, M, zp, zp);
    
    // Initialise x_0 to 0
    mpz_set_ui(rop, 0);
    
    // Compute p^(t-1)
    var1 = t - 1;
    e1 = power(zp, var1);
    
    
    // Compute F^(p^(t-1)) of order p
    qfe.power(T1, F, e1);
    
    //
    for(int k = 0; k < t; k++)
    {
        if (k == 0)
        {
            qfe2.power(T2, M, e1);
            T2.getb(temp3);
            mpz_divexact(temp4, temp3, p);
            mpz_invert(temp4, temp4, p);
        }
        
        else
        {
        
            // Compute (g^(-x_k)*h)^(p^(e-1-k))
            var1 = t - 1 - k;
            ZZ_limbs_set(zrop, rop[0]._mp_d, rop[0]._mp_size);
            e1 = power(zp, var1);
            e2 = zrop * e1;
            qfed2.power(T2, T6, M, e2, e1);
        
        
            // Compute d_k
            T2.getb(temp3);
            mpz_divexact(temp4, temp3, p);
            mpz_invert(temp4, temp4, p);
        }
        
        // Compute x_{k+1}
        mpz_pow_ui(temp1, p, k);
        mpz_mul(temp1, temp1, temp4);
        mpz_add(rop, rop, temp1);
    }
    mpz_mod(rop, rop, f);
}*/



inline void CL::solve(mpz_t rop, QF& GG, const QF& GE, const mpz_t con, const mpz_t prime, const int exp)
{
    // * Local temporaries
    int var1;
    
    // * Initialise
    ZZ_limbs_set(zp, prime[0]._mp_d, prime[0]._mp_size);
    inv(T6, GG);
    qfe.initialize(GG, zp);
    qfe2.initialize(GE, zp);
    qfed2.initialize(T6, GE, zp, zp);
    
    // STEP 1: Initialise x_0 to 0
    mpz_set_ui(rop, 0);
    
    // Compute p^(e-1)
    var1 = exp - 1;
    e1 = power(zp, var1);
    
    // STEP 2: Compute gamma = g^(p^(e-1))
    qfe.power(T1, GG, e1);
    
    // * Extract b from ideal T1 = (a,b,c)
    T1.getb(temp2);
    
    // Set variables for use when computing g^0; ANTL can't handle non-positive exponentiation
    GG.getD(temp1);
    mpz_set_ui(temp3, 1);
    
    // STEP 3:
    for(int k = 0; k < exp; k++)
    {
        // LOOP-STEP 1: Compute (g^(-x_k)*h)^(p^(e-1-k))
        var1 = exp - 1 - k;
        e1 = power(zp, var1);                              // Compute p^(e-1-k)
        ZZ_limbs_set(zrop, rop[0]._mp_d, rop[0]._mp_size); // Set (x_k)
        e2 = zrop * e1;                                    // Compute (x_k) (p^(e-1-k))
        if (zrop == 0)
        {
            qfe2.power(T2, GE, e1);
        }
        else
        {
            qfed2.power(T2, T6, GE, e2, e1);
        }
        T2.getb(temp3); // * Extract b from ideal T2 = (a,b,c)
        
        // LOOP Step 2: Compute d_k(temp4)
        if (mpz_cmp_ui(temp3, 1) != 0)
        {
            mpz_divexact(temp4, temp3, temp2);
            mpz_invert(temp4, temp4, prime);
        }
        else
        {
            mpz_set(temp4, prime);
        }
        mpz_pow_ui(temp1, prime, k);  // Compute p^k
        mpz_mul(temp1, temp1, temp4); // Compute (p^k)(d_k)
        
        // LOOP STEP 3: Compute x_k + 1
        mpz_add(rop, rop, temp1);
    }
    mpz_mod(rop, rop, con);
}





/**
 *
 * Computes message m from ideal M when f = (p1p2...pn)^t
 * Uses generic Pohlig Hellman subroutine
 * Pseudocode/Algorithm source - https://en.wikipedia.org/wiki/Pohlig%E2%80%93Hellman_algorithm
 *
 */
inline void CL::solve(mpz_t rop, QF& F, QF& M, const mpz_t con, const mpz_t* pntn, const mpz_t* pn, const int n, const int t)
{
    int temp;
    qfeM.initialize(M, e1);
    
    //Initial value
    mpz_set_ui(rop, 0);

    //
    for(int i = 0; i < n; i++)
    {
        mpz_divexact(stemp2, con, pntn[i]);
        ZZ_limbs_set(e1, stemp2[0]._mp_d, stemp2[0]._mp_size);
        
        qfeF.power(T3, F, e1);
        qfeM.power(T4, M, e1);
        solve(st, T3, T4, pntn[i], pn[i], t);
        mpz_set(var[i], st);
        
        // CRT steps
        temp = mpz_invert(st, stemp2, pntn[i]);
        mpz_mul(st, st, stemp2);
        mpz_mul(st, st, var[i]);
        mpz_add(rop, rop, st);
    }
    
    mpz_mod(rop, rop, con);
}

#endif // CL_HPP
