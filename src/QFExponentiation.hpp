/**
 *
 * AUTHOR - Parthasarathi Das
 * BRIEF  - A class for ideal exponentiation using ideal objects
 */

#ifndef QFEXPONENTIATION_HPP
#define QFEXPONENTIATION_HPP

#include <gmp.h>
#include <NTL/ZZ.h>
#include "QF.hpp"
#include "../antl/include/ANTL/Exponentiation/ExponentiationNAF.hpp"
#include "../antl/include/ANTL/Exponentiation/ExponentiationWNAF.hpp"
#include "../antl/include/ANTL/Exponentiation/DoubleExponentiation/DoubleExponentiationIL.hpp"

NTL_CLIENT

class QFExponentiation
{
    public:
    QFExponentiation();
    ~QFExponentiation();
    
    void init(const int w);
    void initialize (const QF& A, const ZZ& n);
    void initialize (const QF& A, const QF& B, const ZZ& m, const ZZ& n);
    void power      (QF& C, const QF& A, const ZZ& n);
    void power      (QF& C, const QF& A, const QF& B, const ZZ& m, const ZZ& n);
    void CLpower    (QF& C, const mpz_t con, const mpz_t D, const mpz_t m);
    
    protected:
    int                        window_size;
    mpz_t                      temp1, temp2;
    ExponentiationWNAF<QF>     wnaf;
    DoubleExponentiationIL<QF> il;
        
};

inline QFExponentiation::QFExponentiation()
{
    mpz_inits(temp1, temp2, NULL);
}

inline QFExponentiation::~QFExponentiation()
{
    mpz_clears(temp1, temp2, NULL);
}

inline void QFExponentiation::init(const int w)
{
    window_size = w;
}

inline void QFExponentiation::initialize(const QF& A, const ZZ& n)
{
    wnaf.initialize (A, n, window_size);
}

inline void QFExponentiation::initialize(const QF& A, const QF& B, const ZZ& m, const ZZ& n)
{
    il.initialize (A, B, m, n, window_size, window_size);
}

inline void QFExponentiation::power (QF& C, const QF& A, const ZZ& n)
{
    wnaf.power (C, A, n);
}

inline void QFExponentiation::power (QF &C, const QF &A, const QF &B, const ZZ & m, const ZZ & n)
{
    il.power (C, A, B, m, n);
}

inline void QFExponentiation::CLpower(QF& C, const mpz_t con, const mpz_t D, const mpz_t m)
{
    mpz_invert(temp1, m, con);
    if (mpz_mod_ui(temp2, temp1, 2) == 0)
    {
        mpz_sub(temp1, temp1, con);
    }
    mpz_mul(temp1, temp1, con);
    mpz_mul(temp2, con, con);
    C.assign(temp2, temp1, D);
}

#endif // QFEXPONENTIATION_HPP
