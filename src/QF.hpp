/**
 *
 * AUTHOR - Parthasarathi Das
 * BRIEF  - A class for quadratic ideals
 */

#ifndef QF_HPP
#define QF_HPP

#include <NTL/ZZ_limbs.h>

NTL_CLIENT

class QF
{
    public:
    QF();
    QF(mpz_qform_group_t& group);
    ~QF();
    
    void init            (mpz_qform_group_t& group);
    void assign          (const mpz_t a, const mpz_t b, const mpz_t D);
    void getb            (mpz_t b);
    void setD            (const mpz_t D);
    void getD            (mpz_t D);
    void random_prime_QF ();
    int  next_prime_QF   (int& index, const mpz_t D);
    void findQFprimeto   (const QF& I, mpz_t temp, const mpz_t D, const mpz_t con);
    void clear           (mpz_qform_group_t& group);
    
    friend void assign             (QF& Z, const QF& X);
    friend void clear              (QF& Z);
    friend bool equals             (QF& Z, const QF& X);
    friend void inv                (QF& Z, const QF& X);
    friend void mul                (QF& Z, const QF& X, const QF& Y); // forms composition
    friend void sqr                (QF& Z, const QF& X);
    friend void cub                (QF& Z, const QF& X);
    friend void sendQFtosmallerord (QF& Z, mpz_t temp, const QF& I, const mpz_t DK, const mpz_t D, const mpz_t con);
    
    friend std::ostream& operator<< (std::ostream& out, const QF& A);
    friend bool          operator== (const QF& X, const QF& Y);
    
    
    protected:
    mpz_qform_group_t* gp;
    mpz_qform_t form;
};

// Out of line class member definitions
inline QF::QF ()
{
    mpz_qform_init (gp, &form);
}

inline QF::~QF ()
{
    mpz_qform_clear (gp, &form);
}

inline QF::QF (mpz_qform_group_t& group)
{
    gp = &group;
    mpz_qform_init (gp, &form);
}

inline void QF::init (mpz_qform_group_t& group)
{
    gp = &group;
}

inline void QF::clear (mpz_qform_group_t& group)
{
    mpz_qform_group_clear (&group);
}

inline void QF::assign (const mpz_t a, const mpz_t b, const mpz_t D)
{
    mpz_qform_group_set_discriminant (gp, D);
    mpz_set (form.a, a);
    mpz_set (form.b, b);
    mpz_qform_c (gp, form.c, form.a, form.b);
}

inline void QF::getb (mpz_t b)
{
    mpz_set (b, form.b);
}
inline void QF::setD (const mpz_t D)
{
    mpz_qform_group_set_discriminant (gp, D);
}

inline void QF::getD (mpz_t D)
{
    mpz_set (D, gp->D);
}

inline void QF::random_prime_QF ()
{
    qform_random_primeform (&gp->desc, &form);
}

inline int QF::next_prime_QF (int& index, const mpz_t D)
{
    mpz_qform_group_set_discriminant (gp, D);
    return qform_next_primeform (&gp->desc, &form, index);
}

inline void QF::findQFprimeto(const QF& I, mpz_t temp, const mpz_t D, const mpz_t con)
{    
    mpz_qform_group_set_discriminant (gp, D);
    mpz_gcd (temp, I.form.a, con);
    
    if(mpz_cmp_ui (temp, 1) > 0)
    {
        mpz_gcd (temp, I.form.c, con);
        if(mpz_cmp_ui(temp, 1) > 0)
        {
            mpz_add (form.a, I.form.a, I.form.b);
            mpz_add (form.a, form.a, I.form.c);
            mpz_mul_ui (temp, I.form.a, 2);
            mpz_add (form.b, I.form.b, temp);
            mpz_neg (form.b, form.b);
        }
        else
        {
            mpz_set (form.a, I.form.c);
            mpz_set (form.b, I.form.b);
            mpz_neg (form.b, form.b);
        }
        
        mpz_qform_c (gp, form.c, form.a, form.b);
    }
    else
    {
        mpz_set (form.a,I.form.a);
        mpz_set (form.b,I.form.b);
        mpz_qform_c (gp, form.c, form.a, form.b);
    }
}

// Non class members
inline void assign (QF& Z, const QF& X)
{
    Z.gp = X.gp;
    mpz_qform_set (Z.gp, &Z.form, &X.form);
}

inline void clear (QF& Z)
{
    mpz_qform_clear (Z.gp, &Z.form);
}

inline bool equals (QF& Z, const QF& X)
{
    return mpz_qform_equal (Z.gp, &Z.form, &X.form);
}

inline void inv (QF& Z, const QF& X)
{
    assign(Z, X);
    mpz_qform_inverse (Z.gp, &Z.form);
}

inline void mul (QF& Z, const QF& X, const QF& Y)
{
    Z.gp = X.gp;
    mpz_qform_compose (Z.gp, &Z.form, &X.form, &Y.form);
}

inline void sqr (QF& Z, const QF& X)
{
    Z.gp = X.gp;
    mpz_qform_square (Z.gp, &Z.form, &X.form);
}

inline void cub (QF& Z, const QF& X)
{
    Z.gp = X.gp;
    mpz_qform_cube (Z.gp, &Z.form, &X.form);
}

inline void sendQFtosmallerord(QF& Z, mpz_t temp, const QF& I, const mpz_t Dsmall, const mpz_t Dlarge, const mpz_t con)
{    
    Z.findQFprimeto (I, temp, Dsmall, con); // Here group has D = Dsmall;
    mpz_mul (Z.form.b, Z.form.b, con);
    mpz_mul_2exp (temp, Z.form.a, 1);
    mpz_mod (Z.form.b, Z.form.b, temp);
    mpz_qform_group_set_discriminant (Z.gp, Dlarge); // Now group has D = Dlarge
    mpz_qform_c (Z.gp, Z.form.c, Z.form.a, Z.form.b);
}

inline std::ostream& operator<< (std::ostream& out, const QF& A)
{
    mpz_qform_print (A.gp, &A.form);
    return out;
}

inline bool operator== (const QF& X, const QF& Y)
{
    return mpz_qform_equal (X.gp, &X.form, &Y.form);
}

#endif // QF_HPP
