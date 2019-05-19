#ifndef CLPLAINTEXT_HPP
#define CLPLAINTEXT_HPP

class CLPlaintext
{
    protected:
    mpz_t mm;
    ZZ zm;
    
    public:
    CLPlaintext();
    ~CLPlaintext();
    
    void get(mpz_t m);
    void get(ZZ &m);
    void get(mpz_t m1, ZZ &zm1);
    void set(const mpz_t m);
    
}; // end of class

inline CLPlaintext::CLPlaintext()
{
    mpz_init(mm);
}

inline CLPlaintext::~CLPlaintext()
{
    mpz_clear(mm);
}

inline void CLPlaintext::get(mpz_t m)
{
    mpz_set(m, mm);
}

inline void CLPlaintext::get(ZZ &m)
{
    m = zm;
}

inline void CLPlaintext::get(mpz_t m1, ZZ &zm1)
{
    mpz_set(m1, mm);
    zm1 = zm;
}

inline void CLPlaintext::set(const mpz_t m)
{
    mpz_set(this->mm, m);
    ZZ_limbs_set(zm, m[0]._mp_d, m[0]._mp_size);
}

#endif // CLPLAINTEXT_HPP
