#ifndef BCPPUBLICKEY_HPP
#define BCPPUBLICKEY_HPP

class BCPPublickey
{
    protected:
    mpz_t ln;
    mpz_t lnn;
    mpz_t lg;
    mpz_t lh;
        
    public:
    BCPPublickey()  { mpz_inits(ln, lnn, lg, lh, NULL);  };
    ~BCPPublickey() { mpz_clears(ln, lnn, lg, lh, NULL); };
    void set(const mpz_t n, const mpz_t nn, const mpz_t g, const mpz_t h);
    void get(mpz_t n, mpz_t nn, mpz_t g, mpz_t h);
    void get(mpz_t n, mpz_t nn);
    void get(mpz_t n);
        
}; // end of class


inline void BCPPublickey::set(const mpz_t n, const mpz_t nn, const mpz_t g, const mpz_t h)
{
    mpz_set(ln,  n);
    mpz_set(lnn, nn);
    mpz_set(lg,  g);
    mpz_set(lh,  h);
}


inline void BCPPublickey::get(mpz_t n, mpz_t nn, mpz_t g, mpz_t h)
{
    mpz_set(n,  ln);
    mpz_set(nn, lnn);
    mpz_set(g,  lg);
    mpz_set(h,  lh);
}

inline void BCPPublickey::get(mpz_t n, mpz_t nn)
{
    mpz_set(n,  ln);
    mpz_set(nn, lnn);
}

inline void BCPPublickey::get(mpz_t n)
{
    mpz_set(n,  ln);
}

#endif // CLPUBLICKEY_HPP
