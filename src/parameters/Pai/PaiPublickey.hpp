#ifndef PaiPUBLICKEY_HPP
#define PaiPUBLICKEY_HPP

class PaiPublickey
{
    protected:
    mpz_t ln;
    mpz_t lnn;
    mpz_t lg;
        
    public:
    PaiPublickey()  { mpz_inits(ln, lnn, lg, NULL);  };
    ~PaiPublickey() { mpz_clears(ln, lnn, lg, NULL); };
        
    void set(const mpz_t n, const mpz_t nn, const mpz_t g);
    void get(mpz_t n, mpz_t nn);
    void get(mpz_t n, mpz_t nn, mpz_t g);
    void get(mpz_t n);    
        
}; // end of class


inline void PaiPublickey::set(const mpz_t n, const mpz_t nn, const mpz_t g)
{
    mpz_set(ln,  n);
    mpz_set(lnn, nn);
    mpz_set(lg,  g);
}

inline void PaiPublickey::get(mpz_t n, mpz_t nn)
{
    mpz_set(n,  ln);
    mpz_set(nn, lnn);
}

inline void PaiPublickey::get(mpz_t n, mpz_t nn, mpz_t g)
{
    mpz_set(n,  ln);
    mpz_set(nn, lnn);
    mpz_set(g,  lg);
}

inline void PaiPublickey::get(mpz_t n)
{
    mpz_set(n,  ln);
}

#endif // PaiPUBLICKEY_HPP
