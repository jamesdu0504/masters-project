#ifndef BCPCIPHERTEXT_HPP
#define BCPCIPHERTEXT_HPP

class BCPCiphertext
{
    protected:
    mpz_t lc1;
    mpz_t lc2;
    
    public:
    BCPCiphertext()  { mpz_inits(lc1, lc2, NULL);  };
    ~BCPCiphertext() { mpz_clears(lc1, lc2, NULL); };
    void get(mpz_t c1, mpz_t c2);
    void set(const mpz_t c1, const mpz_t c2);
    
}; // end of class

inline void BCPCiphertext::get(mpz_t c1, mpz_t c2)
{
    mpz_set(c1, lc1);
    mpz_set(c2, lc2);
}

inline void BCPCiphertext::set(const mpz_t c1, const mpz_t c2)
{
    mpz_set(lc1, c1);
    mpz_set(lc2, c2);
}

#endif // BCPCIPHERTEXT_HPP
