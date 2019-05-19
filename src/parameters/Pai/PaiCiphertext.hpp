#ifndef PaiCIPHERTEXT_HPP
#define PaiCIPHERTEXT_HPP

class PaiCiphertext
{
    protected:
        mpz_t lc;
        
    public:
        PaiCiphertext()  { mpz_init(lc);  };
        ~PaiCiphertext() { mpz_clear(lc); };
        
        void get(mpz_t c);
        void set(const mpz_t c);
}; // end of class

inline void PaiCiphertext::get(mpz_t c)
{
    mpz_set(c, lc);
}

inline void PaiCiphertext::set(const mpz_t c)
{
    mpz_set(lc, c);
}

#endif // PaiCIPHERTEXT_HPP
