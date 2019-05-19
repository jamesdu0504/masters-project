#ifndef BCPPLAINTEXT_HPP
#define BCPPLAINTEXT_HPP

class BCPPlaintext
{
    protected:
    mpz_t lm;
    
    public:
    BCPPlaintext()  { mpz_init(lm);  };
    ~BCPPlaintext() { mpz_clear(lm); };
    void get(mpz_t m) { mpz_set(m, lm); }
    void set(const mpz_t m) { mpz_set(lm, m); }
        
    
}; // end of class

#endif // BCPPLAINTEXT_HPP
