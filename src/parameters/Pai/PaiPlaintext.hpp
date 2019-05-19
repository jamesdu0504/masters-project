#ifndef PaiPLAINTEXT_HPP
#define PaiPLAINTEXT_HPP

class PaiPlaintext
{
    protected:
        mpz_t lm;
    public:
        PaiPlaintext()  { mpz_init(lm); };
        ~PaiPlaintext() { mpz_clear(lm); };
        
        void get(mpz_t m) { mpz_set(m, lm); }
        void set(const mpz_t m) { mpz_set(lm, m); }
        
    
}; // end of class


#endif // PaiPLAINTEXT_HPP
