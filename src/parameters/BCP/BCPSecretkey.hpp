#ifndef BCPSECRETKEY_HPP
#define BCPSECRETKEY_HPP

class BCPSecretkey
{
    protected:
    mpz_t la;
        
    public:
    BCPSecretkey() { mpz_init(la); };
    ~BCPSecretkey() {mpz_clear(la); };
    void set(const mpz_t a) { mpz_set(la, a); }
    void get(mpz_t a)       { mpz_set(a, la); }
    
}; // end of class

#endif // BCPSECRETKEY_HPP
