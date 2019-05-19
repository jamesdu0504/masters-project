#ifndef CLSECRETKEY_HPP
#define CLSECRETKEY_HPP

class CLSecretkey
{
    protected:
    ZZ x;
        
    public:
    CLSecretkey() {};
    ~CLSecretkey() {};
        
    void get(ZZ& x)       { x = this->x; }
    void set(const ZZ& x) { this->x = x; }
    
}; // end of class

#endif // CLSECRETKEY_HPP
