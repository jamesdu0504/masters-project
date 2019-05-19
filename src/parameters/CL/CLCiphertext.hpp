#ifndef CLCIPHERTEXT_HPP
#define CLCIPHERTEXT_HPP

#include <vector>

class CLCiphertext
{
    public:
    CLCiphertext();
    ~CLCiphertext();

    void init(mpz_qform_group_t& group);
    void set(const QF& C1, const QF& C2);
    void set(const QF& C1, const vector<QF> &cn, const int n);
    void get(QF& C1, QF& C2);
    void get(QF& C1, vector<QF> &cn, const int n);
    
    friend std::ostream& operator<< (std::ostream& out, const CLCiphertext& clct);
    
    
    protected:
    QF C1;
    QF C2;
    vector<QF> CN;
    int len;
    
}; // end of class

inline CLCiphertext::CLCiphertext()
{
    len = 8;    
}

inline CLCiphertext::~CLCiphertext()
{
    
}

inline void CLCiphertext::init(mpz_qform_group_t& group)
{
    C1.init(group);
    C2.init(group);
    CN.resize(len);
    for (int i = 0; i < len; i++)
    {
        CN[i].init(group);
    }
    
}

// SET
inline void CLCiphertext::set(const QF& C1, const QF& C2)
{
    assign(this->C1, C1);
    assign(this->C2, C2);
}

inline void CLCiphertext::set(const QF& C1, const vector<QF> &cn, const int n)
{
    assign(this->C1, C1);
    for (int i = 0; i < n; i++)
    {
        assign(CN[i], cn[i]);
    }
    
}

// GET
inline void CLCiphertext::get(QF& C1, QF& C2)
{
    assign(C1, this->C1);
    assign(C2, this->C2);
}

inline void CLCiphertext::get(QF& C1, vector<QF> &cn, const int n)
{
    assign(C1, this->C1);
    for (int i = 0; i < n; i++)
    {
        assign(cn[i], CN[i]);
    }
    
}

inline std::ostream& operator<< (std::ostream& out, const CLCiphertext& clct)
{
    cout << "C1 = " << clct.C1;
    cout << "C2 = " << clct.C2;
    return out;
}

#endif // CLCIPHERTEXT_HPP
