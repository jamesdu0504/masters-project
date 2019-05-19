#ifndef CLPUBLICKEY_HPP
#define CLPUBLICKEY_HPP

class CLPublickey
{
    protected:
    mpz_t p;
    mpz_t B;
    mpz_t UB;
    mpz_t con;
    mpz_t DK;
    mpz_t Dcon;
    mpz_t *pn;
    mpz_t *pntn;
    QF F;
    QF G;
    QF H;
    int n;
    int t;
        
        
    public:
    CLPublickey();
    ~CLPublickey();
    
    void init(mpz_qform_group_t& group);

    void set (const int N, const int T);
    void set (const mpz_t UB, const mpz_t con, const mpz_t DK, const mpz_t Dcon, const QF& F, const QF& G, const QF& H);
    void set (const mpz_t UB, const mpz_t con, const mpz_t DK, const mpz_t Dcon, const QF& F, const QF& G, const QF& H, mpz_t *pn, const int n);
    void set (const mpz_t UB, const mpz_t con, const mpz_t DK, const mpz_t Dcon, const QF& F, const QF& G, const QF& H, mpz_t *pn, const int n, const int t, mpz_t *pntn);
    void set (const mpz_t UB, const mpz_t con, const mpz_t DK, const mpz_t Dcon, const QF& F, const QF& G, const QF& H, const mpz_t p, const int t);
    

    void get (int N, int T);
    void get (mpz_t con);
    void get (mpz_t con, mpz_t Dcon);
    void get (mpz_t con, mpz_t DK, mpz_t Dcon);
    void get (mpz_t con, mpz_t Dcon, mpz_t* pn, int n);
    void get (mpz_t con, mpz_t DK, mpz_t Dcon, mpz_t* pn, int n);
    void get (mpz_t con, mpz_t *pn, int n);
    
    void get (mpz_t UB, mpz_t con, mpz_t Dcon, QF& G, QF& H);
    void get (mpz_t UB, mpz_t con, mpz_t Dcon, QF& G, QF& H, mpz_t *pn, int n);
    void get (mpz_t UB, mpz_t con, mpz_t DK, mpz_t Dcon, QF& F, QF& G, QF& H);
    void get (mpz_t UB, mpz_t con, mpz_t DK, mpz_t Dcon, QF& G, QF& H);
    void get (mpz_t UB, mpz_t con, QF& F, QF& G, QF& H); // p^t
    
    void get (QF& F, mpz_t con, mpz_t p, int t);
    void get (QF& F, mpz_t con, mpz_t* pn, int n, int t);
    void get (QF& F, mpz_t con, mpz_t DK, mpz_t Dcon, mpz_t p, int t);
    void get (QF& F, mpz_t con, mpz_t DK, mpz_t Dcon, mpz_t* pn, int n, int t);
        
        
}; // end of class

inline CLPublickey::CLPublickey()
{
    mpz_inits(p, B, UB, con, DK, Dcon, NULL);
}

inline CLPublickey::~CLPublickey()
{
    mpz_clears(p, B, UB, con, DK, Dcon, NULL);
}

// SET
inline void CLPublickey::set (const int N, const int T)
{
    n = N;
    t = T;
}

inline void CLPublickey::set (const mpz_t UB, const mpz_t con, const mpz_t DK, const mpz_t Dcon, const QF& F, const QF& G, const QF& H)
{
    mpz_set(this->UB,   UB);
    mpz_set(this->con,  con);
    mpz_set(this->DK,   DK);
    mpz_set(this->Dcon, Dcon);
    assign (this->F,    F);
    assign (this->G,    G);
    assign (this->H,    H);
}

inline void CLPublickey::init(mpz_qform_group_t& group)
{
    F.init(group);
    G.init(group);
    H.init(group);
}

inline void CLPublickey::set (const mpz_t UB, const mpz_t con, const mpz_t DK, const mpz_t Dcon, const QF& F, const QF& G, const QF& H, mpz_t *pn, const int n)
{
    mpz_set(this->UB,   UB);
    mpz_set(this->con,  con);
    mpz_set(this->DK,   DK);
    mpz_set(this->Dcon, Dcon);
    assign (this->F,    F);
    assign (this->G,    G);
    assign (this->H,    H);
    
    this->pn = pn;
    this->n = n;
}

inline void CLPublickey::set (const mpz_t UB, const mpz_t con, const mpz_t DK, const mpz_t Dcon, const QF& F, const QF& G, const QF& H, const mpz_t p, const int t)
{
    mpz_set(this->UB,   UB);
    mpz_set(this->con,  con);
    mpz_set(this->DK,   DK);
    mpz_set(this->Dcon, Dcon);
    assign (this->F,    F);
    assign (this->G,    G);
    assign (this->H,    H);
    mpz_set(this->p,    p);

    this->t = t;
}

inline void CLPublickey::set (const mpz_t UB, const mpz_t con, const mpz_t DK, const mpz_t Dcon, const QF& F, const QF& G, const QF& H, mpz_t *pn, const int n, const int t, mpz_t *pntn)
{
    mpz_set(this->UB,   UB);
    mpz_set(this->con,  con);
    mpz_set(this->DK,   DK);
    mpz_set(this->Dcon, Dcon);
    assign (this->F,    F);
    assign (this->G,    G);
    assign (this->H,    H);
    
    this->pn = pn;
    this->n = n;
    this->t = t;
    this->pntn = pntn;
}

// GET
inline void CLPublickey::get (int N, int T)
{
    N = n;
    T = t;
}

inline void CLPublickey::get (mpz_t con)
{
    mpz_set(con, this->con);
}

inline void CLPublickey::get (mpz_t con, mpz_t Dcon)
{
    mpz_set(con,  this->con);
    mpz_set(Dcon, this->Dcon);
}

inline void CLPublickey::get (mpz_t con, mpz_t DK, mpz_t Dcon)
{
    mpz_set(con,  this->con);
    mpz_set(DK,   this->DK);
    mpz_set(Dcon, this->Dcon);
}

inline void CLPublickey::get (mpz_t con, mpz_t *pn, int n)
{
    mpz_set(con, this->con);

    pn = this->pn;
    n = this->n;
}

inline void CLPublickey::get (mpz_t con, mpz_t Dcon, mpz_t* pn, int n)
{
    mpz_set(con,  this->con);
    mpz_set(Dcon, this->Dcon);
    pn = this->pn;
    n = this->n;
}

inline void CLPublickey::get (mpz_t UB, mpz_t con, mpz_t Dcon, QF& G, QF& H)
{
    mpz_set(UB,   this->UB);
    mpz_set(con,  this->con);
    mpz_set(Dcon, this->Dcon);
    assign (G,    this->G);
    assign (H,    this->H);
}

inline void CLPublickey::get (mpz_t UB, mpz_t con, mpz_t Dcon, QF& G, QF& H, mpz_t *pn, int n)
{
    mpz_set(UB,   this->UB);
    mpz_set(con,  this->con);
    mpz_set(Dcon, this->Dcon);
    assign (G,    this->G);
    assign (H,    this->H);

    pn = this->pn;
    n = this->n;
}

inline void CLPublickey::get(mpz_t UB, mpz_t con, QF& F, QF& G, QF& H)
{
    mpz_set(UB,  this->UB);
    mpz_set(con, this->con);
    assign (F,   this->F);
    assign (G,   this->G);
    assign (H,   this->H);
}

inline void CLPublickey::get(QF& F, mpz_t con, mpz_t p, int t)
{
    assign (F,   this->F);
    mpz_set(con, this->con);
    mpz_set(p,   this->p);

    t = this->t;
}

inline void CLPublickey::get (QF& F, mpz_t con, mpz_t* pn, int n, int t)
{
    assign (F,   this->F);
    mpz_set(con, this->con);

    pn = this->pn;
    n = this->n;
    t = this->t;
}

inline void CLPublickey::get (mpz_t UB, mpz_t con, mpz_t DK, mpz_t Dcon, QF& G, QF& H)
{
    mpz_set(UB,   this->UB);
    mpz_set(con,  this->con);
    mpz_set(DK,   this->DK);
    mpz_set(Dcon, this->Dcon);
    assign (G,    this->G);
    assign (H,    this->H);
}

inline void CLPublickey::get (mpz_t UB, mpz_t con, mpz_t DK, mpz_t Dcon, QF& F, QF& G, QF& H)
{
    mpz_set(UB,   this->UB);
    mpz_set(con,  this->con);
    mpz_set(DK,   this->DK);
    mpz_set(Dcon, this->Dcon);
    assign (F,    this->F);
    assign (G,    this->G);
    assign (H,    this->H);
}

inline void CLPublickey::get (mpz_t con, mpz_t DK, mpz_t Dcon, mpz_t* pn, int n)
{
    mpz_set(con,  this->con);
    mpz_set(DK,   this->DK);
    mpz_set(Dcon, this->Dcon);
    
    pn = this->pn;
    n = this->n;
}

inline void CLPublickey::get (QF& F, mpz_t con, mpz_t DK, mpz_t Dcon, mpz_t p, int t)
{
    assign (F,    this->F);
    mpz_set(con,  this->con);
    mpz_set(DK,   this->DK);
    mpz_set(Dcon, this->Dcon);
    mpz_set(p,    this->p);
    
    t = this->t;
}

inline void CLPublickey::get (QF& F, mpz_t con, mpz_t DK, mpz_t Dcon, mpz_t* pn, int n, int t)
{
    assign (F,    this->F);
    mpz_set(con,  this->con);
    mpz_set(DK,   this->DK);
    mpz_set(Dcon, this->Dcon);
    
    pn = this->pn;
    n = this->n;
    t = this->t;
}

#endif // CLPUBLICKEY_HPP
