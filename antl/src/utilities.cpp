/**
 * @file utilities.cpp
 * @author Michael Jacobson
 * @remark specialized implementations of utility routines.
 */

#include <ANTL/utilities.hpp>

//
// Polynomial Utility Methods
//

// Convert base-q integer representation to polynomial
void
get_poly_modq (GF2X & p, const ZZ & X, const ZZ & q)
{
  GF2 coeff;
  ZZ a, b, temp;
  long i = 0;

  clear (p);
  a = X;
  while (!IsZero (a))
    {
      DivRem (temp, b, a, q);
      a = temp;

      conv (coeff, b);

      SetCoeff (p, i, coeff);
      ++i;
    }
}

void
get_poly_modq (GF2EX & p, const ZZ & X, const ZZ & q)
{
  GF2E pp, one;
  GF2E xx;
  ZZ a, b, temp;
  long i = 0;

  GF2X x2;
  SetX(x2);
  conv(xx,x2);

  clear (p);
  a = X;
  while (!IsZero (a))
    {
      DivRem (temp, b, a, q);
      a = temp;

      set (one);
      clear (pp);
      while (b > 0)
	{
	  if (IsOdd (b))
	    pp += one;

	  one *= xx;
	  b >>= 1;
	}

      SetCoeff (p, i, pp);
      ++i;
    }
}

void
get_poly_modq (ZZ_pX & p, const ZZ & X, const ZZ & q)
{
  ZZ_p coeff;
  ZZ a, b, temp;
  long i = 0;

  clear (p);
  a = X;
  while (!IsZero (a))
    {
      DivRem (temp, b, a, q);
      a = temp;

      conv (coeff, b);

      SetCoeff (p, i, coeff);
      ++i;
    }
}



void
get_poly_modq (zz_pX & p, const ZZ & X, const ZZ & q)
{
  zz_p coeff;
  ZZ a, b, temp;
  long i = 0;

  clear (p);
  a = X;
  while (!IsZero (a))
    {
      DivRem (temp, b, a, q);
      a = temp;

      conv (coeff, b);

      SetCoeff (p, i, coeff);
      ++i;
    }
}

void
get_poly_modq (ZZ_pEX & p, const ZZ & X, const ZZ & q)
{
  ZZ_pX repcoeff;
  ZZ pp = ZZ_p::modulus ();
  ZZ_pE coeff;
  ZZ a, b, temp;
  long i = 0;

  clear (p);
  a = X;
  while (!IsZero (a))
    {
      DivRem (temp, b, a, q);
      a = temp;

      get_poly_modq (repcoeff, b, pp);
      conv (coeff, repcoeff);

      SetCoeff (p, i, coeff);
      ++i;
    }
}

void
get_poly_modq (zz_pEX & p, const ZZ & X, const ZZ & q)
{
  zz_pX repcoeff;
  ZZ pp = to_ZZ (zz_p::modulus ());
  zz_pE coeff;
  ZZ a, b, temp;
  long i = 0;

  clear (p);
  a = X;
  while (!IsZero (a))
    {
      DivRem (temp, b, a, q);
      a = temp;

      get_poly_modq (repcoeff, b, pp);
      conv (coeff, repcoeff);

      SetCoeff (p, i, coeff);
      ++i;
    }
}



//
// eval_poly - evaluate polynomial at q
//

template <>
ZZ eval_poly <ZZ> (const ZZ & p, const ZZ & q)
{
  return p;
}



template <>
ZZ eval_poly <ZZ_pEX> (const ZZ_pEX & p, const ZZ & q)
{
  ZZ hval;
  static vec_ZZ qvec;
  static ZZ pp;
  long i,j,da;

  // MV - Added the || ... because what if q changes? (It was an issue)
  if (qvec.length() == 0 || (qvec.length() > 1 && qvec[1] != q)) {
    qvec.SetLength(1);
    qvec[0] = 1;
  }
  // make sure static array of powers of q has enough elements!
  if (qvec.length() == 0) {
    pp= ZZ_p::modulus();
    qvec.SetLength(1);
    qvec[0] = 1;
  }

  if (qvec.length() <= deg(p)) {
    j = qvec.length();
    qvec.SetLength(deg(p)+1);

    for (i=j; i<=deg(p); ++i)
      qvec[i] = qvec[i-1]*q;
  }
  pp = ZZ_p::modulus(); //Moved here because it wasn't quite working.
  clear(hval);
  da = deg (p);
  // cout << "ZZ_pEX qvec = " << qvec << endl;
  // cout << "pp = " << pp << endl;
  for (i=0; i<=da; ++i)
    hval += eval_poly(rep(coeff(p,i)),pp) * qvec[i];

  return hval;
}


template <>
ZZ eval_poly <zz_pEX> (const zz_pEX & p, const ZZ & q)
{
  ZZ hval;
  static vec_ZZ qvec;
  static ZZ pp;
  long i,j,da;


  // make sure static array of powers of q has enough elements!
  if (qvec.length() == 0) {
    pp= to_ZZ(zz_p::modulus());
    qvec.SetLength(1);
    qvec[0] = 1;
  }

  if (qvec.length() <= deg(p)) {
    j = qvec.length();
    qvec.SetLength(deg(p)+1);

    for (i=j; i<=deg(p); ++i)
      qvec[i] = qvec[i-1]*q;
  }

  clear(hval);
  da = deg (p);

  for (i=0; i<=da; ++i)
    hval += eval_poly(rep(coeff(p,i)),pp) * qvec[i];

  return hval;
}


template <>
ZZ eval_poly <GF2EX> (const GF2EX & p, const ZZ & q)
{
  GF2X crep;
  ZZ X, qx, cval, x2;
  static vec_ZZ qvec;
  register long i, j;

  if (qvec.length() == 0) {
    qvec.SetLength(1);
    qvec[0] = 1;
  }

  // make sure static array of powers of q has enough elements!
  if (qvec.length() <= deg(p)) {
    j = qvec.length();
    qvec.SetLength(deg(p)+1);

    for (i=j; i<=deg(p); ++i)
      qvec[i] = qvec[i-1]*q;
  }

  // compute X
  X = 0;
  qx = 1;
  for (i=0; i<=deg(p); ++i) {
    crep = rep(coeff(p,i));
    if (!IsZero(crep)) {
      // compute cval
      cval = 0;
      x2 = 1;
      for (j=0; j<=deg(crep); ++j) {
        if (!IsZero(coeff(crep,j)))
          cval += x2;
        x2 <<= 1;
      }

      X += cval*qvec[i];
    }
  }

  return X;
}




//
// Quadratic residuosity functions
//

/*
 * Function: Jacobi_base
 * Purpose: tests to see if a is a square mod p
 */

long
Jacobi_base(const ZZ & a, const ZZ & n)
{
  ZZ aa, nn;
  long t, k, d;
  
  nn = n;
  aa = a;
  t = 1;

  while (aa != 0) {
    k = MakeOdd(aa);
    d = trunc_long(nn, 3);
    if ((k & 1) && (d == 3 || d == 5)) t = -t;

    if (trunc_long(aa, 2) == 3 && (d & 3) == 3) t = -t;
    swap(aa, nn);
    rem(aa, aa, nn);
  }

  if (nn == 1)
    return t;
  else
    return 0;
}


long
Jacobi_base (const long & a, const long & n)
{
  long t,k,d;
  long nn,aa,temp;

  nn = n;
  aa = a;

  t = 1;

  while (aa != 0) {
    k = 0;
    while (!(aa & 1)) {
      aa >>= 1;
      ++k;
    }

    d = (nn & 7);
    if ((k & 1) && (d == 3 || d == 5)) t = -t;

    if ((aa & 3) == 3 && (d & 3) == 3) t = -t;
    temp = aa;
    aa = nn;
    nn = temp;
    aa %= nn;
  }

  if (nn == 1)
    return t;
  else
    return 0;
}



long
Jacobi(const ZZ & a, const ZZ & n)
{
  ZZ temp;
  rem(temp,a,n);
  return Jacobi_base(temp,n);
}



long
Jacobi(const long & a, const long & n)
{
  long temp = a % n;
  if (temp < 0)  temp += n;
  return Jacobi_base(temp,n);
}



long
Jacobi (const ZZ_pX & a, const ZZ_pX & n)
{
  ZZ_pX temp;
  rem(temp,a,n);
  if (IsZero(temp))  return 0;

  ZZ_pX P1, P2;
  ZZ c,q;
  q = ZZ_p::modulus ();

  return Jacobi_base (a, n, P1, P2, c, q);
}

long
Jacobi (const zz_pX & a, const zz_pX & n)
{
  zz_pX temp;
  rem(temp,a,n);
  if (IsZero(temp))  return 0;

  zz_pX P1, P2;
  long c, q;
  q = (long) zz_p::modulus ();

  return Jacobi_base (a, n, P1, P2, c, q);
}

long
Jacobi (const ZZ_pEX & a, const ZZ_pEX & n)
{
  ZZ_pEX temp;
  rem(temp,a,n);
  if (IsZero(temp))  return 0;

  ZZ_pEX P1, P2;
  ZZ_pX c;
  ZZ_pX q = ZZ_pE::modulus ();

  return Jacobi_base (temp, n, P1, P2, c, q);
}

long
Jacobi (const zz_pEX & a, const zz_pEX & n)
{
  zz_pEX temp;
  rem(temp,a,n);
  if (IsZero(temp))  return 0;

  zz_pEX P1, P2;
  zz_pX c;
  zz_pX q = zz_pE::modulus ();

  return Jacobi_base (temp, n, P1, P2, c, q);
}


//
// modular square root functions
//

// get_qdp - computes q^deg(p), where q is the finite field size

void get_qdp(ZZ & q, const ZZ_pX & p)
{
  power(q,ZZ_p::modulus(),deg(p));
}

void get_qdp(ZZ & q, const zz_pX & p)
{
  power(q,zz_p::modulus(),deg(p));
}

void get_qdp(ZZ & q, const ZZ_pEX & p)
{
  power(q,ZZ_pE::cardinality(),deg(p));
}

void get_qdp(ZZ & q, const zz_pEX & p)
{
  power(q,zz_pE::cardinality(),deg(p));
}




//
// GF2EX versions of quadratic residuosity and modular square root funcionts
//

/*
 * Function: trace
 * Purpose: Compute the trace of 'a' (mod p), i.e.
 *          trace = sum from i=0 to v-1 of a^2^i (mod p)
 *          where deg(p) = v;
 * Inputs: GF2EX & t - a pointer to an instance of GF2EX which
 *                     contains the value of trace at the end of the function.
 *         GF2EX & a - The instance of GF2EX to find the trace of.
 *         GF2EXModulus & pmod -  The modulus under which the
 *                                calculation is preformed.
 * Outputs: NONE
 */

inline void
trace (GF2X & t, const GF2X & a, const GF2X & pmod)
{
  GF2X temp, temp2;
  long i, bound;

  bound = deg (pmod);

  clear (t);
  rem(temp,a,pmod);
  for (i = 0; i < bound; ++i)
    {
      t += temp;

      SqrMod (temp2, temp, pmod);
      temp = temp2;
    }
}

/*
 * Function: Jacobi(h,f,n)
 * Purpose: Compute the Artin symbol [D/P]. Determines the splitting behaviour
 *          of p.
 * Inputs: const GF2EX & h - polynomial of at most degree = genus
 *                           same as the hx in a hyperelliptic curve.
 *         const GF2EX & f - polynomial of degree 2g+1, again same as f
 *                           in a hyperelliptic curve.
 *         const GF2EX & p - polynomial to test.
 * Outputs: long - 0  if p | h
 *                 -1 if y^2 + hy - f = 0 (mod p) has no solutions
 *                 1  if y^2 + hy - f = 0 (mod p) has 2 solutions
 */

long
Jacobi (const GF2X & h, const GF2X & f, const GF2X & n)
{
  GF2X hn;
  rem(hn,h,n);
  if (IsZero (hn)) return 0;

  // Compute a = f h^-2.  Solubility of y^2 + y + a = 0 (mod p) is
  // equivalent to solubility of y^2 + hy -f = 0 (mod p)

  GF2X temp, tr, a;

  InvMod (temp, hn, n);
  sqr(a,temp);
  MulMod(a,a,f,n);

  // sol exists if tr(a) = sum_{i=0}^{vt-1} a^2^i = 0 (mod p)
  trace (tr, a, n);

  if (IsZero (tr))
    return 1;
  else
    return -1;
}



/*
 * Function: trace
 * Purpose: Compute the trace of 'a' (mod p), i.e.
 *          trace = sum from i=0 to v-1 of a^2^i (mod p)
 *          where deg(p) = v;
 * Inputs: GF2EX & t - a pointer to an instance of GF2EX which
 *                     contains the value of trace at the end of the function.
 *         GF2EX & a - The instance of GF2EX to find the trace of.
 *         GF2EXModulus & pmod -  The modulus under which the
 *                                calculation is preformed.
 * Outputs: NONE
 */

inline void
trace (GF2EX & t, const GF2EX & a, const GF2EXModulus & pmod)
{
  GF2EX temp, temp2;
  long i, bound;

  bound = deg (pmod) * GF2E::degree ();

  clear (t);
  rem(temp,a,pmod);
  for (i = 0; i < bound; ++i)
    {
      t += temp;

      SqrMod (temp2, temp, pmod);
      temp = temp2;
    }
}

/*
 * Function: Jacobi(h,f,n)
 * Purpose: Compute the Artin symbol [D/P]. Determines the splitting behaviour
 *          of p.
 * Inputs: const GF2EX & h - polynomial of at most degree = genus
 *                           same as the hx in a hyperelliptic curve.
 *         const GF2EX & f - polynomial of degree 2g+1, again same as f
 *                           in a hyperelliptic curve.
 *         const GF2EX & p - polynomial to test.
 * Outputs: long - 0  if p | h
 *                 -1 if y^2 + hy - f = 0 (mod p) has no solutions
 *                 1  if y^2 + hy - f = 0 (mod p) has 2 solutions
 */

long
Jacobi (const GF2EX & h, const GF2EX & f, const GF2EX & n)
{
  GF2EX hn, fn, an;
  rem(hn,h,n);
  if (IsZero (hn)) return 0;

  // Compute a = f h^-2.  Solubility of y^2 + y + a = 0 (mod p) is
  // equivalent to solubility of y^2 + hy -f = 0 (mod p)

  GF2EX temp, tr, a;

  InvMod (temp, hn, n);
  sqr(a,temp);
  rem(fn, f,n); 	//FIXME MDV When eventually called by lower_bound_hR,
  rem(an, a, n);	//MulMod complains about bad args, the degree of
			// f, a was too big.
  MulMod(a,an,fn,n);

  // sol exists if tr(a) = sum_{i=0}^{vt-1} a^2^i = 0 (mod p)
  trace (tr, a, n);

  if (IsZero (tr))
    return 1;
  else
    return -1;
}




//
// char 2 ressol function
//

long ressol (GF2EX & x, const GF2EX & h, const GF2EX & f, const GF2EX & p)
{
  GF2EX temp, temp2, tr, a, c, d, atemp, ctemp;
  ZZ e;
  long i, bound;

  bound = deg (p) * GF2E::degree ();

  if (IsZero (h % p))
    {
      e = 1 << (bound - 1);
      PowerMod (x, f % p, e, p);

      return 0;
    }

  // Compute a = f h^-2.  Find a solution s of y^2 + y + a = 0 (mod p)

  InvMod (temp, h % p, p);
  a = (f * temp * temp) % p;

  // sol exists if tr(a) = sum_{i=0}^{vt-1} a^2^i = 0 (mod p)
  trace (tr, a, p);

  if (!IsZero (tr))
    return -1;
  else
    {
      if (bound % 2)
	{
	  clear (x);
	  temp = a;
	  for (i = 0; i <= (bound - 1) / 2; ++i)
	    {
	      x += temp;
	      x %= p;

	      PowerMod (temp2, temp, 4, p);
	      temp = temp2;
	    }
	}
      else
	{
	  // randomly find c such that tr(c) = 1
          GF2EXModulus pmod;
          build(pmod,p);

          random(c, deg(p));
	  trace (tr, c, pmod);
	  while (!IsOne (tr))
	    {
              random(c, deg(p));
	      trace (tr, c, pmod);
	    }

	  clear (x);
	  atemp = a;
	  d = atemp;
	  ctemp = c;

	  for (i = 0; i < bound; ++i)
	    {
	      x += (d*ctemp);
	      x %= p;

	      ctemp = (ctemp*ctemp) % p;
	      atemp = (atemp*atemp) % p;

	      d += atemp;
	      d %= p;
	    }
	}

      x *= h;
      x %= p;

      // take x or x+h, whichever is smallest when evaluated at q
      GF2EX x2;
      ZZ xval, x2val;

      x2 = x+h;
      x2 %= p;

      xval = eval_poly(x,GF2E::cardinality());
      x2val = eval_poly(x2,GF2E::cardinality());

      if (x2val < xval)
        x = x2;

      return 1;
    }
}
