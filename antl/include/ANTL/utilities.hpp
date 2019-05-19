/**
 * @file utilities.hpp
 * @author Michael Jacobson
 * @brief Utilitiy routines supporting arithmetic with reduced ideals of
 * number and function fields.  The main functionalities are
 *   - eval_poly and get_poly_modq: routines for converting to and from
 *     integers and polynomials
 *   - Jacobi: quadratic residuosity
 *   - ressol: roots of modular quadratic equations
 */

#ifndef UTILITIES_H
#define UTILITIES_H

#include <gmp.h>
#include <NTL/ZZ.h>
#include <NTL/GF2EX.h>
#include <NTL/ZZ_pX.h>
#include <NTL/lzz_pX.h>
#include <NTL/ZZ_pEX.h>
#include <NTL/lzz_pEX.h>

#include <NTL/ZZ_pXFactoring.h>
#include <NTL/lzz_pXFactoring.h>
#include <NTL/ZZ_pEXFactoring.h>
#include <NTL/lzz_pEXFactoring.h>
#include <NTL/GF2XFactoring.h>
#include <NTL/GF2EXFactoring.h>



// We use the NTL namespace everywhere. Rather than have a using directive in
// every file, we just put it here, for convenience and clarity.
NTL_CLIENT


//
// Utility routines for polynomials.
//

/**
 * @fn get_poly_modq
 * Computes a polynomial mod q corresponding to the integer X (q-adic
 * representation of X).
 * 
 * @param p Reference to an instance of GF2EX to hold the newly
 * created polynomial.
 * @param X The integer to create the polynomial from.
 * @param q The modulus of the polynomial creation.
 */
void get_poly_modq (GF2X & p, const ZZ & X, const ZZ & q);
void get_poly_modq (GF2EX & p, const ZZ & X, const ZZ & q);
void get_poly_modq (ZZ_pX & p, const ZZ & X, const ZZ & q);
void get_poly_modq (ZZ_pEX & p, const ZZ & X, const ZZ & q);
void get_poly_modq (zz_pX & p, const ZZ & X, const ZZ & q);
void get_poly_modq (zz_pEX & p, const ZZ & X, const ZZ & q);


/**
 * @brief Computes the integer X corresponding to the q-adic
 * representation of p.
 */
template < class T >
ZZ eval_poly (const T & p, const ZZ & q);
template <> ZZ eval_poly <ZZ_pEX> (const ZZ_pEX & p, const ZZ & q);
template <> ZZ eval_poly <zz_pEX> (const zz_pEX & p, const ZZ & q);
template <> ZZ eval_poly <GF2EX> (const GF2EX & p, const ZZ & q);



//
// quadratic residuosity functions
//

/* Jacobi functions - assumes a is reduced mod n */
long Jacobi_base (const ZZ & a, const ZZ & n);
long Jacobi_base (const long & a, const long & n);

/* Jacobi functions - no preconditions on a and n */
long Jacobi(const ZZ & a, const ZZ & n);
long Jacobi(const long & a, const long & n);
long Jacobi(const ZZ_pX & a, const ZZ_pX & n);
long Jacobi(const zz_pX & a, const zz_pX & n);
long Jacobi(const ZZ_pEX & a, const ZZ_pEX & n);
long Jacobi(const zz_pEX & a, const zz_pEX & n);

/* GF2X quadratic character */
/**
 * @brief Computes the trace of a (mod p), i.e.,
 * 		tr = sum_{i=0}^{vt-1} a^2^i (mod p)
 *  where deg(p) = v.
 */
void trace (GF2X & t, const GF2X & a, const GF2X & pmod);
long Jacobi (const GF2X & h, const GF2X & f, const GF2X & n);

/* GF2EX quadratic character */
/**
 * @brief Computes the trace of a (mod p), i.e.,
 * 		tr = sum_{i=0}^{vt-1} a^2^i (mod p)
 *  where deg(p) = v.
 */
void trace (GF2EX & t, const GF2EX & a, const GF2EXModulus & pmod);
long Jacobi (const GF2EX & h, const GF2EX & f, const GF2EX & n);



//
// modular square root functions
//

//
// get_qdp (computes a square root of Delta mod p)
//   - used by templated ressol)
//
void get_qdp(ZZ & q, const ZZ_pX & p);
void get_qdp(ZZ & q, const zz_pX & p);
void get_qdp(ZZ & q, const ZZ_pEX & p);
void get_qdp(ZZ & q, const zz_pEX & p);

/**
 * @brief Computes x, a a solution of x^2 = a (mod p).
 * @return
 * 		 0  if the equation has a unique solution (p | a)
 * 		-1  if x^2 = a (mod p) has no solutions
 * 		 1  if x^2 = a (mod p) has 2 solutions
 */
template < class T >
long ressol (T & x, const T & a, const T & p);



/**
 * @brief Computes x, a a solution of x^2 + hx -f = 0 (mod p).
 * @return
 * 		 0  if the equation has a unique solution (p | h)
 * 		-1  if x^2 + hx -f = 0 (mod p) has no solutions
 * 		 1  if x^2 + hx -f = 0 (mod p) has 2 solutions
 */
long ressol (GF2EX & x, const GF2EX & h, const GF2EX & f, const GF2EX & p);



/**
 * @brief Useful routines for access to the size of a finite field
 */

template <class> ZZ CARDINALITY(void);

template <> inline ZZ CARDINALITY < ZZ > (void) {
  return to_ZZ(0);
}

template <> inline ZZ CARDINALITY < ZZ_p > (void) {
  return ZZ_p::modulus ();
}

template <> inline ZZ CARDINALITY < zz_p > (void) {
  return to_ZZ (zz_p::modulus ());
}

template <> inline ZZ CARDINALITY < ZZ_pX > (void) {
  return ZZ_p::modulus ();
}

template <> inline ZZ CARDINALITY < zz_pX > (void) {
  return to_ZZ (zz_p::modulus ());
}

template <> inline ZZ CARDINALITY < ZZ_pE > (void) {
  return ZZ_pE::cardinality ();
}

template <> inline ZZ CARDINALITY < zz_pE > (void) {
  return zz_pE::cardinality ();
}

template <> inline ZZ CARDINALITY < ZZ_pEX > (void) {
  return ZZ_pE::cardinality ();
}

template <> inline ZZ CARDINALITY < zz_pEX > (void) {
  return zz_pE::cardinality ();
}

template <> inline ZZ CARDINALITY < GF2E > (void) {
  return GF2E::cardinality ();
}

template <> inline ZZ CARDINALITY < GF2EX > (void) {
  return GF2E::cardinality ();
}


// Unspecialized template definitions.
#include "../src/utilities_impl.hpp"

#endif // guard
