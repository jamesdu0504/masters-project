/**
 * @file common.hpp
 * @author Michael Jacobson
 * @brief General-purpose methods to interface with NTL.  All library files should include this file.
 */

#ifndef ANTL_COMMON_H
#define ANTL_COMMON_H

#include <cmath>
#include <NTL/ZZ.h>
#include <NTL/RR.h>
#include <NTL/GF2EX.h>
#include <NTL/ZZ_pX.h>
#include <NTL/lzz_pX.h>
#include <NTL/ZZ_pEX.h>
#include <NTL/lzz_pEX.h>


// We use the NTL namespace everywhere. Rather than have a using directive in
// every file, we just put it here, for convenience and clarity.
NTL_CLIENT


/**
 * @brief Suggested beginning of the definition of a function-like macro.
 * 
 * This helps to prevent heisenbugs caused by macros that have tricky semantics
 * when placed in, say, an if-block. The optimizer strips out the do-while.
 */
#define ANTL_BEGIN_MACRO_FUNCTION	do {

/**
 * @brief Used to close a macro definition started with \c ANTL_BEGIN_MACRO. 
 */
#define ANTL_END_MACRO_FUNCTION		} while(0)

// Include macros for debugging output, and general output.
#include "debug.hpp"

/**
 * @brief Null macro used to mark a by-reference argument pass.
 *
 * This macro improves human readability by indicating that a method call
 * may change an argument in the caller scope.
 */
#define REF

// Forward declarations for prototypes in this file.
namespace NTL {
  //
  // NTL-like methods, for completeness.
  //
  inline void clear (long &X)	        { X = 0; }
  inline void set (long &X)          { X = 1; }
  inline void clear (float &i)       { i = 0.0f; }
  inline void set (float &i)         { i = 1.0f; }
  inline void clear (double &i)      { i = double(0); }
  inline void set (double &i)	        { i = double(1); }
  inline void clear (quad_float & i) { i = double(0); }
  inline void set (quad_float & i)	{ i = double(1); }


  inline long IsOne (const long &X)         { return (X == 1); }
  inline long IsZero (const long &X)        { return (X == 0); }
  inline long IsOdd (const long &X)         { return (X & 1); }

  inline long IsOne (const float &X)         { return (X == 1.0f); }
  inline long IsZero (const float &X)        { return (X == 0.0f); }

  inline long IsOne (const double &X)         { return (X == double(1)); }
  inline long IsZero (const double &X)        { return (X == double(0)); }

  inline long IsOne (const quad_float &X)         { return (X == double(1)); }
  inline long IsZero (const quad_float &X)        { return (X == double(0)); }

  // JLM: needed to compile quadratic_form. I have not tested them
  // for correctness.
  inline long deg(long X)              { return 0; }
  inline long deg(const ZZ & X)        { return 0; }
  inline long LeadCoeff(long X)        { return X; }
  inline ZZ LeadCoeff(const ZZ & X)    { return X; }
  inline void MakeMonic(long & X)      { X = 1; }
  inline void MakeMonic(ZZ & X)        { X = to_ZZ(1); }

  // assign(C,A) method for NTL classes that only support the assignment operator
  // (for compatibility with templated classes like exponentiation)
  template < class T >
  inline void assign(T &C, const T &A) { C = A; }

  // procedural arithmetic operations for standard types, to increase compatibility with NTL
  inline void add ( long &C, const long &A, const long &B )     { C = A + B; }
  inline void add ( float &C, const float &A, const float &B )     { C = A + B; }
  inline void add ( double &C, const double &A, const double &B )     { C = A + B; }

  inline void sub ( long &C, const long &A, const long &B )     { C = A - B; }
  inline void sub ( float &C, const float &A, const float &B )     { C = A - B; }
  inline void sub ( double &C, const double &A, const double &B )     { C = A - B; }

  inline void mul ( long &C, const long &A, const long &B )     { C = A * B; }
  inline void mul ( float &C, const float &A, const float &B )     { C = A * B; }
  inline void mul ( double &C, const double &A, const double &B )     { C = A * B; }

  inline void div ( long &C, const long &A, const long &B )     { C = A / B; }
  inline void div ( float &C, const float &A, const float &B )     { C = A / B; }
  inline void div ( double &C, const double &A, const double &B )     { C = A / B; }

  inline void sqr ( long &C, const long &A)     { C = A * A; }
  inline void sqr ( float &C, const float &A)     { C = A * A; }
  inline void sqr ( double &C, const double &A)     { C = A * A; }

}

namespace ANTL {
  inline int SqrRoot (const int & a) { return (int) ::floor(::sqrt((double) a)); }

  void DivRem (long &q, long &r, long a, long b);
  inline long SqrRoot (const long &a) { return (long) ::floor(::sqrt((double) a)); }


  //
  // cmath-like methods, for completeness.
  // Most of these exist because of problems where gcc doesn't find the
  // cmath methods, presumably because of namespace issues).
  //
  inline float sqrt(float x) { return (float)std::sqrt((double)x); }
  inline double sqrt(double x) { return std::sqrt(x); }
  inline long abs(long x) { return std::abs(x); }
  inline float abs(float x) { return (float)std::abs((double)x); }
  inline double abs(double x) { return std::abs(x); }
  inline quad_float abs(const quad_float & x) { return NTL::fabs(x); }
  inline float exp(float x) { return (float)std::exp((double) x); }
  inline double exp(double x) { return std::exp(x); }
  inline float log(float x) { return (float)std::log((double)x); }
  inline double log(double x) { return std::log(x); }



  // finite field cardinality macros
  template <class> ZZ CARDINALITY(void);

  template <> inline ZZ CARDINALITY < ZZ > (void) {
    return to_ZZ(0);
  }

  template <> inline ZZ CARDINALITY < long > (void) {
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



  //
  // Template-friendly methods for common conversions & constants.
  // e.g.: ZZ x = to<ZZ>(myLongVariable);
  //


  template <class T> T to(const int& a)          { return T(a); }
  template <class T> T to(const long& a)         { return T(a); }
  template <class T> T to(const float& a)        { return T(a); }
  template <class T> T to(const double& a)       { return T(a); }
  template <class T> T to(const ZZ& a)           { return T(a); }
  template <class T> T to(const RR& a)           { return T(a); }
  template <class T> T to(const quad_float & a)  { return T(a); }
  template <class T> T to(const GF2EX & a)       { return T(a); }
  template <class T> T to(const zz_pX & a)       { return T(a); }
  template <class T> T to(const ZZ_pX & a)       { return T(a); }
  template <class T> T to(const zz_pEX & a)      { return T(a); }
  template <class T> T to(const ZZ_pEX & a)      { return T(a); }

  template<> inline int to<int>(const int& a)        { return a; }
  template<> inline int to<int>(const long& a)       { return to_int(a); }
  template<> inline int to<int>(const float& a)      { return to_int(a); }
  template<> inline int to<int>(const double& a)     { return to_int(a); }
  template<> inline int to<int>(const ZZ& a)         { return to_int(a); }
  template<> inline int to<int>(const RR& a)         { return to_int(a); }
  template<> inline int to<int>(const quad_float& a) { return to_int(a); }

  template<> inline long to<long>(const int& a)        { return to_long(a); }
  template<> inline long to<long>(const long& a)       { return a; }
  template<> inline long to<long>(const float& a)      { return to_long(a); }
  template<> inline long to<long>(const double& a)     { return to_long(a); }
  template<> inline long to<long>(const ZZ& a)         { return to_long(a); }
  template<> inline long to<long>(const RR& a)         { return to_long(a); }
  template<> inline long to<long>(const quad_float& a) { return to_long(a); }

  template<> inline float to<float>(const int& a)        { return to_float(a); }
  template<> inline float to<float>(const long& a)       { return to_float(a); }
  template<> inline float to<float>(const float& a)      { return a; }
  template<> inline float to<float>(const double& a)     { return to_float(a); }
  template<> inline float to<float>(const ZZ& a)         { return to_float(a); }
  template<> inline float to<float>(const RR& a)         { return to_float(a); }
  template<> inline float to<float>(const quad_float& a) { return to_float(a); }

  template<> inline double to<double>(const int& a)        { return to_double(a); }
  template<> inline double to<double>(const long& a)       { return to_double(a); }
  template<> inline double to<double>(const float& a)      { return to_double(a); }
  template<> inline double to<double>(const double& a)     { return a; }
  template<> inline double to<double>(const ZZ& a)         { return to_double(a); }
  template<> inline double to<double>(const RR& a)         { return to_double(a); }
  template<> inline double to<double>(const quad_float& a) { return to_double(a); }

  template<> inline ZZ to<ZZ>(const int& a)        { return to_ZZ(a); }
  template<> inline ZZ to<ZZ>(const long& a)       { return to_ZZ(a); }
  template<> inline ZZ to<ZZ>(const float& a)      { return to_ZZ(a); }
  template<> inline ZZ to<ZZ>(const double& a)     { return to_ZZ(a); }
  template<> inline ZZ to<ZZ>(const ZZ& a)         { return a; }
  template<> inline ZZ to<ZZ>(const RR& a)         { return to_ZZ(a); }
  template<> inline ZZ to<ZZ>(const quad_float& a) { return to_ZZ(a); }

  template<> inline RR to<RR>(const int& a)        { return to_RR(a); }
  template<> inline RR to<RR>(const long& a)       { return to_RR(a); }
  template<> inline RR to<RR>(const float& a)      { return to_RR(a); }
  template<> inline RR to<RR>(const double& a)     { return to_RR(a); }
  template<> inline RR to<RR>(const ZZ& a)         { return to_RR(a); }
  template<> inline RR to<RR>(const RR& a)         { return a; }
  template<> inline RR to<RR>(const quad_float& a) { return to_RR(a); }

  template<> inline quad_float to<quad_float>(const int& a)        { return to_quad_float(a); }
  template<> inline quad_float to<quad_float>(const long& a)       { return to_quad_float(a); }
  template<> inline quad_float to<quad_float>(const float& a)      { return to_quad_float(a); }
  template<> inline quad_float to<quad_float>(const double& a)     { return to_quad_float(a); }
  template<> inline quad_float to<quad_float>(const ZZ& a)         { return to_quad_float(a); }
  template<> inline quad_float to<quad_float>(const RR& a)         { return to_quad_float(a); }
  template<> inline quad_float to<quad_float>(const quad_float& a) { return a; }

  template<> inline GF2EX to<GF2EX>(const int& a)        { return to_GF2EX(to_ZZ(a)); }
  template<> inline GF2EX to<GF2EX>(const long& a)       { return to_GF2EX(a); }
  template<> inline GF2EX to<GF2EX>(const float& a)      { return to_GF2EX(to_ZZ(a)); }
  template<> inline GF2EX to<GF2EX>(const double& a)     { return to_GF2EX(to_ZZ(a)); }
  template<> inline GF2EX to<GF2EX>(const ZZ& a)         { return to_GF2EX(a); }
  template<> inline GF2EX to<GF2EX>(const RR& a)         { return to_GF2EX(to_ZZ(a)); }
  template<> inline GF2EX to<GF2EX>(const quad_float& a) { return to_GF2EX(to_ZZ(a)); }
  template<> inline GF2EX to<GF2EX>(const GF2EX& a) { return a; }

  template<> inline ZZ_pX to<ZZ_pX>(const ZZ_pX& a) { return a; }
  template<> inline ZZ_pEX to<ZZ_pEX>(const ZZ_pEX& a) { return a; }
  template<> inline zz_pX to<zz_pX>(const zz_pX& a) { return a; }
  template<> inline zz_pEX to<zz_pEX>(const zz_pEX& a) { return a; }

  template<> inline zz_pX to<zz_pX>(const int & a) { return zz_pX(a,0); }
} // ANTL

/*
// Unspecialized template definitions.
#include "../../src/common_impl.hpp"
*/

#endif // guard
