/**
 * @file qo_thresholds.hpp
 * @author Laurent Imbert
 * @brief  Various thresholds
 *
 */

#ifndef QO_THRESHOLDS_H
#define QO_THRESHOLDS_H

#include <NTL/GF2EX.h>
#include <NTL/ZZ_pX.h>
#include <NTL/lzz_pX.h>
#include <NTL/ZZ_pEX.h>
#include <NTL/lzz_pEX.h>

NTL_CLIENT

#define NUM_FIELDS 11
#define NUM_SIZES  7

template < typename T >
struct Thresholds {
  static int mulexact_crossover[NUM_FIELDS][NUM_SIZES];
  static int sqrexact_crossover[NUM_FIELDS][NUM_SIZES];

  static int pseudo_xgcd_crossover[NUM_FIELDS];
  static int pseudo_xgcd_left_crossover[NUM_FIELDS];
  static int pseudo_xgcd_partial_crossover[NUM_FIELDS][NUM_SIZES];

  static int hxgcd_inner_crossover;

  static int half_xgcd_crossover[NUM_FIELDS];
  static int half_xgcd_left_crossover[NUM_FIELDS];
  static int half_xgcd_partial_crossover[NUM_FIELDS][NUM_SIZES];

  static int reduce_crossover[NUM_FIELDS][NUM_SIZES];
  static int multiply_crossover[NUM_FIELDS];
  static int square_crossover[NUM_FIELDS];
  static int cube_crossover_cantor[NUM_FIELDS];
  static int cube_crossover_mulsqr[NUM_FIELDS];
  static int cube_crossover_cantor2[NUM_FIELDS];

  static int multiply_plain_crossover[NUM_FIELDS];
  static int square_plain_crossover[NUM_FIELDS];
  static int cube_plain_crossover_cantor[NUM_FIELDS];
  static int cube_plain_crossover_mulsqr[NUM_FIELDS];

  static int get_field_idx();


  //
  // Exact/partial multiplication thresholds
  //

  static int get_mulexact_idx(long size, long bound) { 
    int idx = NumBits(size-bound)-1;
    if (idx < 0)
      idx = 0;
    else if (idx >= NUM_SIZES)
      idx = NUM_SIZES-1;
    return idx;
  }

  static int get_mulexact_crossover(long size, long bound) {
    return mulexact_crossover[get_field_idx()][get_mulexact_idx(size,bound)];
  }

  static int get_sqrexact_crossover(long size, long bound) {
    return sqrexact_crossover[get_field_idx()][get_mulexact_idx(size,bound)];
  }



  //
  // XGCD, XGCD_LEFT, XGCD_PARTIAL thresholds
  //

  static int get_xgcd_partial_idx(long size, long bound) { 
    int idx = NumBits(size-bound)-1;
    if (idx < 0)
      idx = 0;
    else if (idx >= NUM_SIZES)
      idx = NUM_SIZES-1;
    return idx;
  }


  static int get_pseudo_xgcd_crossover() {
    return pseudo_xgcd_crossover[get_field_idx()];
  }

  static int get_pseudo_xgcd_left_crossover() {
    return pseudo_xgcd_left_crossover[get_field_idx()];
  }

  static int get_pseudo_xgcd_partial_crossover(long size, long bound) {
    /*
    cout << "PSEUDO CROSS:  size = " << size << ", bound = " << bound << endl;
    cout << "PC: fieldidx = " << get_field_idx() << ", idx = " << get_xgcd_partial_idx(size,bound) << ", cross = " << pseudo_xgcd_partial_crossover[get_field_idx()][get_xgcd_partial_idx(size,bound)] << endl;
    cout << "NB:  " << NumBits(size/bound)-1 << endl;
    */
    return pseudo_xgcd_partial_crossover[get_field_idx()][get_xgcd_partial_idx(size,bound)];
  }




  static int get_half_xgcd_crossover() {
    return half_xgcd_crossover[get_field_idx()];
  }

  static int get_half_xgcd_left_crossover() {
    return half_xgcd_left_crossover[get_field_idx()];
  }

  static int get_half_xgcd_partial_crossover(long size, long bound) {
    /*
    cout << "HALF CROSS:  size = " << size << ", bound = " << bound << endl;
    cout << "HC: fieldidx = " << get_field_idx() << ", idx = " << get_xgcd_partial_idx(size,bound) << ", cross = " << half_xgcd_partial_crossover[get_field_idx()][get_xgcd_partial_idx(size,bound)] << endl;
    */
    return half_xgcd_partial_crossover[get_field_idx()][get_xgcd_partial_idx(size,bound)];
  }



  //
  // ideal arithmetic thresholds
  //

   static int get_reduce_idx(long degree) { 
    int idx = NumBits(degree) - 2;
    if (idx < 0)
      idx = 0;
    else if (idx >= NUM_SIZES)
      idx = NUM_SIZES-1;
    return idx;
  }


  static int get_reduce_crossover(long degree) {
    return half_xgcd_partial_crossover[get_field_idx()][get_reduce_idx(degree)];
  }

  static int get_multiply_crossover() {
    return multiply_crossover[get_field_idx()];
  }

  static int get_square_crossover() {
    return square_crossover[get_field_idx()];
  }

  static int get_cube_crossover_cantor() {
    return cube_crossover_cantor[get_field_idx()];
  }

  static int get_cube_crossover_mulsqr() {
    return cube_crossover_mulsqr[get_field_idx()];
  }

  static int get_cube_crossover_cantor2() {
    return cube_crossover_cantor2[get_field_idx()];
  }


  //
  // ideal arithmetic thresholds (plain xgcd, mulexact)
  //

  static int get_multiply_plain_crossover() {
    return multiply_plain_crossover[get_field_idx()];
  }

  static int get_square_plain_crossover() {
    return square_plain_crossover[get_field_idx()];
  }

  static int get_cube_plain_crossover_cantor() {
    return cube_plain_crossover_cantor[get_field_idx()];
  }

  static int get_cube_plain_crossover_mulsqr() {
    return cube_plain_crossover_mulsqr[get_field_idx()];
  }

};

#endif // guard
