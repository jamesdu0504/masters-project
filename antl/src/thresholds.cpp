/**
 * @file thresholds.cpp
 * @author Michael Jacobson
 * @remark threshold values controlling various algorithm choices
 */

#include <ANTL/thresholds.hpp>


template <> 
int Thresholds<GF2EX>::get_field_idx() 
{
  int idx = NumBits(GF2E::degree()) - 2;
  if (idx < 0)
    idx = 0;
  else if (idx >= NUM_FIELDS)
    idx = NUM_FIELDS-1;
  return idx;
}

template <> 
int Thresholds<ZZ_pX>::get_field_idx() 
{
  int idx = NumBits(ZZ_p::modulus()) - 2;
  if (idx < 0)
    idx = 0;
  else if (idx >= NUM_FIELDS)
    idx = NUM_FIELDS-1;
  return idx;
}

template <> 
int Thresholds<zz_pX>::get_field_idx() 
{
  int idx = NumBits(zz_p::modulus()) - 2;
  if (idx < 0)
    idx = 0;
  else if (idx >= NUM_FIELDS)
    idx = NUM_FIELDS-1;
  return idx;
}

template <> 
int Thresholds<ZZ_pEX>::get_field_idx() 
{
  int idx = NumBits(ZZ_pE::degree()) - 2;
  if (idx < 0)
    idx = 0;
  else if (idx >= NUM_FIELDS)
    idx = NUM_FIELDS-1;
  return idx;
}

template <> 
int Thresholds<zz_pEX>::get_field_idx() 
{
  int idx = NumBits(zz_pE::degree()) - 2;
  if (idx < 0)
    idx = 0;
  else if (idx >= NUM_FIELDS)
    idx = NUM_FIELDS-1;
  return idx;
}



// These thresholds control the calls in MulExact
template<> int Thresholds<GF2EX>::mulexact_crossover[NUM_FIELDS][NUM_SIZES] = 
{ {3, 4, 13, 35, 100, 0, 0},
  {3, 4, 14, 30, 100, 0, 0},
  {2, 4, 10, 25, 70, 0, 0},
  {2, 4, 7, 20, 45, 0, 0},
  {2, 4, 5, 11, 35, 90, 0},
  {2, 3, 5, 9, 20, 60, 0},
  {2, 3, 6, 9, 35, 70, 0},
  {2, 3, 7, 15, 35, 80, 0},
  {2, 3, 7, 15, 40, 80, 0},
  {2, 3, 7, 15, 35, 90, 0}, 
  {2, 4, 7, 15, 40, 90, 0}};

template<> int Thresholds<ZZ_pX>::mulexact_crossover[NUM_FIELDS][NUM_SIZES] = 
{ {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0}};

template<> int Thresholds<zz_pX>::mulexact_crossover[NUM_FIELDS][NUM_SIZES] = 
{ {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0}};

template<> int Thresholds<ZZ_pEX>::mulexact_crossover[NUM_FIELDS][NUM_SIZES] = 
{ {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0}};

template<> int Thresholds<zz_pEX>::mulexact_crossover[NUM_FIELDS][NUM_SIZES] = 
{ {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0}};



// These thresholds control the calls in SqrExact
template<> int Thresholds<GF2EX>::sqrexact_crossover[NUM_FIELDS][NUM_SIZES] = 
{ {100000000,100000000,100000000,100000000,100000000,100000000,100000000},
  {100000000,100000000,100000000,100000000,100000000,100000000,100000000},
  {100000000,100000000,100000000,100000000,100000000,100000000,100000000},
  {100000000,100000000,100000000,100000000,100000000,100000000,100000000},
  {100000000,100000000,100000000,100000000,100000000,100000000,100000000},
  {100000000,100000000,100000000,100000000,100000000,100000000,100000000},
  {100000000,100000000,100000000,100000000,100000000,100000000,100000000},
  {100000000,100000000,100000000,100000000,100000000,100000000,100000000},
  {100000000,100000000,100000000,100000000,100000000,100000000,100000000},
  {100000000,100000000,100000000,100000000,100000000,100000000,100000000},
  {100000000,100000000,100000000,100000000,100000000,100000000,100000000}};

template<> int Thresholds<ZZ_pX>::sqrexact_crossover[NUM_FIELDS][NUM_SIZES] = 
{ {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0}};

template<> int Thresholds<zz_pX>::sqrexact_crossover[NUM_FIELDS][NUM_SIZES] = 
{ {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0}};

template<> int Thresholds<ZZ_pEX>::sqrexact_crossover[NUM_FIELDS][NUM_SIZES] = 
{ {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0}};

template<> int Thresholds<zz_pEX>::sqrexact_crossover[NUM_FIELDS][NUM_SIZES] = 
{ {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0}};




// These thresholds control the calls in HXGCD, thus is less sensitive
template<> int Thresholds<GF2EX>::hxgcd_inner_crossover = 18;
template<> int Thresholds<ZZ_pX>::hxgcd_inner_crossover = 18;
template<> int Thresholds<zz_pX>::hxgcd_inner_crossover = 18;
template<> int Thresholds<ZZ_pEX>::hxgcd_inner_crossover = 18;
template<> int Thresholds<zz_pEX>::hxgcd_inner_crossover = 18;



// These thresholds control whether plain or pseudodivision xgcd is used
template<> int Thresholds<GF2EX>::pseudo_xgcd_crossover[NUM_FIELDS] = {0,0,0,0,0,0,15,15,15,15,15};
template<> int Thresholds<ZZ_pX>::pseudo_xgcd_crossover[NUM_FIELDS] = {0,0,0,0,0,0,0,0,0,0,0};
template<> int Thresholds<zz_pX>::pseudo_xgcd_crossover[NUM_FIELDS] = {0,0,0,0,0,0,0,0,0,0,0};
template<> int Thresholds<ZZ_pEX>::pseudo_xgcd_crossover[NUM_FIELDS] = {0,0,0,0,0,0,0,0,0,0,0};
template<> int Thresholds<zz_pEX>::pseudo_xgcd_crossover[NUM_FIELDS] = {0,0,0,0,0,0,0,0,0,0,0};

template<> int Thresholds<GF2EX>::pseudo_xgcd_left_crossover[NUM_FIELDS] = {0,0,0,0,0,0,20,20,20,20,20};
template<> int Thresholds<ZZ_pX>::pseudo_xgcd_left_crossover[NUM_FIELDS] = {0,0,0,0,0,0,0,0,0,0,0};
template<> int Thresholds<zz_pX>::pseudo_xgcd_left_crossover[NUM_FIELDS] = {0,0,0,0,0,0,0,0,0,0,0};
template<> int Thresholds<ZZ_pEX>::pseudo_xgcd_left_crossover[NUM_FIELDS] = {0,0,0,0,0,0,0,0,0,0,0};
template<> int Thresholds<zz_pEX>::pseudo_xgcd_left_crossover[NUM_FIELDS] = {0,0,0,0,0,0,0,0,0,0,0};

template<> int Thresholds<GF2EX>::pseudo_xgcd_partial_crossover[NUM_FIELDS][NUM_SIZES] = 
{ {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 6, 10, 0, 0, 0},
  {0, 0, 9, 13, 0, 0, 0},
  {0, 0, 9, 12, 0, 0, 0},
  {0, 0, 9, 13, 0, 0, 0},
  {0, 0, 9, 14, 0, 0, 0}};

template<> int Thresholds<ZZ_pX>::pseudo_xgcd_partial_crossover[NUM_FIELDS][NUM_SIZES] = 
{ {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0}};

template<> int Thresholds<zz_pX>::pseudo_xgcd_partial_crossover[NUM_FIELDS][NUM_SIZES] = 
{ {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0}};

template<> int Thresholds<ZZ_pEX>::pseudo_xgcd_partial_crossover[NUM_FIELDS][NUM_SIZES] = 
{ {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0}};

template<> int Thresholds<zz_pEX>::pseudo_xgcd_partial_crossover[NUM_FIELDS][NUM_SIZES] = 
{ {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0}};




// These thresholds control whether iterative of half xgcd is used
template<> int Thresholds<GF2EX>::half_xgcd_crossover[NUM_FIELDS] = {30,30,30,30,30,30,60,60,60,60,60};
template<> int Thresholds<ZZ_pX>::half_xgcd_crossover[NUM_FIELDS] = {0,0,0,0,0,0,0,0,0,0,0};
template<> int Thresholds<zz_pX>::half_xgcd_crossover[NUM_FIELDS] = {0,0,0,0,0,0,0,0,0,0,0};
template<> int Thresholds<ZZ_pEX>::half_xgcd_crossover[NUM_FIELDS] = {0,0,0,0,0,0,0,0,0,0,0};
template<> int Thresholds<zz_pEX>::half_xgcd_crossover[NUM_FIELDS] = {0,0,0,0,0,0,0,0,0,0,0};

template<> int Thresholds<GF2EX>::half_xgcd_left_crossover[NUM_FIELDS] = {60,60,60,60,60,60,210,210,210,210,210};
template<> int Thresholds<ZZ_pX>::half_xgcd_left_crossover[NUM_FIELDS] = {0,0,0,0,0,0,0,0,0,0,0};
template<> int Thresholds<zz_pX>::half_xgcd_left_crossover[NUM_FIELDS] = {0,0,0,0,0,0,0,0,0,0,0};
template<> int Thresholds<ZZ_pEX>::half_xgcd_left_crossover[NUM_FIELDS] = {0,0,0,0,0,0,0,0,0,0,0};
template<> int Thresholds<zz_pEX>::half_xgcd_left_crossover[NUM_FIELDS] = {0,0,0,0,0,0,0,0,0,0,0};

template<> int Thresholds<GF2EX>::half_xgcd_partial_crossover[NUM_FIELDS][NUM_SIZES] = 
{ {30, 30, 140, 0, 0, 70, 70},
  {30, 60, 0, 0, 0, 60, 70},
  {30, 100, 0, 0, 0, 60, 70},
  {30, 60, 0, 0, 0, 60, 70},
  {30, 90, 0, 0, 0, 60, 70},
  {30, 90, 200, 190, 200, 70, 70},
  {30, 90, 120, 0, 0, 0, 160},
  {50, 70, 0, 0, 0, 160, 120},
  {30, 30, 0, 0, 0, 150, 110},
  {30, 30, 0, 30, 0, 150, 110},
  {30, 40, 160, 0, 0, 110, 110}};

template<> int Thresholds<ZZ_pX>::half_xgcd_partial_crossover[NUM_FIELDS][NUM_SIZES] = 
{ {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0}};

template<> int Thresholds<zz_pX>::half_xgcd_partial_crossover[NUM_FIELDS][NUM_SIZES] = 
{ {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0}};

template<> int Thresholds<ZZ_pEX>::half_xgcd_partial_crossover[NUM_FIELDS][NUM_SIZES] = 
{ {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0}};

template<> int Thresholds<zz_pEX>::half_xgcd_partial_crossover[NUM_FIELDS][NUM_SIZES] = 
{ {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0}};




template<> int Thresholds<GF2EX>::reduce_crossover[NUM_FIELDS][NUM_SIZES] =
{ {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0}};

template<> int Thresholds<ZZ_pX>::reduce_crossover[NUM_FIELDS][NUM_SIZES] =
{ {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0}};

template<> int Thresholds<zz_pX>::reduce_crossover[NUM_FIELDS][NUM_SIZES] = 
{ {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0}};

template<> int Thresholds<ZZ_pEX>::reduce_crossover[NUM_FIELDS][NUM_SIZES] = 
{ {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0}};

template<> int Thresholds<zz_pEX>::reduce_crossover[NUM_FIELDS][NUM_SIZES] = 
{ {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0}};


template<> int Thresholds<GF2EX>::multiply_crossover[NUM_FIELDS] = {10, 10, 7, 7, 6, 6, 6, 7, 8, 9, 9};
template<> int Thresholds<ZZ_pX>::multiply_crossover[NUM_FIELDS] = {0,0,0,0,0,0,0,0,0,0,0};
template<> int Thresholds<zz_pX>::multiply_crossover[NUM_FIELDS] = {0,0,0,0,0,0,0,0,0,0,0};
template<> int Thresholds<ZZ_pEX>::multiply_crossover[NUM_FIELDS] = {0,0,0,0,0,0,0,0,0,0,0};
template<> int Thresholds<zz_pEX>::multiply_crossover[NUM_FIELDS] = {0,0,0,0,0,0,0,0,0,0,0};

template<> int Thresholds<GF2EX>::square_crossover[NUM_FIELDS] = {5, 5, 5, 6, 6, 6, 6, 5, 5, 5, 7};
template<> int Thresholds<ZZ_pX>::square_crossover[NUM_FIELDS] = {0,0,0,0,0,0,0,0,0,0,0};
template<> int Thresholds<zz_pX>::square_crossover[NUM_FIELDS] = {0,0,0,0,0,0,0,0,0,0,0};
template<> int Thresholds<ZZ_pEX>::square_crossover[NUM_FIELDS] = {0,0,0,0,0,0,0,0,0,0,0};
template<> int Thresholds<zz_pEX>::square_crossover[NUM_FIELDS] = {0,0,0,0,0,0,0,0,0,0,0};

template<> int Thresholds<GF2EX>::cube_crossover_cantor[NUM_FIELDS] = {5,6,6,5,5,5,5,5,0,0,0};
template<> int Thresholds<ZZ_pX>::cube_crossover_cantor[NUM_FIELDS] = {0,0,0,0,0,0,0,0,0,0,0};
template<> int Thresholds<zz_pX>::cube_crossover_cantor[NUM_FIELDS] = {0,0,0,0,0,0,0,0,0,0,0};
template<> int Thresholds<ZZ_pEX>::cube_crossover_cantor[NUM_FIELDS] = {0,0,0,0,0,0,0,0,0,0,0};
template<> int Thresholds<zz_pEX>::cube_crossover_cantor[NUM_FIELDS] = {0,0,0,0,0,0,0,0,0,0,0};

template<> int Thresholds<GF2EX>::cube_crossover_mulsqr[NUM_FIELDS] = {12,11,11,13,9,9,9,9,0,0,0};
template<> int Thresholds<ZZ_pX>::cube_crossover_mulsqr[NUM_FIELDS] = {0,0,0,0,0,0,0,0,0,0,0};
template<> int Thresholds<zz_pX>::cube_crossover_mulsqr[NUM_FIELDS] = {0,0,0,0,0,0,0,0,0,0,0};
template<> int Thresholds<ZZ_pEX>::cube_crossover_mulsqr[NUM_FIELDS] = {0,0,0,0,0,0,0,0,0,0,0};
template<> int Thresholds<zz_pEX>::cube_crossover_mulsqr[NUM_FIELDS] = {0,0,0,0,0,0,0,0,0,0,0};

template<> int Thresholds<GF2EX>::cube_crossover_cantor2[NUM_FIELDS] = {55,45,40,30,45,30,0,0,0,0,0};
template<> int Thresholds<ZZ_pX>::cube_crossover_cantor2[NUM_FIELDS] = {0,0,0,0,0,0,0,0,0,0,0};
template<> int Thresholds<zz_pX>::cube_crossover_cantor2[NUM_FIELDS] = {0,0,0,0,0,0,0,0,0,0,0};
template<> int Thresholds<ZZ_pEX>::cube_crossover_cantor2[NUM_FIELDS] = {0,0,0,0,0,0,0,0,0,0,0};
template<> int Thresholds<zz_pEX>::cube_crossover_cantor2[NUM_FIELDS] = {0,0,0,0,0,0,0,0,0,0,0};



template<> int Thresholds<GF2EX>::multiply_plain_crossover[NUM_FIELDS] = {8, 8, 8, 9, 9, 9, 10, 9, 9, 9, 10};
template<> int Thresholds<ZZ_pX>::multiply_plain_crossover[NUM_FIELDS] = {0,0,0,0,0,0,0,0,0,0,0};
template<> int Thresholds<zz_pX>::multiply_plain_crossover[NUM_FIELDS] = {0,0,0,0,0,0,0,0,0,0,0};
template<> int Thresholds<ZZ_pEX>::multiply_plain_crossover[NUM_FIELDS] = {0,0,0,0,0,0,0,0,0,0,0};
template<> int Thresholds<zz_pEX>::multiply_plain_crossover[NUM_FIELDS] = {0,0,0,0,0,0,0,0,0,0,0};

template<> int Thresholds<GF2EX>::square_plain_crossover[NUM_FIELDS] = {6, 6, 6, 7, 8, 9, 7, 8, 8, 9, 10};
template<> int Thresholds<ZZ_pX>::square_plain_crossover[NUM_FIELDS] = {0,0,0,0,0,0,0,0,0,0,0};
template<> int Thresholds<zz_pX>::square_plain_crossover[NUM_FIELDS] = {0,0,0,0,0,0,0,0,0,0,0};
template<> int Thresholds<ZZ_pEX>::square_plain_crossover[NUM_FIELDS] = {0,0,0,0,0,0,0,0,0,0,0};
template<> int Thresholds<zz_pEX>::square_plain_crossover[NUM_FIELDS] = {0,0,0,0,0,0,0,0,0,0,0};

template<> int Thresholds<GF2EX>::cube_plain_crossover_cantor[NUM_FIELDS] = {7,7,7,7,7,7,7,7,7,7,7};
template<> int Thresholds<ZZ_pX>::cube_plain_crossover_cantor[NUM_FIELDS] = {0,0,0,0,0,0,0,0,0,0,0};
template<> int Thresholds<zz_pX>::cube_plain_crossover_cantor[NUM_FIELDS] = {0,0,0,0,0,0,0,0,0,0,0};
template<> int Thresholds<ZZ_pEX>::cube_plain_crossover_cantor[NUM_FIELDS] = {0,0,0,0,0,0,0,0,0,0,0};
template<> int Thresholds<zz_pEX>::cube_plain_crossover_cantor[NUM_FIELDS] = {0,0,0,0,0,0,0,0,0,0,0};

template<> int Thresholds<GF2EX>::cube_plain_crossover_mulsqr[NUM_FIELDS] = {12,13,13,13,10,10,10,10,10,10,10};
template<> int Thresholds<ZZ_pX>::cube_plain_crossover_mulsqr[NUM_FIELDS] = {0,0,0,0,0,0,0,0,0,0,0};
template<> int Thresholds<zz_pX>::cube_plain_crossover_mulsqr[NUM_FIELDS] = {0,0,0,0,0,0,0,0,0,0,0};
template<> int Thresholds<ZZ_pEX>::cube_plain_crossover_mulsqr[NUM_FIELDS] = {0,0,0,0,0,0,0,0,0,0,0};
template<> int Thresholds<zz_pEX>::cube_plain_crossover_mulsqr[NUM_FIELDS] = {0,0,0,0,0,0,0,0,0,0,0};
