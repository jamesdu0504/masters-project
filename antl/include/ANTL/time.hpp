/**
 * @file time.hpp
 * @author Michael Jacobson
 * @brief Utilitiy routine for timing.
 */

#ifndef QO_TIME_H
#define QO_TIME_H

// -----------------------------------------------------------------------------
//   Returns the delay between t1 and t2 in milliseconds
// -----------------------------------------------------------------------------
double delay_msec (clock_t t1, clock_t t2)
{
  const double k = 1000.0/((double)CLOCKS_PER_SEC);
  double diff = difftime(t2, t1);
  return k * diff;
}
// -----------------------------------------------------------------------------

#endif
