/*
 * ExponentiationYao.hpp
 *
 *  Created on: Sep 16, 2015
 *      Author: Anton Mosunov
 */

#ifndef EXPONENTIATION_YAO_HPP_H
#define EXPONENTIATION_YAO_HPP_H

#include <vector>
#include "Exponentiation.hpp"

namespace ANTL
{

  /**
   * @brief Class for Yao's exponentiation algorithm
   * @remarks Given an element A and the power n, works by precomputing A, A^2, ... , A^{2^m}
   * for certain m>=0 first, m <= lg(n), and then evaluating A^n. Designed for evaluating
   * multiple powers.
   */
  template < class T >
  class ExponentiationYao : public Exponentiation<T>
  {
  protected:
	  vector<T> powersOf2;  /**< vector containing A, A^2, ... , A^{2^m} for some m >= 0 */

  public:
    ExponentiationYao() {};
    ~ExponentiationYao() {};

    /**
     * @brief Computes A^n by first precomputing A, A^2, ... , A^{2^m} for certain m >= 0.
     * Note that after previous exponentiations some of the powers might already be
     * computed and stored in memory. In case if there's not enough powers, the algorithm
     * computes the rest and stores the resulting powers to memory, ready for use in
     * subsequent exponentiations.
     * @param[out] C result of computing A^n
     * @param[in] A base for exponentiation
     * @param[in] n exponent
     */
    void power (T &C, const T &A, const ZZ &n);
  };

} // ANTL

// Unspecialized template definitions.
#include "../../../src/Exponentiation/ExponentiationYao_impl.hpp"

#endif // EXPONENTIATION_YAO_HPP_H
