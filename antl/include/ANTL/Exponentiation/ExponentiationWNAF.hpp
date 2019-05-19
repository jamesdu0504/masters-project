/**
 * @file ExponentiationWNAF.hpp
 * @author Michael Jacobson/Reginald Lybbert
 * @brief class for left-to-right WNAF exponentiation
 */

#ifndef EXPONENTIATION_WNAF_H
#define EXPONENTIATION_WNAF_H

#include <vector>
#include "Exponentiation.hpp"

namespace ANTL
{
    template < class T >
    class ExponentiationWNAF : public Exponentiation<T>
    {
        protected:
        vector<short> e;
        vector<T> precomp;
        vector<T> precompInv;
        short w;
        
        public:
        ExponentiationWNAF() {};
        ~ExponentiationWNAF() {};
        
        void initialize (const T &A, const ZZ &n, const short w);
        void power      (T &C, const T &A, const ZZ &n);
        void clear      ();
    };

} // ANTL

#include "../../../src/Exponentiation/ExponentiationWNAF_impl.hpp"

#endif // EXPONENTIATION_WNAF_H
