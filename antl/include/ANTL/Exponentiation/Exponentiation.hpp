/**
 * @file Exponentiation.hpp
 * @author Michael Jacobson
 * @brief virtual superclass for generic exponentiation
 */

#ifndef EXPONENTIATION_H
#define EXPONENTIATION_H

#include "../common.hpp"

namespace ANTL 
{
    /**
     * @brief Virtual superclass for generic exponentiation
     * @remarks This virtual class defines a method for generic exponentiation.
     * Particular exponentiation algorithms are defined and implemented in concrete
     * subclasses.  The base type is templated, and must have as a minimum the following
     * functions defined:
     *   - assign(T &C, const T &A): sets C = A
     *   - mul(T &C, const T &A, const T &B): computes C = AB
     *   - sqr(T &C, const T &A): computes C = AA
     * Concrete subclasses may have further requirements (eg. a cube function) that are
     * described in the documentation of each concrete class.
     */
    
    template <class T>
    class Exponentiation
    {
        public:
        Exponentiation() {};
        virtual ~Exponentiation() {};
        
        /**
         * @brief Computes A^n.  Virtual method, to be overridden by concrete classes corresponding
         * to different exponentiation methods.
         * @param[out] C result of computing A^n
         * @param[in] A base for exponentiation
         * @param[in] n exponent
         *
         * @pre n >= 0
         */
        virtual void power (T &C, const T &A, const ZZ & n) = 0;
    };
} // ANTL

#endif // EXPONENTIATION_H

