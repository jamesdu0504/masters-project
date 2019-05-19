/**
 @author
    Parthasarathi Das
 
 @desc
    Generates conductors f = (p_1...p_N)^t and DK = -(p_1...p_N)q for inputs (N, t, f bit-length, DK bit-length)
 
 @prereqsuisites
    Creat a directory-subdirectory structure where directory name is the desired conductor bit length and contains
    subdirectories where the subdirectory names are desired DK bit lengths.
    Example: Directory name: 80
                Subdirectory 1 name: 1828
                Subdirectory 2 name: 3598
                Subdirectory 3 name: 5972
 
 @note
    When f bit-length is more than DK bit-length, a f bit-length should occur exactly once for some DK bit-length
    Example: Directory 1 name: 3072
                Subdirectory name: 1828
             Directory 2 name: <Anything but 3072>
                Subdirectory name: <Some DK bit-length>
 */



#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <gmp.h>
#include "../../src/measurement.h"
#include "../../src/maxheaders.hpp"



/**
 * For each prime p present (count) in an array of primes (p_array) and a prime q (q)
 *  - Checks if (p/q) = (q/p) = 1 or -1 (legendre_symbol)
 *  - Returns 1 only if all primes p of the array satisfy the Legendre symbol condition with q
 */
int check_legendre(const mpz_t* p_array, const mpz_t q, const int count, const int legendre_symbol)
{
    int j = 0;
    while (j < count)
    {
        if (!((mpz_legendre(p_array[j], q) == legendre_symbol) && (mpz_legendre(q, p_array[j]) == legendre_symbol)))
            return 0;
        j++;
    }
    return 1;
}



/**
 * Generates and stores parameters
 *  - Conductor f = (p1...pN)^t
 *  - Conductor factorisation
 *  - Product of primes in the conductor
 *  - Fundamental discriminant DK
 * for different values of N and t for a conductor size (fbits) at a security level (dkbits)
 */
int gen (const int N, const int t, const int fbits, const int dkbits)
{
    const int POSITIVE_ONE    =   1;
    const int NEGATIVE_ONE    =  -1;
    const int MIN_P_BITS      =  10;
    const int PARAMETER_COUNT =  10;
    const int RAND_BITS       = 100;
    
    int dk_p_count;
    int pbits;
    int lbits;
    int factbits;
    int temp, check, count;
    
    char f_file  [PATH_MAX];
    char fa_file [PATH_MAX];
    char pr_file [PATH_MAX];
    char dk_file [PATH_MAX];
    
    size_t           buffer;
    clock_t          starttime, starttime1;
    mpz_t            p, q, f, dk, seed, product, temp1, temp2, *p_array;
    gmp_randstate_t  rands;
    
    FILE *inpr;
    FILE *inf;
    FILE *infa;
    
    FILE *outpr;
    FILE *outf;
    FILE *outfa;
    FILE *outdk;
    
    
    // Initialise
    if (fbits > dkbits)
    {
        factbits   = dkbits;
        dk_p_count = N;
    }
    else
    {
        factbits   = (int) ceil (((double) fbits) / t);
        dk_p_count = N + 1;
    }
    
    //
    pbits = (int) ceil (((double) factbits) / N);
    
    //
    mpz_inits (p, q, f, dk, seed, product, temp1, temp2, NULL);
    p_array = mpz_init_array (dk_p_count);
    gmp_randinit_default (rands);
    
    
    
    // Exit if pbits < MIN_P_BITS
    if (pbits <= MIN_P_BITS)
        return 1;
    
    
    // Set output file names
    snprintf(f_file,  PATH_MAX, "%d/%d_%d_f.txt",  fbits, N, t);
    snprintf(fa_file, PATH_MAX, "%d/%d_%d_fa.txt", fbits, N, t);
    snprintf(pr_file, PATH_MAX, "%d/%d_%d_pr.txt", fbits, N, t);
    snprintf(dk_file, PATH_MAX, "%d/%d/%d_%d_dk.txt", fbits, dkbits, N, t);
    
    
    // Compute primes of conductor if conductor file does not exist
    // Generate conductors f for the lowest security level and re-use to generate only DK for higher security levels
    if ( !( access ( f_file, F_OK) != -1 ) ) // condition for non-existence
    {
        outf  = fopen (f_file,  "a");
        outfa = fopen (fa_file, "a");
        outpr = fopen (pr_file, "a");
        outdk = fopen (dk_file, "a");
        
        for ( int count = 1; count < PARAMETER_COUNT + 1; count++)
        {
            mpz_rrandomb (seed, rands, RAND_BITS);
            gmp_randseed (rands, seed);
                        
            // Generate the first prime in DK
            mpz_random_prime (rands, p, pbits);
            if (dk_p_count == 1)
            {
                mpz_pow_ui(temp1, p, t);
                while ( !( (mpz_mod_ui(temp2, p, 4) == 3) && (mpz_sizeinbase(temp1, 2) >= fbits) ) )
                {
                    mpz_random_prime (rands, p, pbits);
                    mpz_pow_ui(temp1, p, t);
                }
                mpz_set (dk, p);
            }
            
            // Update the array and the product
            mpz_set (p_array[0], p);
            mpz_set (product, p);
            
            // Generate remaining primes of DK - 1
            for (int i = 1; i < dk_p_count - 1; i++)
            {
                check = 0;
                while (check == 0)
                {
                    mpz_random_prime(rands, p, pbits);
                    check = check_legendre(p_array, p, i, POSITIVE_ONE);
                }
                mpz_set(p_array[i], p);
                mpz_mul(product, product, p);
            }
            
            
            // When fbits < dkbits, check generated f size and regenerate a prime of f if size of f < fbits
            if (fbits < dkbits)
            {
                // Check size first
                mpz_pow_ui(temp1, product, t);
                if (mpz_sizeinbase(temp1, 2) < fbits)
                {
                    mpz_divexact(product, product, p); // remove the last p from product
                    mpz_pow_ui(temp1, product, t);
                    lbits = (int) ceil(((double) fbits - mpz_sizeinbase(temp1, 2)) / t);
                    while (mpz_sizeinbase(temp1, 2) < fbits)
                    {
                        check = 0;
                        while (check == 0)
                        {
                            mpz_random_prime(rands, p, lbits);
                            check = check_legendre(p_array, p, dk_p_count - 2, POSITIVE_ONE);
                        }
                        mpz_mul(temp2, product, p);
                        mpz_pow_ui(temp1, temp2, t);
                    }
                    mpz_set(p_array[dk_p_count - 2], p);
                    mpz_mul(product, product, p);
                }
            }
            
            
            // Generate last prime of DK
            if (dk_p_count > 1)
            {
                lbits = dkbits - mpz_sizeinbase(product, 2) + 1;
                check = 0;
                while (check == 0)
                {
                    mpz_random_prime(rands, p, lbits);
                
                    //Check for fundamental property
                    mpz_mul(temp1, product, p);
                    if ((mpz_mod_ui(temp2, temp1, 4) == 3))
                        check = check_legendre(p_array, p, dk_p_count - 1, NEGATIVE_ONE);
                    else
                        check = 0;
                }
                mpz_set(p_array[dk_p_count - 1], p);
                mpz_mul(dk, product, p);
            }
            
            
            // Set conductor and write product to files
            if (N < dk_p_count)
            {
                mpz_pow_ui(f, product, t);
                buffer = mpz_out_str (outpr, 16, product);
            }
            else
            {
                mpz_pow_ui(f, dk, t);
                buffer = mpz_out_str (outpr, 16, dk);
            }
        
                
            // Write the rest to separate files as well
            for (int i = 0; i < N; i++)
            {
                buffer = mpz_out_str (outfa, 16, p_array[i]);
                fprintf(outfa, "\t");
            }
            buffer = mpz_out_str (outf,  16, f);
            buffer = mpz_out_str (outdk, 16, dk);
                
            fprintf(outf,  "\n");
            fprintf(outfa, "\n");
            fprintf(outpr, "\n");
            fprintf(outdk, "\n");

            // Display information on screen
            std::cout << mpz_sizeinbase(f, 2) << "\t" << mpz_sizeinbase(dk, 2)  << "\t" << N << "\t" << t << "\t" << count << "\n";
        }
        
        // Close files
        fclose(outf);
        fclose(outfa);
        fclose(outpr);
        fclose(outdk);
    }
    
    // ELSE CASE - conductor files exist; just compute DK for higher security levels
    else
    {
        // Open parameter files to read and compute DK
        inf   = fopen (f_file,  "r");
        infa  = fopen (fa_file, "r");
        inpr  = fopen (pr_file, "r");
        outdk = fopen (dk_file, "a");
        
        for ( int count = 1; count < PARAMETER_COUNT + 1; count++ )
        {
            // Read parameters
            mpz_inp_str(f, inf, 16);
            mpz_inp_str(product, inpr, 16);
            for (int i = 0; i < N; i++)
                mpz_inp_str(p_array[i], infa, 16);
            
            // Generate last prime of DK
            lbits = dkbits - mpz_sizeinbase(product, 2) + 1;
            check = 0;
            while (check == 0)
            {
                mpz_random_prime(rands, p, lbits);
                
                //Check for fundamental property
                mpz_mul(temp1, product, p);
                if ((mpz_mod_ui(temp2, temp1, 4) == 3))
                    check = check_legendre(p_array, p, N, NEGATIVE_ONE);
                else
                    check = 0;
            }
            mpz_set(p_array[dk_p_count - 1], p);
            mpz_mul(dk, product, p);
            
            // Write to file
            buffer = mpz_out_str(outdk,  16, dk);
            fprintf(outdk, "\n");
            
            // Display information on screen
            std::cout << mpz_sizeinbase(f, 2) << "\t" << mpz_sizeinbase(dk, 2)  << "\t" << N << "\t" << t << "\t" << count << "\n";
        }
        
        // Close files
        fclose(inf);
        fclose(infa);
        fclose(inpr);
        fclose(outdk);
    }
    return 0;
}



int main(int argc, char** argv)
{
    if (argc != 3)
    {
        std::cout << "Error: Incorrect number of user input(s)" << std::endl;
        std::cout << "Usage: [./t] [Conductor length in bits] [space] [DK length in bits]" << std::endl;
        exit(0);
    }
    
    int buf;
    int fbits  = strtol(argv[1], NULL, 10);
    int dkbits = strtol(argv[2], NULL, 10);
    
    int nmin = 1;
    int nmax = 5;
    
    int tmin = (int) ceil (((double) fbits) / dkbits);
    int tmax = 5;
    if (fbits > dkbits)
        tmax = tmin;
    
    
    // Generate various parameters for each fsize and psize
    for (int i = nmin; i < nmax + 1; i++)
    {
        for (int j = tmin; j < tmax + 1; j++)
        {
            buf = gen(i, j, fbits, dkbits);
            std::cout << std::endl;
        }
    }
    
    return 0;
}
