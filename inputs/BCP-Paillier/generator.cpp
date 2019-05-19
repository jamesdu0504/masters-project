/**
 * AUTHOR - Parthasarathi Das
 * BRIEF  - Generates semiprimes N = pq for use in BCP and Paillier cryptosystems
 */

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <gmp.h>
#include <time.h>
#include "../include/cheaders.hpp"

/**
 * Function that writes p, q and N = pq to files
 */
void gen(const int msize)
{
    // Local variables
    int              randbits, bits, parameter_count;
    char             p_file[PATH_MAX], q_file[PATH_MAX], N_file[PATH_MAX];
    size_t           buffer;
    mpz_t            p, q, N, seed, temp;
    gmp_randstate_t  rands;
    
    randbits = 100;
    parameter_count = 10;
    mpz_inits(p, q, N, seed, temp, NULL);
    gmp_randinit_default(rands);
    
    // Set output file(s)
    snprintf(p_file, PATH_MAX, "%d_p.txt", msize);
    snprintf(q_file, PATH_MAX, "%d_q.txt", msize);
    snprintf(N_file, PATH_MAX, "%d_N.txt", msize);
    FILE *outp = fopen(p_file,  "w");
    FILE *outq = fopen(q_file,  "w");
    FILE *outN = fopen(N_file,  "w");
    
    // Generate <parameter_count> distinct semiprimes N
    for ( int temp_count = 1; temp_count < parameter_count + 1; temp_count++ )
    {
        
        // Generate a safe prime p
        bits = ceil(msize / 2);
        mpz_random_prime(rands, p, bits);
        while (mpz_mod_ui(temp, p, 12) != 11)
            mpz_random_prime(rands, p, bits);
        
        
        // Generate a safe prime q != p
        bits = msize - mpz_sizeinbase(p, 2);
        mpz_random_prime(rands, q, bits);
        while ( !( (mpz_mod_ui(temp, q, 12) == 11) && (mpz_cmp(p, q) != 0) ) )
            mpz_random_prime(rands, q, bits);
        
        // Set N
        mpz_mul(N, p, q);
        
        // Write to files
        buffer = mpz_out_str(outp,  16, p);
        buffer = mpz_out_str(outq,  16, q);
        buffer = mpz_out_str(outN,  16, N);
        fprintf(outp, "\n");
        fprintf(outq, "\n");
        fprintf(outN, "\n");
        
        // Seed the random state
        mpz_rrandomb(seed, rands, randbits);
        gmp_randseed(rands, seed);
        
        // Display run information
        std::cout << msize << "\t" << mpz_sizeinbase(N, 2) << "\t" << temp_count << std::endl;
    }
    
    // Close parameter files
    fclose(outp);
    fclose(outq);
    fclose(outN);
    
    // Clear variables
    mpz_clears(p, q, N, seed, temp, NULL);
    gmp_randclear(rands);
}

int main(int argc, char** argv)
{
    // Check input(s)
    if (argc != 2)
    {
        std::cout << "Error: Check the number of input(s)" << std::endl;
        std::cout << "Usage: RSA modulus size in bits" << std::endl;
        exit(0);
    }
    
    int msize = strtol(argv[1], NULL, 10);
    
    // Generate N
    gen(msize);
    
    return 0;
}
