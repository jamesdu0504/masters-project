/**
 
 @author Parthasarathi Das
 @desc   Generates messages of a given bit length
 
 */


#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <gmp.h>

int main(int argc, char** argv)
{
    // Check input(s)
    if (argc != 2)
    {
        std::cout << "Error: Check the number of input(s)" << std::endl;
        std::cout << "Usage: Message size in bits" << std::endl;
        exit(0);
    }
    
    
    // Declare local variables
    int             mcount;
    int             mbits;
    char            m_file[PATH_MAX];
    mpz_t           m;
    FILE            *outm;
    size_t          buffer;
    gmp_randstate_t rands;
    
    
    // Initialise
    mpz_inits (m, NULL);
    gmp_randinit_default (rands);
    
    mbits  = strtol(argv[1], NULL, 10);
    mcount = 1000;
    
    
    // Set output file
    snprintf(m_file, PATH_MAX, "%d.txt", mbits);
    outm = fopen(m_file, "w");
    
    
    // Generate mcount messages
    for (int i = 0; i < mcount; i++)
    {
        mpz_urandomb(m, rands, mbits);
        buffer = mpz_out_str(outm, 16, m);
        fprintf(outm, "\n");
    }
    
    
    // Close output file
    fclose(outm);
    
    // Clear local variables
    mpz_clears(m, NULL);
    gmp_randclear(rands);
    
    return 0;
}
