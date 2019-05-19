#include <unistd.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <time.h>
#include <gmp.h>
#include <mpfr.h>
#include <iostream>
#include <fstream>
#include <vector>


#include "../../../src/measurement.h"
#include "../../../src/maxheaders.hpp"
#include "../../../src/QF.hpp"
#include "../../../src/systems/CL/Basic.hpp"


int main(int argc, char** argv)
{
    int               fbits           = strtol(argv[1], NULL, 10);
    int               dkbits          = strtol(argv[2], NULL, 10);
    int               ebits           = 0;
    int               fmaxbits        = 0;
    int               parameter_count = 10;
    int               message_count   = 10;
    int               count           = 0;
    int               mcounter        = 0;
    int               window_size     = 8;
    int               N               = 5;
    int               t               = 5;
    char              *type           = argv[3];
    
    char              infile_m  [PATH_MAX];
    char              infile_f  [PATH_MAX];
    char              infile_fa [PATH_MAX];
    char              infile_pr [PATH_MAX];
    char              infile_dk [PATH_MAX];
    
    char              outfile_ed  [PATH_MAX];
    char              outfile_sc  [PATH_MAX];
    char              outfile_red [PATH_MAX];
    char              outfile_rsc [PATH_MAX];
    
    ofstream          edout;
    ofstream          scout;
    ofstream          redout;
    ofstream          rscout;
    
    CLPublickey       clpk;
    CLSecretkey       clsk;
    CLPlaintext       clpt, *cldt, *cldtsum, *cldtscal;
    CLCiphertext      clct1, clct2, *clct, *clctsum, *clctscal;
    
    mpz_t             m, dm, twom, alpha, alpham, c1, c2, c3, c4, DK, prod, con, gcd, seed, temp, *pN;
    RR                factor;
    QF                CC1, CC2;
    vector<QF>        CC1N, CC2N;
    mpz_qform_t       form1;
    gmp_randstate_t   rands;
    mpz_qform_group_t group1;
    
    
    // Initialise variables
    if (strcmp(type, "short") == 0)
    {
        switch(dkbits)
        {
            case(1828):
                ebits = 2 * 128;
                break;
            case(3598):
                ebits = 2 * 192;
                break;
            case(5972):
                ebits = 2 * 256;
                break;

            case(3072):
                ebits = 2 * 128;
                break;
            case(7680):
                ebits = 2 * 192;
                break;
            case(15360):
                ebits = 2 * 256;
                break;

            default:
                cout << "Error: Check inputs\n";
                exit(0);
                break;
        }
    }
    else if (strcmp(type, "full") == 0) {}
    else
    {
        cout << "Error: Incorrect exponent size argument\n";
        exit(0);
    }

    //
    Basic obj(window_size, fbits, dkbits, ebits);
    
    //
    mpz_inits            (m, dm, DK, prod, con, gcd, seed, temp, NULL);
    mpz_inits            (c1, c2, c3, c4, twom, alpha, alpham, NULL);
    
    mpz_qform_init       (&group1, &form1);
    mpz_qform_group_init (&group1);
    
    CC1.init             (group1);
    CC2.init             (group1);
    
    clct1.init           (group1);
    clct2.init           (group1);
    
    gmp_randinit_default (rands);
    
    pN = mpz_init_array(N);
    
    CC1N.resize(N);
    CC2N.resize(N);
    
    for (int i = 0; i < N; i++)
    {
        CC1N[i].init (group1);
        CC2N[i].init (group1);
    }
    
    // Set fmaxbits
    fmaxbits = (dkbits / 2) - 2;
    
    // Set the name and path of .txt files to write in timing information
    snprintf (outfile_ed,  PATH_MAX, "%s/%s/%s_ed.txt",  argv[1], argv[2], argv[3]);
    snprintf (outfile_sc,  PATH_MAX, "%s/%s/%s_sc.txt",  argv[1], argv[2], argv[3]);
    snprintf (outfile_red, PATH_MAX, "%s/%s/%s_red.txt", argv[1], argv[2], argv[3]);
    snprintf (outfile_rsc, PATH_MAX, "%s/%s/%s_rsc.txt", argv[1], argv[2], argv[3]);
    
    edout  .open (outfile_ed);
    scout  .open (outfile_sc);
    redout .open (outfile_red);
    rscout .open (outfile_rsc);
    
    
    
    //
    for (int i = 1; i < N + 1; i++)
    {
        for (int j = 1; j < t + 1; j++)
        {
            snprintf (infile_m,  PATH_MAX, "../../../inputs/messages/%s.txt",       argv[1]               );
            snprintf (infile_f,  PATH_MAX, "../../../inputs/CL/%s/%d_%d_f.txt" ,    argv[1],          i, j);
            snprintf (infile_fa, PATH_MAX, "../../../inputs/CL/%s/%d_%d_fa.txt",    argv[1],          i, j);
            snprintf (infile_pr, PATH_MAX, "../../../inputs/CL/%s/%d_%d_pr.txt",    argv[1],          i, j);
            snprintf (infile_dk, PATH_MAX, "../../../inputs/CL/%s/%s/%d_%d_dk.txt", argv[1], argv[2], i, j);
            
            //ed time: E D D1 E2 D2 E3 D3
            //sc time: S C    S2 C2 S3 C3
            
            Vec<ZZ> edticks;
            Vec<ZZ> scticks;
            Vec<RR> edtime;
            Vec<RR> sctime;
            
            edticks .SetLength (7);
            edtime  .SetLength (7);
            scticks .SetLength (6);
            sctime  .SetLength (6);
            

            // if file exists, begin timing computation
            if ( access (infile_f, F_OK) != -1 )
            {
                FILE *inf  = fopen(infile_f,  "r");
                FILE *infa = fopen(infile_fa, "r");
                FILE *inpr = fopen(infile_pr, "r");
                FILE *indk = fopen(infile_dk, "r");
                FILE *inm  = fopen(infile_m,  "r");
                
                count = 0;
                for (int temp_count = 0; temp_count < parameter_count; temp_count++)
                {
                    // read in product of primes of f, f, DK
                    mpz_inp_str(prod, inpr, 16);
                    mpz_inp_str(con, inf, 16);
                    mpz_inp_str(DK, indk, 16);
                    
                    // read conductor factorisation file
                    for (int k = 0; k < i; k++)
                        mpz_inp_str(pN[k], infa, 16);
                    
                    // send parameters to keygen
                    obj.initialise (i, j, con, pN, DK);
                    obj.keygen (clpk, clsk);
                    
                    
                    // encrypt and decrpyt each message
                    mcounter = 0;
                    for (int mcount = 0; mcount < message_count; mcount++)
                    {
                        // Get message
                        mpz_inp_str(m, inm, 16);
                        
                        // Check gcd
                        mpz_gcd(gcd, m, con);
                        if ((mpz_cmp_ui(gcd, 1) == 0) && mpz_cmp(m, con) < 0)
                        {
                            mcounter++;
                            clpt.set(m);
                            
                            if (fbits <= fmaxbits)
                            {
                                // 1. Modular inverse decryption
                                MEASURE ( clct = &obj.encrypt(clpt, clpk); );
                                edticks[0] = edticks[0] + (etime);
                                MEASURE ( cldt = &obj.decrypt(*clct, clpk, clsk); );
                                edticks[1] = edticks[1] + (etime);
                                cldt->get(dm);
                                assert (mpz_cmp(m, dm) == 0);
                                
                                if (i > 1)
                                {
                                    // 2. CRT 1 decryption
                                    MEASURE ( cldt = &obj.decryptCRT(*clct, clpk, clsk); );
                                    edticks[2] = edticks[2] + (etime);
                                    cldt->get(dm);
                                    assert (mpz_cmp(m, dm) == 0);
                                    
                                    
                                    // 3. CRT 2 encryption and decryption
                                    MEASURE ( clct = &obj.encrypt2(clpt, clpk); );
                                    edticks[3] = edticks[3] + (etime);
                                    MEASURE ( cldt = &obj.decrypt2(*clct, clpk, clsk); );
                                    edticks[4] = edticks[4] + (etime);
                                    cldt->get(dm);
                                    assert (mpz_cmp(m, dm) == 0);
                                    
                                    
                                    // 4. CRT 3 encryption and decryption
                                    MEASURE ( clct = &obj.encrypt3(clpt, clpk); );
                                    edticks[5] = edticks[5] + (etime);
                                    MEASURE ( cldt = &obj.decrypt3(*clct, clpk, clsk); );
                                    edticks[6] = edticks[6] + (etime);
                                    cldt->get(dm);
                                    assert (mpz_cmp(m, dm) == 0);
                                }
                            }
                            
                            if (fbits > dkbits)
                            {
                                // 1. Pohlig Hellman decryption
                                MEASURE ( clct = &obj.djsencrypt(clpt, clpk); );
                                edticks[0] = edticks[0] + (etime);
                                MEASURE ( cldt = &obj.djsdecrypt(*clct, clpk, clsk); );
                                edticks[1] = edticks[1] + (etime);
                                cldt->get(dm);
                                assert (mpz_cmp(m, dm) == 0);
                                
                                
                                if (i > 1)
                                {
                                    // 2. CRT 3 Encryption and Decryption
                                    MEASURE ( clct = &obj.encrypt3(clpt, clpk); );
                                    edticks[5] = edticks[5] + (etime);
                                    MEASURE ( cldt = &obj.decrypt3(*clct, clpk, clsk); );
                                    edticks[6] = edticks[6] + (etime);
                                    cldt->get(dm);
                                    assert (mpz_cmp(m, dm) == 0);
                                }
                            }
                            
                        } // gcd and m size check
                        
                    } // end of message loop
                    
                    count = count + mcounter;
                    
                } // end of parameter loop
                
                
                // Close parameter files
                fclose(inf);
                fclose(infa);
                fclose(inpr);
                fclose(indk);
                fclose(inm);
                
                // Print timing information on screen
                cout << count << "\tmessages checked for " << i << " " << j << endl;

                // Compute factor before writing timing information to files
                factor = 1000 / (frequency * GHz * count);
                
                // Write raw E D timing information to file
                redout << i << "\t" ;
                redout << j << "\t" ;
                for (int k = 0; k < 7; k++)
                {
                    edtime[k] = conv<RR>(edticks[k]) * factor;
                    redout << round(edtime[k]) << "\t";
                }
                redout << endl;
                
            } // end of file check condition

        } //j
    
    }//i

    
    // Close time files
    redout.close();
    
    // Deallocate memory
    mpz_clears            (m, dm, DK, prod, con, seed, gcd, temp, NULL);
    mpz_clears            (c1, c2, c3, c4, twom, alpha, alpham, NULL);
    mpz_clear_array       (pN, N);
    gmp_randclear         (rands);
    mpz_qform_clear       (&group1, &form1);
    mpz_qform_group_clear (&group1);
    
    return 0;
}
