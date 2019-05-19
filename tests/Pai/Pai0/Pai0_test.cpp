/**
 * AUTHOR - Parthasarathi Das
 * BRIEF  - Measures the runtimes the Paillier cryptosystem (basic) and writes timing information to files
 */

#include <iostream>
#include <fstream>
#include <math.h>
#include <time.h>
#include <gmp.h>
#include <string.h>
#include <unistd.h>
#include <NTL/ZZ.h>
#include <NTL/RR.h>
#include "../../../include/measurement.h"
#include "../../../include/maxheaders.hpp"
#include "../../../include/Systems/Pai/Pai.hpp"
#include "../../../include/Systems/Pai/Pai0.hpp"
#include <vector>

NTL_CLIENT

int main(int argc, char** argv)
{
    int               mesbits         = strtol(argv[1], NULL, 10);
    int               mbits           = strtol(argv[2], NULL, 10);
    int               switcher        = 0;
    int               parameter_count = 10;
    int               message_count   = 1;
    char              *type           = argv[3];
    char              message_file[PATH_MAX], N_file[PATH_MAX], p_file[PATH_MAX], q_file[PATH_MAX], time_file[PATH_MAX];
    ofstream          timeout;
    RR                factor;
    mpz_t             N, N2, p, q, m, dm, c1, c2, mod, twom, scalar, scalarm;
    gmp_randstate_t   rands;
    PaiPublickey      pk;
    PaiSecretkey      sk;
    PaiPlaintext      pt, *dt, *dtsum, *dtscal;
    PaiCiphertext     *ct, ct1, ct2, *ctsum, *ctscal;
    
    mpz_inits            (N, N2, p, q, m, dm, c1, c2, mod, twom, scalar, scalarm, NULL);
    gmp_randinit_default (rands);
    
    if      (strcmp(type, "crt") == 0) {switcher = 1;}
    else if (strcmp(type, "std") == 0) {}
    else
    {
        cout << "Error: Invalid CRT option\n";
        cout << "Enter 'CRT' for CRT decryption\n";
        exit(0);
    }
    
    Pai0 b(mbits, switcher);
    
    
    // Set the name and path of .txt files for reading and  writing
    FILE *inm;
    FILE *inN;
    FILE *inp;
    FILE *inq;

    snprintf(message_file,  PATH_MAX, "../../../gen/genmod/%s_m.txt", argv[1]);
    snprintf(N_file,        PATH_MAX, "../../../gen/genmod/%s_N.txt", argv[2]);
    snprintf(p_file,        PATH_MAX, "../../../gen/genmod/%s_p.txt", argv[2]);
    snprintf(q_file,        PATH_MAX, "../../../gen/genmod/%s_q.txt", argv[2]);
    snprintf(time_file,     PATH_MAX, "%s_%s_0_%s_time.txt",          argv[1], argv[2], argv[3]);

    
    // Open file for writing
    timeout.open(time_file);
    
    // DS for storing timing information
    Vec<ZZ> ticks;
    Vec<RR> runtime;
    ticks.SetLength(4);
    runtime.SetLength(4);
    
    
    // open message file
    inm  = fopen(message_file,  "r");
    inN  = fopen(N_file,  "r");
    inp  = fopen(p_file,  "r");
    inq  = fopen(q_file,  "r");
    
    // Start
    for (int temp_count = 0; temp_count < parameter_count; temp_count++)
    {
        mpz_inp_str(mod, inN, 16);
        mpz_inp_str(p,   inp, 16);
        mpz_inp_str(q,   inq, 16);
        
        b.init(mod, p, q);
        b.keygen(pk, sk);
        
        // Encrypt and decrpyt for message_count messages
        for (int mcount = 0; mcount < message_count; mcount++)
        {
            // Set plaintext
            mpz_inp_str(m, inm, 16);
            pt.set(m);
            
            // Generate homomorphic parameters
            pk.get(N, N2);
            mpz_mul_ui(twom, m, 2);
            mpz_mod(twom, twom, N);
            mpz_urandomm(scalar, rands, N);
            mpz_mul(scalarm, scalar, m);
            mpz_mod(scalarm, scalarm, N);
            
            // Measure function runtimes
            MEASURE ( ct = &b.encrypt(pt, pk); );
            ticks[0] = ticks[0] + (etime);
            MEASURE ( dt = &b.decrypt(*ct, pk, sk); );
            ticks[1] = ticks[1] + (etime);
            dt->get(dm);
            if (mpz_cmp(m, dm) == 0) { cout << "Success ED\t" << temp_count + 1 << " " << mcount + 1 << endl; }
            else                     { cout << "Failure ED\t" << temp_count + 1 << " " << mcount + 1 << endl; exit(0); }
            
            ct->get(c1);
            ct1.set(c1);
            ct2.set(c1);
            
            MEASURE ( ctsum = &b.evalsum(pk, ct1, ct2); );
            ticks[2] = ticks[2] + (etime);
            dtsum = &b.decrypt(*ctsum, pk, sk);
            dtsum->get(dm);
            mpz_mod(dm, dm, N);
            if (mpz_cmp(twom, dm) == 0) { cout << "Success EDSum\t" << temp_count + 1 << " " << mcount + 1 << endl; }
            else                        { cout << "Failure EDSum\t" << temp_count + 1 << " " << mcount + 1 << endl; exit(0); }
            
            
            MEASURE ( ctscal = &b.evalscal(pk, ct1, scalar); );
            ticks[3] = ticks[3] + (etime);
            dtscal = &b.decrypt(*ctscal, pk, sk);
            dtscal->get(dm);
            mpz_mod(dm, dm, N);
            if (mpz_cmp(scalarm, dm) == 0) { cout << "Success EDSca\t" << temp_count + 1 << " " << mcount + 1 << endl; }
            else                           { cout << "Failure EDSca\t" << temp_count + 1 << " " << mcount + 1 << endl; exit(0); }
            
        } // end-of-message-count-loop
        
        cout << endl;
        
    } // end-of-parameter-for-loop
    
    
    
    // Compute factor before writing timing information to files
    factor = 1000 / (frequency * GHz * parameter_count * message_count);
    
    
    // Write timing information to file
    for (int k = 0; k < 4; k++)
    {
        runtime[k] = conv<RR>(ticks[k]) * factor;
        timeout << " & " << (runtime[k]);
    }
    timeout << "\\\\";
    
    
    timeout.close();
    fclose(inm);
    fclose(inN);
    fclose(inp);
    fclose(inq);
    
    mpz_clears (N, N2, p, q, m, dm, c1, c2, mod, twom, scalar, scalarm, NULL);
    gmp_randclear (rands);
    
    return 0;
}
