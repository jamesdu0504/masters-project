# Author: Parthasarathi Das
#
#
# The program reads N / t tables and finds
# 1) the smallest encryption value(s) among E, E2 and E3     and prints the N, t, D and D1 (or D2 or D3) values for each occurence of the smallest value
# 2) the smallest decryption value(s) among D, D1, D2 and D3 and prints the N, t and E     (or E2 or E3) values for each occurence of the smallest value
# 3) the smallest evalsum    value(s) among S, S2 and S3     and prints the N, t and C     (or C2 or C3) values for each occurence of the smallest value
# 4) the smallest evalscal   value(s) among C, C2 and C3     and prints the N, t and S     (or S2 or S3) values for each occurence of the smallest value


import sys




'''
This function takes a tuple as input and returns all the occurences of the smallest element of the tuple and
their indices in the tuple
'''
def smallest_run_times(tuple_list, function_index_list_1, function_index_list_2):
    
    # Create a list of the smallest values of each tuple of each function
    smallest_values_list_1 = []
    smallest_values_list_2 = []
    for i in function_index_list_1:
        smallest_values_list_1.append((min([x for x in tuple_list[i] if x != 0])))
    for i in function_index_list_2:
        smallest_values_list_2.append((min([x for x in tuple_list[i] if x != 0])))



    '''
    # Find and store all the indices of each occurence of the smallest value of each tuple in tuple list
    
    smallest_value_index_list_1 = []
    smallest_value_index_list_2 = []
    for index,value in enumerate(function_index_list_1):
        smallest_value_index_list_1.append([i for i,j in enumerate(tuple_list[value]) if j == smallest_values_list_1[index]])
    for index,value in enumerate(function_index_list_2):
        smallest_value_index_list_2.append([i for i,j in enumerate(tuple_list[value]) if j == smallest_values_list_2[index]])

    print (smallest_value_index_list_1)
    print (smallest_value_index_list_2)
    '''


    # Find and store all the indices of each occurence of the smallest timing value of each function
    val1 = min(smallest_values_list_1)
    val2 = min(smallest_values_list_2)

    function1_tuple_index = []
    function2_tuple_index = []
    
    list1 = []
    list2 = []
    for index,value in enumerate(smallest_values_list_1):
        if (value == val1):
            function1_tuple_index.append(function_index_list_1[index])
            temp = function_index_list_1[index]
            list1.append([i for i,j in enumerate(tuple_list[temp]) if j == val1])

    for index,value in enumerate(smallest_values_list_2):
        if (value == val2):
            function2_tuple_index.append(function_index_list_2[index])
            temp = function_index_list_2[index]
            list2.append([i for i,j in enumerate(tuple_list[temp]) if j == val2])
    '''
    print (smallest_values_list_1)
    print (smallest_values_list_2)
    print ()
    print (list1)
    print (list2)
    print ()
    print (function1_tuple_index)
    print (function2_tuple_index)
    '''
    return list1, list2, function1_tuple_index, function2_tuple_index


def encdecinfo(small_list1, small_list2, funi_list1, funi_list2, crypto_seq, tuplelist):
    print ("The smallest encryption values are:")
    for index, value in enumerate(funi_list1):
        for innervalue in small_list1[index]:
            print (crypto_seq[0], end = " ")
            print (tuplelist[0][innervalue], end = "\t")
            print (crypto_seq[1], end = " ")
            print (tuplelist[1][innervalue], end = "\t")
            
            if (value == 2):
                print (crypto_seq[2], end = " ")
                print (tuplelist[2][innervalue], end = "\t")
                print (crypto_seq[3], end = " ")
                print (tuplelist[3][innervalue], end = "\t")
                print (crypto_seq[4], end = " ")
                print (tuplelist[4][innervalue], end = "\t")
            
            if (value == 5):
                print (crypto_seq[5], end = " ")
                print (tuplelist[5][innervalue], end = "\t")
                print (crypto_seq[6], end = " ")
                print (tuplelist[6][innervalue], end = "\t")
            
            if (value == 7):
                print (crypto_seq[7], end = " ")
                print (tuplelist[7][innervalue], end = "\t")
                print (crypto_seq[8], end = " ")
                print (tuplelist[8][innervalue], end = "\t")
            print ()
        #print ()
    
    
    
    print ("The smallest decryption values are:")
    for index, value in enumerate(funi_list2):
        for innervalue in small_list2[index]:
            print (crypto_seq[0], end = " ")
            print (tuplelist[0][innervalue], end = "\t")
            print (crypto_seq[1], end = " ")
            print (tuplelist[1][innervalue], end = "\t")
            
            if (value == 3):
                print (crypto_seq[2], end = " ")
                print (tuplelist[2][innervalue], end = "\t")
                print (crypto_seq[3], end = " ")
                print (tuplelist[3][innervalue], end = "\t")
            
            if (value == 4):
                print (crypto_seq[2], end = " ")
                print (tuplelist[2][innervalue], end = "\t")
                print (crypto_seq[4], end = " ")
                print (tuplelist[4][innervalue], end = "\t")
            
            if (value == 6):
                print (crypto_seq[5], end = " ")
                print (tuplelist[5][innervalue], end = "\t")
                print (crypto_seq[6], end = " ")
                print (tuplelist[6][innervalue], end = "\t")
            
            if (value == 8):
                print (crypto_seq[7], end = " ")
                print (tuplelist[7][innervalue], end = "\t")
                print (crypto_seq[8], end = " ")
                print (tuplelist[8][innervalue], end = "\t")
            print ()
    #print ()


def homoinfo(small_list1, small_list2, funi_list1, funi_list2, homo_seq, tuplelist):
    print ("The smallest EvalSum values are:")
    for index, value in enumerate(funi_list1):
        for innervalue in small_list1[index]:
            print (homo_seq[0], end = " ")
            print (tuplelist[0][innervalue], end = "\t")
            print (homo_seq[1], end = " ")
            print (tuplelist[1][innervalue], end = "\t")
            
            if (value == 2):
                print (homo_seq[2], end = " ")
                print (tuplelist[2][innervalue], end = "\t")
                print (homo_seq[3], end = " ")
                print (tuplelist[3][innervalue], end = "\t")
            
            
            if (value == 4):
                print (homo_seq[4], end = " ")
                print (tuplelist[4][innervalue], end = "\t")
                print (homo_seq[5], end = " ")
                print (tuplelist[5][innervalue], end = "\t")
            
            if (value == 6):
                print (homo_seq[6], end = " ")
                print (tuplelist[6][innervalue], end = "\t")
                print (homo_seq[7], end = " ")
                print (tuplelist[7][innervalue], end = "\t")
            print ()
        #print ()
    
    
    
    print ("The smallest EvalScal values are:")
    for index, value in enumerate(funi_list2):
        for innervalue in small_list2[index]:
            print (homo_seq[0], end = " ")
            print (tuplelist[0][innervalue], end = "\t")
            print (homo_seq[1], end = " ")
            print (tuplelist[1][innervalue], end = "\t")
            
            if (value == 3):
                print (homo_seq[2], end = " ")
                print (tuplelist[2][innervalue], end = "\t")
                print (homo_seq[3], end = " ")
                print (tuplelist[3][innervalue], end = "\t")
            
            if (value == 5):
                print (homo_seq[4], end = " ")
                print (tuplelist[4][innervalue], end = "\t")
                print (homo_seq[5], end = " ")
                print (tuplelist[5][innervalue], end = "\t")
            
            if (value == 7):
                print (homo_seq[6], end = " ")
                print (tuplelist[6][innervalue], end = "\t")
                print (homo_seq[7], end = " ")
                print (tuplelist[7][innervalue], end = "\t")
            print ()
        #print ()




'''
Each row of the cryptographic timings file contains N t E D D1 E2 D2 E3 D3 values (total 9)
Each row of the homomorphic   timings file contains N t S C    S2 C2 S3 C3 values (total 8)
'''
def main():
    
    scheme_no    = sys.argv[1]
    conductor    = sys.argv[2]
    discriminant = sys.argv[3]
    explen_type  = sys.argv[4]
    
    scheme_list  = ['Basic', 'BasicPlus', 'Variant', 'VariantPlus']
    explen_list  = ['full', 'short']
    crypto_seq   = ['N', 't', 'E ', 'D ', 'D1', 'E2', 'D2', 'E3', 'D3']
    homo_seq     = ['N', 't', 'S ', 'C ', 'S2', 'C2', 'S3', 'C3']


    if (scheme_no == '1'):
        scheme_name = 'Basic'
    elif (scheme_no == '2'):
        scheme_name = 'BasicPlus'
    elif (scheme_no == '3'):
        scheme_name = 'Variant'
    elif (scheme_no == '4'):
        scheme_name = 'VariantPlus'
    else:
        print ("Error: First argument should belong in the range [1, 4]")
        sys.exit()


    if (explen_type == 'f'):
        explen = 'full'
    elif (explen_type == 's'):
        explen = 'short'
    else:
        print ("Error: Fourth argument should be either 'f' for full sized exponents or 's' for short sized exponents")
        sys.exit()


    enc_index_list = [2,5,7]
    dec_index_list = [3,4,6,8]
    sum_index_list = [2,4,6]
    sca_index_list = [3,5,7]
    

    with open("%s/%s/%s/%s_rawedtime.txt" % (scheme_name, conductor, discriminant, explen)) as inp:
        tuplelist = list( zip( *( [int(i) for i in line.strip().split('\t')] for line in inp) ) )
        
    # Set the first two tuple elements of timelist to N-tup and t-tup variables
    Ntup = tuplelist[0]
    ttup = tuplelist[1]
        
    # Find the indices of the smallest running times of encryption and decryption functions
    small_list1, small_list2, funi_list1, funi_list2 = smallest_run_times(tuplelist, enc_index_list, dec_index_list)

    #
    encdecinfo(small_list1, small_list2, funi_list1, funi_list2, crypto_seq, tuplelist)

    with open("%s/%s/%s/%s_rawsctime.txt" % (scheme_name, conductor, discriminant, explen)) as inp:
        tuplelist = list( zip( *( [int(i) for i in line.strip().split('\t')] for line in inp) ) )
    
    # Set the first two tuple elements of timelist to N-tup and t-tup variables
    Ntup = tuplelist[0]
    ttup = tuplelist[1]
    
    # Find the indices of the smallest running times of encryption and decryption functions
    small_list1, small_list2, funi_list1, funi_list2 = smallest_run_times(tuplelist, sum_index_list, sca_index_list)
    
    #
    homoinfo(small_list1, small_list2, funi_list1, funi_list2, homo_seq, tuplelist)




main()
