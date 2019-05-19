/**
 * @file MultiExponentiationPippenger_impl.hpp
 * @author Reginald Lybbert
 * @brief class for Pippenger's Multi-exponentiation Algorithm, as described in
 * "Pippenger's Multiproduct and Multiexponentiation Algorithms"; Ryan Henry; 2010
 * 
 *   cacr.uwaterloo.ca/techreports/2010/cacr2010-26.pdf
 */

#include<algorithm>
#include<math.h>

using namespace ANTL;



template< class T >
void
MultiExponentiationPippenger<T>::decompose( vector<T> &x_prime, vector<vector<long> > &y_prime,
                                              const vector<T> &x, const vector<vector<ZZ> > &y, const long r, const long b ){
    x_prime.clear();
    x_prime.resize(r*x.size());
    for(long i = 0; i < x.size(); i++){
        for(long j = 0; j < r; j++){
            if(j == 0){
                assign(x_prime[r*i],x[i]);
            }else{
                assign(x_prime[r*i + j],x_prime[r*i + j - 1]);
                for(long k = 0; k < b; k++){
                    sqr(x_prime[r*i + j],x_prime[r*i+j]);
                }
            }
        }
    }

    y_prime.clear();
    y_prime.resize(b*y.size());
    for(long i = 0; i < y.size(); i++){
        for(long j = 0; j < y[i].size(); j++){     //we require y[i].size() = x.size()
            for(long k = 0; k < r; k++){
                for(long l = 0; l < b; l++){
                    if(bit(y[i][j],k*b + l) == 1){
                       y_prime[i*b + (b - l - 1)].push_back(j*r+k);
                    }
                }
            }
        }
    }

}

template< class T >
void
MultiExponentiationPippenger<T>::combine( vector<T> &C, const vector<T> &x, const long &r){

        C.clear();
        C.resize(x.size()/r);
        long q = -1;
	for(long i = 0; i < x.size(); i++){
            if(i%r == 0){
               q++;
               assign(C[q],x[i]);
            }else{
               sqr(C[q],C[q]);
               mul(C[q],C[q],x[i]);
            }
        }
}



template< class T >
void
MultiExponentiationPippenger<T>::multiprod( vector<T> &C, const vector<T> &x, const vector<vector<long> > &y){

  long ell, c;
  vector<long> alpha,beta;
  getParams(ell,c,alpha,beta,x,y);
  computeMultiProd(C,0,ell,c,alpha,beta,x,y);
}


template< class T >
void
MultiExponentiationPippenger<T>::getParams( long &ell, long &c, vector<long> &alpha, vector<long> &beta, 
                                              const vector<T> &x, const vector<vector<long> > &y){
    
    ell = 4;       //temporary values...
    c = 3;         //TODO: Fix these parameters
    alpha.clear();
    alpha.push_back(20);
    alpha.push_back(10);
    alpha.push_back(6);
    alpha.push_back(4);
    beta.clear();
    beta.push_back(2);
    beta.push_back(2);
    beta.push_back(2);
    beta.push_back(2);
    if(y.empty()){
       ell = 0;
       return;
    }

}


template< class T >
void
MultiExponentiationPippenger<T>::computeMultiProd( vector<T> &C, const long i, const long &ell, const long &c, const vector<long> &alpha, 
                                                    const vector<long> &beta, const vector<T> &x, const vector<vector<long> > &y){

  vector<T> x_prime,x_doubleprime;
  vector<vector<long> > y_prime,y_doubleprime;
  if(i == ell){
     naiveMultiply(C,x,y);
  }else if(i == 0){
     inputPartition(x_prime,y_prime,x,y,c);
     computeMultiProd(C,1,ell,c,alpha,beta,x_prime,y_prime);
  }else if(i%2 == 1){
     outputClump(y_prime,y_doubleprime,x,y,alpha[i],beta[i]);
     computeMultiProd(x_doubleprime,i+1,ell,c,alpha,beta,x,y_prime);
     x_doubleprime.insert(x_doubleprime.begin(),x.begin(),x.end());
     computeMultiProd(C,i+1,ell,c,alpha,beta,x_doubleprime, y_doubleprime);
  }else{
     inputClump(x_prime,y_prime,x,y,alpha[i],beta[i]);
     computeMultiProd(C,i+1,ell,c,alpha,beta,x_prime,y_prime);
  }
}


template< class T >
void
MultiExponentiationPippenger<T>::inputPartition(vector<T> &x_prime, vector<vector<long> > &y_prime, const vector<T> &x,
                                                      const vector<vector<long> > &y, const long &c){
   x_prime.clear();
   y_prime.clear();
   y_prime.resize(y.size());
   long numOfPartitions = (long) ceil(((double) x.size())/c);
   vector<T> P_i, X_i;
   vector<vector<long> > Y_i;
   long new_index = 0;
   long indices[(1 << c)];   //think of this as a hash table under the perfect hash: [a,b,c] goes to 2^a+2^b+2^c
   long indicesIndex;
   vector<long> y_test;
   for(int i = 0; i < numOfPartitions; i++){
      fill_n(indices, (1<<c), -1);
      P_i.clear();
      X_i.clear();
      Y_i.clear();
      for(int ic = i*c; ic < (i+1)*c && ic < x.size(); ic++){
          P_i.push_back(x[ic]);
      }
      for(int y_index = 0; y_index < y.size(); y_index++){
         indicesIndex = 0;
         y_test.clear();
         for(int inner_y_index = 0; inner_y_index < y[y_index].size(); inner_y_index++){ 
             if(y[y_index][inner_y_index] >= i*c && y[y_index][inner_y_index] < (i+1)*c){
                 y_test.push_back(y[y_index][inner_y_index] - i*c);
                 indicesIndex += (1 << y_test.back());
             }
         }
         if(indices[indicesIndex] == -1){
             Y_i.push_back(y_test);
             indices[indicesIndex] = new_index;
             new_index++;
         }
         y_prime[y_index].push_back(indices[indicesIndex]);          
      }

      naiveMultiply(X_i,P_i,Y_i);
    
      for(long t = 0; t < X_i.size(); t++){ 
          x_prime.push_back(X_i[t]);
      }
   }

}

//This algorithm is broken...
template< class T >
void
MultiExponentiationPippenger<T>::outputClump(vector<vector<long> > &y_prime, vector<vector<long> > &y_doubleprime, const vector<T> &x, 
                                                       const vector<vector<long> > &y, const long alpha, const long beta){

   long numOfPartitions = (long) ceil(((double) y.size())/alpha);
   long clumpCounter; //counts up to beta to find clumps
   vector<long>::const_iterator it;  //temporary iterator
   long indicesIndex;
   vector<long> singles;
   long currIndex = 0;
   long indices[1<<alpha]; 

   y_doubleprime.clear();
   y_doubleprime.resize(y.size());



   for(long i = 0; i < numOfPartitions; i++){
       fill_n(indices,(1<<alpha),-1);
       for(int x_index = 0; x_index < x.size(); x_index++){
           clumpCounter = 0;
           indicesIndex = 0;
           for(int y_index = i*alpha; y_index < (i+1)*alpha && y_index < y.size(); y_index++){  //look through partition i
               

               it = find(y[y_index].begin(),y[y_index].end(),x_index);

               if(it != y[y_index].end()){
                    clumpCounter++;
                    indicesIndex += (1 << (y_index - i*alpha));
                    singles.push_back(y_index);

                    if(clumpCounter == beta){
                        clumpCounter = 0;

                        if(indices[indicesIndex] == -1){
                            indices[indicesIndex] = currIndex;
                            currIndex++;
                            y_prime.resize(currIndex);

                            for(long new_index = 0; new_index < singles.size(); new_index++){
                                y_doubleprime[singles[new_index]].push_back(currIndex + x.size() - 1);
                            }
                        }

                        singles.clear();

                        y_prime[indices[indicesIndex]].push_back(x_index); 
                        indicesIndex = 0;                    
                    }
               }
            }

            for(long new_index = 0; new_index < singles.size(); new_index++){
                y_doubleprime[singles[new_index]].push_back(x_index);
            }
            singles.clear();

       } 
    }
}


template< class T >
void
MultiExponentiationPippenger<T>::inputClump(vector<T> &x_prime, vector<vector<long> > &y_prime, const vector<T> &x,
                                                  const vector<vector<long> > &y, const long alpha, const long beta){

    long currIndex = 0;  //The next index in x_prime in which to place a clump or single
    vector<long> clumpIndices;  //The indices in X_prime that hold clumps, as opposed to singles    
    
    vector<long> currClump;	//where to store the current clump
    currClump.resize(beta);
    long currClumpIndex;    

    long clumpHash[(1<<alpha)];   //keep track of which clumps (and singles) have been defined, and where.
    long clumpHashIndex = 0;

    vector<vector<long> > clumps;  // a place to store the clumps, and pass in to compute them.
    vector<T> clumpComputed;    // a place to store the computed values of the clumps  

    x_prime.clear();
    y_prime.clear();
    y_prime.resize(y.size());     //Number of outputs stays constant


    for(long i = 0; i < x.size(); i += alpha){   //for each partition
        fill_n(clumpHash, (1<<alpha), -1);       //clean out the hash 
        for(long j = 0; j < y.size(); j++){      //for each output
            currClumpIndex = 0;			 //restart counting clumps
            clumpHashIndex = 0;
            for(long k = 0; k < y[j].size(); k++){  //for each index in the output
                if((0 <= y[j][k] - i) && (y[j][k] - i) < alpha){  //if the index is in the partition
                      currClump[currClumpIndex] = y[j][k];      //add the index to the current clump
                      clumpHashIndex += (1 << (y[j][k]-i));
                      currClumpIndex++;
                      if(currClumpIndex == beta){           //if the clump is full
                         if(clumpHash[clumpHashIndex] == -1){   //if the clump has not been found elsewhere
                            clumpHash[clumpHashIndex] = currIndex;
                            clumpIndices.push_back(currIndex);	   //Save a place in x_prime for the clump
                            clumps.push_back(currClump);           //Put in clumps, so it is ready to be computed
                            x_prime.push_back(x[0]);		   //garbage placeholder
                            currIndex++;
                         }
                         y_prime[j].push_back(clumpHash[clumpHashIndex]);  //store clumps index in x_prime in the appropriate output
                         currClumpIndex = 0;				   //restart building the clump
                         clumpHashIndex = 0;
                      }
                }
            }
            for(long l = 0; l < currClumpIndex; l++){  //for everything in the incomplete clump
                if(clumpHash[(1<<(currClump[l] - i))] == -1){ //if the single has not been found elsewhere
                   clumpHash[(1<<(currClump[l] - i))] = currIndex;
                   x_prime.push_back(x[currClump[l]]);		//put the single in x_prime 
                   currIndex++;
                }
                y_prime[j].push_back(clumpHash[(1<<(currClump[l] - i))]);   //store the single's index in the appropriate output
            }
        }
    }


    naiveMultiply(clumpComputed, x, clumps);    //Compute the values of the clumps
   
    for(int i = 0; i < clumpComputed.size(); i++){
        x_prime[clumpIndices[i]] = clumpComputed[i];   //Put the computed values in the stored spaces in x_prime
    }
}



template< class T >
void
MultiExponentiationPippenger<T>::naiveMultiply(vector<T> &C, const vector<T> &x, const vector<vector<long> > &y){

  C.clear();
  T y_i;
  for(vector<vector<long> >::const_iterator S_i = y.begin(); S_i != y.end(); S_i++){
      id(y_i);
      for(vector<long>::const_iterator x_j = (*S_i).begin(); x_j != (*S_i).end(); x_j++){
          mul(y_i,y_i,x[(*x_j)]);
      } 
      C.push_back(y_i);
  }
}


template< class T >
void
MultiExponentiationPippenger<T>::power(vector<T> &C, const vector<T> &A, const vector<vector<ZZ> > &n){
 
  ZZ k;
  ZZ e;
  k = 0;
  for(vector<vector<ZZ> >::const_iterator y_i = n.begin(); y_i != n.end(); y_i++){
      e = *max_element((*y_i).begin(), (*y_i).end());
      if(e > k){
          k = e;
      }
  }


  long a,b;
  a = (long) ceil(sqrt( ((double) A.size())*NumBits(k+1)/n.size()));
  b = (long) ceil(sqrt( ((double) n.size())*NumBits(k+1)/A.size()));
  vector<T> x_prime, x_doubleprime;
  vector<vector<long> > y_prime,y_doubleprime;
  if(n.size() >= A.size()){
      decompose(x_prime,y_prime,A,n,a,b);
      multiprod(x_doubleprime,x_prime,y_prime);
      combine(C,x_doubleprime,a);
  }else{
      decompose(x_prime,y_prime,A,n,b,a);
      multiprod(x_doubleprime,x_prime,y_prime);
      combine(C,x_doubleprime,b);
  }
}
