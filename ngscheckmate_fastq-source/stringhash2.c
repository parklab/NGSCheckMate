#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "stringhash2.h"

// returns a hash value (bucket index) given a string and length (k) and number of keys (n)
// returns -1 if string is shorter than the specified length
long hashfunc(char* str,int k,long n){
   int i;
   long hashval=0;
   if(strlen(str)<k) return(-1); // invalid
   for(i=0;i<k;i++) hashval = hashval*5 + base2int(str[i]);
   return(hashval%n);
}

int base2int(char b){
  switch(b){
    case 'A':
      return(0);
    case 'C':
      return(1);
    case 'G':
      return(2);
    case 'T':
      return(3);
    default:
      return(4);
  }
}


hash* newhash(long n){
  hash* h;
  h=(hash*)malloc(1*sizeof(hash));
  //fprintf(stderr,"hashbucket size=%d (before adjustment)\n",n); // DEBUGGING
  n=adjust_n(n);
  //fprintf(stderr,"hashbucket size=%d (after adjustment)\n",n); // DEBUGGING
  h->bucket = create_bucket(n);
  h->n=n;
  return(h);
}

// create a bucket array given number of keys (n) and initialize to NULL.
node** create_bucket(long n){
   return((node**)calloc(n,sizeof(node*)));
}


// find the minimal prime number that is larger than or equal to n. (to minimize collision)
long adjust_n(long n){
  long i;
  n*=BUCKET_EXPANSION_FACTOR;
  for(i=2;i<n;i++){
    if(n%i==0) { n++; i=1; }
  }
  return(n);
}


// create a node and initialize it.
node* create_node(char* key, int k, int2 val){
   node* p = (node*)malloc(sizeof(node));
   p->key=(char*)malloc((k+1)*sizeof(char));
   strncpy(p->key,key,k);
   p->val=val;
   p->next=NULL;
   return(p);
}


// return 0 if successful.
// return 1 if not.
int insert_into_hash(char* key, int k, int2 val, hash* h){
   int index;
   node* p;
   if(h->bucket!=NULL){
      index = hashfunc(key,k,h->n);
      if(h->bucket[index]==NULL) h->bucket[index]=create_node(key,k,val); 
      else {  // collision.
         //fprintf(stderr,"collision.\n"); // DEBUGGING
         p=h->bucket[index];
         while(p->next!=NULL) p=p->next;
         p->next=create_node(key,k,val);
      }
      return(0);
   }
   else return(1);
}


// return an int2 variable with c=0 and i=-1 if either the bucket is not allocated or if the key doesn't exist.
// otherwise, returns the value associated with the key
int2 search_hash(char* key, int k, hash* h){
  long index;
  node* p;
  int2 err;
  err.c=0; err.i=-1;
  if(h->bucket==NULL) return(err);
  else {
    index=hashfunc(key,k,h->n);
    p=h->bucket[index];
    while(p!=NULL && strncmp(p->key,key,k)!=0) p=p->next;  // even if there is one entry associated with the index, we still need to check the key because the query string may happen to have the same hash modulo.
    if(p!=NULL) return(p->val);
    else return(err); // key doesn't exist
  }
}


// delete hash
void delete_hash(hash* h){
  int i;
  node *p, *q;
  for(i=0;i<h->n;i++){
    p=h->bucket[i];
    while(p!=NULL) { q=p->next; free(p->key); free(p); p=q; }
  }
  free(h->bucket); h->bucket=NULL;
  free(h);
}



