#ifndef STRINGHASH2
#define STRINGHASH2

#define BUCKET_EXPANSION_FACTOR 5 //The bucket size will be expanded by this factor and the probability of collision will be reduced by this factor as well.


typedef struct integer2
{
   char c; // 0 or 1
   int i; // index
} int2;


typedef struct hashchainnode
{
   char* key;
   int2 val; // value as an integer pointer
   struct hashchainnode* next;
} node;

typedef struct hashstruct
{
  node** bucket;
  long n;
} hash;


long hashfunc(char*,int,long);
int base2int(char);
hash* newhash(long);
node** create_bucket(long);
long adjust_n(long);
node* create_node(char*, int, int2);
int insert_into_hash(char*, int, int2, hash*);
int2 search_hash(char*, int, hash*);
void delete_hash(hash*);

#endif

