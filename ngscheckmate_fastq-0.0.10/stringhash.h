
/* This stringhash is a simple tree structure where each node corresponds to a letter of a key string and a value is a nonnegative integer (-1 is used for errors, so the stored value must be nonnegative). It does not support individual key deletion. Over-writing an existing key is supported. */

#ifndef INTEGER2
#define INTEGER2
typedef struct integer2
{
   char c; // 0 or 1
   int i; // index
} int2;
#endif

#ifndef HASHTREENODE
#define HASHTREENODE
typedef struct hashtreenode
{
   char letter;
   int2 val;
   int size;
   struct hashtreenode** next;
} treenode;
#endif

treenode* insert_key(char*,int2,treenode*);
int2 search_treehash(char*,treenode*,int);
treenode* new_treenode(char);
treenode* new_lasttreenode(char,int2);
treenode* new_tree(void);
treenode* new_branch(treenode*);
void delete_treenodeptrarray(treenode**, int);
void delete_treenode(treenode*);

