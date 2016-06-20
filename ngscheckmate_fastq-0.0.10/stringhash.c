#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "stringhash.h"



treenode* new_treenode(char c)
{
   treenode *newnode;
   newnode = (treenode*)malloc(sizeof(treenode));
   newnode->letter=c;
   newnode->size=0;
   newnode->val.c=0;
   newnode->val.i=-1;
   newnode->next=NULL;
   //fprintf(stdout,"node created with letter %c\n",c);  /* DEBUGGING */
   return(newnode);
}

treenode* new_lasttreenode(char c,int2 val)
{
   treenode *newnode;
   newnode = (treenode*)malloc(sizeof(treenode));
   newnode->letter=c;
   newnode->size=0;
   newnode->val=val;
   newnode->next=NULL;
   //fprintf(stdout,"last node created with letter %c and value %d\n",c,val);  /* DEBUGGING */
   return(newnode);
}


treenode* new_tree(void)
{
   treenode *newnode;
   newnode = (treenode*)malloc(sizeof(treenode));
   newnode->letter=0;  //NULL letter for the root
   newnode->size=0;
   newnode->val.c=0;
   newnode->val.i=-1;
   newnode->next=NULL;
   //fprintf(stdout,"new tree created.\n");  /* DEBUGGING */
   return(newnode);
}

treenode* new_branch(treenode* node)
{
   if(node->next==NULL) node->next = (treenode**)calloc(1,sizeof(treenode*));
   else { 
       node->next = realloc(node->next,(node->size+1)*sizeof(treenode*));
       node->next[node->size]=NULL;
   }
   node->size++;
   return(node);
}

treenode* insert_key(char* s,int2 val,treenode* node)
{
  int i;
  
  if(strlen(s)==1) { // reached last character
    if(node->size>0){
      for(i=0;i<node->size;i++){
         if(s[0]==node->next[i]->letter) { node->next[i]->val=val; return(node);}
      }
    }
    //either node size==0 or new key letter
    node = new_branch(node);
    node->next[node->size-1]=new_lasttreenode(s[0],val); 
    return(node); 
  }
  else { // not the last character
    if(node->size>0){
      for(i=0;i<node->size;i++){
         if(s[0]==node->next[i]->letter) {  node->next[i] = insert_key(s+1,val,node->next[i]); return(node); }
      }
    }
    //either node size==0 or new key letter
    node = new_branch(node); 
    node->next[node->size-1]=new_treenode(s[0]);  
    node->next[node->size-1]=insert_key(s+1,val,node->next[node->size-1]);
    return(node);  
  }
}



void delete_treenodeptrarray(treenode** arr, int size){
   int j;
   if(arr!=NULL){
    for(j=0;j<size;j++) delete_treenode(arr[j]);
    free(arr);
   }
}

void delete_treenode(treenode* node){
   if(node!=NULL){
     delete_treenodeptrarray(node->next,node->size);
     free(node);
   }
}



/* returns value associated with a given key string.
Returns -1 if the key doesn't exists. */
int2 search_treehash(char* s,treenode* node,int k){
    int j;
    int2 val,errval;
    errval.c=0; errval.i=-1;

    if(node==NULL) return(errval);
    if(strlen(s)==0 || k==0) return(errval); /* length 0 string cannot be stored in this structure */
    if(strlen(s)==1 || k==1) {
       for(j=0;j<node->size;j++) {
          if(node->next[j]->letter==s[0]) return(node->next[j]->val); /* this is the value we're looking for */ 
       }
    }
    else {
       for(j=0;j<node->size;j++) {
          if(node->next[j]->letter==s[0]) return( search_treehash(s+1,node->next[j],--k) ); 
       }
    }
    return(errval); /* the key string doesn't exist. */
}


