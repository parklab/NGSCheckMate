#include "stringhash.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
int main(int argc, char* argv[])
{
    treenode* tname_tree=new_tree();
    char line[100];
    FILE* file;
    int tid=0;

    char* filename = (char*)argv[1];

    file = fopen((const char*)filename,"r");
    /*Check for validity of the file.*/
    if(file == 0)
    {
       printf("can't open file.\n");
       exit(1);
    }
    while(fgets(line,sizeof(line),file)) {
       line[strlen(line)-1]='\0'; 
       fprintf(stderr,"%s\n",line);
       insert_key(line,tid++,tname_tree); 
    }
    fclose(file);

    delete_treenode(tname_tree); tname_tree=NULL;

    return 0;
}

