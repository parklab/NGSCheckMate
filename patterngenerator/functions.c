#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include "functions.h"


unsigned char* convert2bi(char* str){
  unsigned char *bstr;
  char bb;
  int i,k,n;

  bstr=calloc(BPATTERNLEN,sizeof(unsigned char));

  n=strlen(str);
  k=0;
  for(i=0;i<n+3;i++){
     if(i>=n) bb=0;
     else if(str[i]=='A') bb=0;
     else if(str[i]=='C') bb=1;
     else if(str[i]=='G') bb=2;
     else bb=3;
     bstr[k]|=bb;
     if(i%4!=3) bstr[k]<<=2; // first three times, skip the fourth one.
     else k++; // only the fourth one.
     if(k>=BPATTERNLEN) break;
  } 
  return(bstr);
}


void convert2str(unsigned char* bstr, char** pstr, char** prcstr){
  unsigned char bb;
  int n,k,i;
  (*pstr)=malloc((PATTERNLEN+1)*sizeof(char));
  (*prcstr)=malloc((PATTERNLEN+1)*sizeof(char));

  n=PATTERNLEN-1;
  k=0;
  for(i=0;i<BPATTERNLEN;i++){
    do{
     bb=bstr[i]&484; //11000000
     bb>>=6;
     bstr[i]<<=2;
     (*pstr)[k]=base[bb];
     (*prcstr)[n-k]=rcbase[bb];
     k++;
    }while(k%4!=0 && k<PATTERNLEN);
  }
  (*pstr)[PATTERNLEN]='\0';
  (*prcstr)[PATTERNLEN]='\0';

}


Patternline parse_patternfileline (char* line){
       char tmpstr[SINGLE_LINE_MAX],pattern[SINGLE_LINE_MAX];
       int numtab,j,i,index;
       Patternline patternline;
       char* pattern_bi;

       numtab=0;j=0;
       for(i=0;i<=strlen(line);i++){
          if(line[i]=='\t'||line[i]=='\0') {
            numtab++;
            tmpstr[j]='\0';
            if(numtab==1) strcpy(pattern,tmpstr);
            else if(numtab==2) index = atoi(tmpstr);
            j=0;
          }
          else { tmpstr[j]=line[i]; j++; }
       }
       if(numtab<3) { fprintf(stderr,"Error: input pattern file must have three columns.\n"); exit(1); }

       pattern_bi=convert2bi(pattern);
       for(i=0;i<BPATTERNLEN;i++)patternline.pattern[i]=pattern_bi[i];
       patternline.index = index;
       free(pattern_bi);

       return(patternline);
}


void convert_pattern_file(char* patternfilename, char* outpatternfilename){
    FILE* pattern_file=0;
    FILE* fout=0;
    char line[SINGLE_LINE_MAX];
    Patternline patternline;
    int k;
    int n; // number of patterns

    // open for the first time to get the line number. 
    pattern_file = fopen((const char*)patternfilename,"r");
    /*Check for validity of the file.*/
    if(pattern_file == 0)
    {
       fprintf(stderr,"can't open pattern file.\n");
       exit(1);
    }

    while (fgets(line, SINGLE_LINE_MAX, pattern_file) != NULL) n++;
    fclose(pattern_file);


    // reopen
    pattern_file = fopen((const char*)patternfilename,"r");
    /*Check for validity of the file.*/
    if(pattern_file == 0)
    {
       fprintf(stderr,"can't open pattern file.\n");
       exit(1);
    }


    fout = fopen((const char*)outpatternfilename , "wb" );
    if(fout==0){
         fprintf(stderr,"can't write to a binary pattern file.\n");
         exit(1);
    }

    fwrite(&n,1,sizeof(int),fout); // the first element is n.

    k=0;
    while (fgets(line, SINGLE_LINE_MAX, pattern_file) != NULL) {
       if(k==0||k==2||k==4){
         patternline = parse_patternfileline(line);
         fwrite(&patternline , 1 , sizeof(Patternline) , fout );
       }
       k++;
       if(k==24) k=0;
    }
    fclose(pattern_file);
    fclose(fout);
}

void printusage_converter(char* program)
{
   printf("Usage : %s input_patternfile(txt) output_patternfile(binary)\n\n",program);
}

void printusage_reader(char* program)
{
   printf("Usage : %s input_patternfile(binary) output_patternfile(txt)\n\n",program);
}


/////////////////////////////////////
// functions for patterntestreader //
/////////////////////////////////////


void convert_pattern_file_back(char* patternfilename, char* outpatternfilename){
    FILE* pattern_file=0;
    FILE* fout=0;
    char line[SINGLE_LINE_MAX];
    Patternline patternline;
    int k,offset,j;
    char original_b;
    char* patternstr=NULL;
    char* rcpatternstr=NULL;
    int n;
 
    // 1. read pattern and decoy files to get the number of keys
    if(strlen(patternfilename)==0){ pattern_file = stdin; }
    else {
      pattern_file = fopen((const char*)patternfilename,"rb");
      /*Check for validity of the file.*/
      if(pattern_file == 0)
      {
         fprintf(stderr,"can't open binary pattern file.\n");
         exit(1);
      }
    }

    fout = fopen((const char*)outpatternfilename , "w" );
    if(fout==0){
         fprintf(stderr,"can't write to a text pattern file.\n");
         exit(1);
    }

    k=0;
    fread(&n,sizeof(int),1,pattern_file); // discard n

    while(fread(&patternline,sizeof(Patternline),1,pattern_file) != 0){
       convert2str(patternline.pattern,&patternstr,&rcpatternstr);
       fprintf(fout,"%s\t%d\t%d\n",patternstr,patternline.index,0);
       fprintf(fout,"%s\t%d\t%d\n",rcpatternstr,patternline.index,0);
       if(k==0) offset=(PATTERNLEN-1)/2;
       else if(k==1) offset=PATTERNLEN-1;
       else offset=0;
       original_b = patternstr[offset];
       for(j=0;j<4;j++){
          if(base[j]!=original_b) {
             patternstr[offset]=base[j];
             rcpatternstr[PATTERNLEN-offset-1]=rcbase[j];
             fprintf(fout,"%s\t%d\t%d\n",patternstr,patternline.index,1);
             fprintf(fout,"%s\t%d\t%d\n",rcpatternstr,patternline.index,1);
          }
       }
       k++;
       if(k==3) k=0;
       if(patternstr!=NULL) free(patternstr); patternstr=NULL;
       if(rcpatternstr!=NULL) free(rcpatternstr); rcpatternstr=NULL;
    }
    fclose(pattern_file);
    fclose(fout);
}


