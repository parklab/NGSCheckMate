#include "ngscheckmate.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <zlib.h>
#include "stringhash2.h"
#include "patternconvert.h"

pthread_mutex_t mutex1 = PTHREAD_MUTEX_INITIALIZER;


hash* read_patternfile_construct_hash(char* patternfilename){
    FILE* pattern_file=0;
    FILE* fout=0;
    char line[SINGLE_LINE_MAX];
    Patternline patternline;
    int k,offset,j;
    char original_b;
    char* patternstr=NULL;
    char* rcpatternstr=NULL;
    hash* h;
    int n;
 

    // 1. read pattern file to get the number of keys
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


    fread(&n,sizeof(int),1,pattern_file);
    h=newhash(n);

    max_index=0;
    k=0;
    while(fread(&patternline,sizeof(Patternline),1,pattern_file) != 0){
       convert2str(patternline.pattern,&patternstr,&rcpatternstr);
       store_each_pattern(patternstr,0,patternline.index,h);  // reference allele
       store_each_pattern(rcpatternstr,0,patternline.index,h);  // reference allele, rc
       if(k==0) offset=(PATTERNLEN-1)/2;
       else if(k==1) offset=PATTERNLEN-1;
       else offset=0;
       original_b = patternstr[offset];
       for(j=0;j<4;j++){
          if(base[j]!=original_b) {
             patternstr[offset]=base[j];
             rcpatternstr[PATTERNLEN-offset-1]=rcbase[j];
             store_each_pattern(patternstr,1,patternline.index,h);  // alt allele
             store_each_pattern(rcpatternstr,1,patternline.index,h);  // alt allele, rc
          }
       }
       if(patternline.index>max_index) max_index=patternline.index;
       k++;
       if(k==3) k=0;
       if(patternstr!=NULL) free(patternstr); patternstr=NULL;
       if(rcpatternstr!=NULL) free(rcpatternstr); rcpatternstr=NULL;
    }
    fclose(pattern_file);
    return(h);
}


void store_each_pattern (char* pattern, char ref_or_alt, int index, hash* h){
       int2 hashval;
       hashval.c=ref_or_alt;
       hashval.i=index;
       insert_into_hash(pattern,strlen(pattern),hashval,h);
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



long*** build_count_array(void)
{
   int i;
   int ti;
   long*** count_arrays = (long***)malloc(nthread*sizeof(long**));
   for(ti=0;ti<nthread;ti++){
     count_arrays[ti]=(long**)malloc(2*sizeof(long*));
     for(i=0;i<2;i++) count_arrays[ti][i]=(long*)calloc((max_index+2),sizeof(long));
   }
   return(count_arrays);
}

void get_readlength_nReads (char* fastqfilename)
{
    FILE* fastq_file=0;
    gzFile* gz_fastq_file=0;
    char line[SINGLE_LINE_MAX],seq[SINGLE_LINE_MAX];
    char gzmode=0;

    if(strlen(fastqfilename)==0){ fastq_file = stdin; }
    else {
      gz_fastq_file = gzopen((const char*)fastqfilename,"r");
      if(gz_fastq_file !=0){ gzmode=1; }
      else {
        fastq_file = fopen((const char*)fastqfilename,"r");
        if(fastq_file ==0)
        {
           fprintf(stderr,"can't open fastq file.\n");
           exit(1);
        }
      }
    }

    if(gzmode){
       // first read
       if(gzgets(gz_fastq_file, line, SINGLE_LINE_MAX) != NULL && gzgets(gz_fastq_file, seq, SINGLE_LINE_MAX) != NULL){
            if(seq[strlen(seq)-1]=='\n') seq[strlen(seq)-1]='\0'; /*chomp*/
            read_length = strlen(seq);
            gzgets(gz_fastq_file, line, SINGLE_LINE_MAX); // ignore 3rd line
            gzgets(gz_fastq_file, line, SINGLE_LINE_MAX); // ignore 4th line
       }
       nReads=1;
       // all the other reads
       while (gzgets(gz_fastq_file, line, SINGLE_LINE_MAX) != NULL) {
            gzgets(gz_fastq_file, seq, SINGLE_LINE_MAX); // ignore 2nd line
            gzgets(gz_fastq_file, line, SINGLE_LINE_MAX); // ignore 3rd line
            gzgets(gz_fastq_file, line, SINGLE_LINE_MAX); // ignore 4th line
            nReads++;  // just count
       }
       gzclose(gz_fastq_file);

    }else{
       // first read
       if(fgets(line, SINGLE_LINE_MAX, fastq_file) != NULL && fgets(seq, SINGLE_LINE_MAX, fastq_file) != NULL){
              if(seq[strlen(seq)-1]=='\n') seq[strlen(seq)-1]='\0'; /*chomp*/
              read_length = strlen(seq);
            fgets(line, SINGLE_LINE_MAX, fastq_file); // ignore 3rd line
            fgets(line, SINGLE_LINE_MAX, fastq_file); // ignore 4th line
       }
       nReads=1;
       // all the other reads
       while (fgets(line, SINGLE_LINE_MAX, fastq_file) != NULL) {
            fgets(seq, SINGLE_LINE_MAX, fastq_file); // ignore 2nd line
            fgets(line, SINGLE_LINE_MAX, fastq_file); // ignore 3rd line
            fgets(line, SINGLE_LINE_MAX, fastq_file); // ignore 4th line
            nReads++;  // just count
       }
       fclose(fastq_file);
    }
}


double guess_subsampling_rate(double desired_depth, double reference_length, char pe){
    double estimated_depth,subsampling_rate;

    if(pe==0) estimated_depth = (double)(read_length * (double)nReads/1.0E6) / (reference_length/1.0E6);
    else estimated_depth = (double)(read_length * 2 * (double)nReads/1.0E6) / (reference_length/1.0E6);
    subsampling_rate = desired_depth / estimated_depth;

    fprintf(stderr,"read length = %d, nReads= %d, reference_length= %e, estimated depth without subsampling =%f, subsampling rate=%f, reads sampled (approximate) =%d\n",read_length,nReads,reference_length,estimated_depth,subsampling_rate,(int)((double)nReads*subsampling_rate)); // DEBUGGING

    if(subsampling_rate>1.0){
      if( nodeptherror==1) { 
        fprintf(stderr,"Warning: desired depth is higher than estimated depth, resulting in subsampling rate higher than 1.0. Resetting it to 1.0 - note that the depth will be lower than the desired depth.\n");
        subsampling_rate=1.0;
      }else{
        fprintf(stderr,"Error: desired depth is higher than estimated depth, resulting in subsampling rate higher than 1.0. Please enter a lower depth.\n");
        exit(1);
      }
    } 

    return(subsampling_rate);
}


void* read_fastq_thread (void* arg){
  read_fastq_args a = *((read_fastq_args*)arg);
  if(strlen(a.fastqfilename2)==0) read_fastq(a.fastqfilename,a.h,a.count_array,a.patternlength,a.subsampling_rate,a.startpos,a.endpos);
  else read_fastq_PE(a.fastqfilename,a.fastqfilename2,a.h,a.count_array,a.patternlength,a.subsampling_rate,a.startpos,a.endpos);
}


// read within [startpos,endpos)
void read_fastq (char* fastqfilename, hash* h, long** count_array, int patternlength, double subsampling_rate, long startpos, long endpos)
{
    FILE* fastq_file=0;
    gzFile* gz_fastq_file=0;
    char line[SINGLE_LINE_MAX],seq[SINGLE_LINE_MAX];
    int i,last_i;
    int2 val;
    long n;
    char gzmode=0;

    if(strlen(fastqfilename)==0){ fastq_file = stdin; }
    else {
      gz_fastq_file = gzopen((const char*)fastqfilename,"r");
      if(gz_fastq_file !=0){ gzmode=1; }
      else {
        fastq_file = fopen((const char*)fastqfilename,"r");
        if(fastq_file ==0)
        {
           fprintf(stderr,"can't open fastq file.\n");
           exit(1);
        }
      }
    }

    if(gzmode){
       // move to startpos
       n=0;
       while(n<startpos) {
         gzgets(gz_fastq_file, line, SINGLE_LINE_MAX);
         gzgets(gz_fastq_file, line, SINGLE_LINE_MAX);
         gzgets(gz_fastq_file, line, SINGLE_LINE_MAX);
         gzgets(gz_fastq_file, line, SINGLE_LINE_MAX);
         n++;
       }
   
       if(subsampling_rate==1.){
         while (gzgets(gz_fastq_file, line, SINGLE_LINE_MAX) != NULL && n<endpos) {
            if(gzgets(gz_fastq_file, seq, SINGLE_LINE_MAX)!=NULL){
              if(seq[strlen(seq)-1]=='\n') seq[strlen(seq)-1]='\0'; /*chomp*/
   
              last_i=strlen(seq)-patternlength;
              for(i=0;i<=last_i;i++){
                val=search_hash(seq+i,patternlength,h);
                if(val.i!=-1) { count_array[val.c][val.i]++; break; }
              }
              gzgets(gz_fastq_file, line, SINGLE_LINE_MAX); // ignore 3rd line
              gzgets(gz_fastq_file, line, SINGLE_LINE_MAX); // ignore 4th line
            }
            n++;
         }
       }else{
         while (gzgets(gz_fastq_file, line, SINGLE_LINE_MAX) != NULL && n<endpos) {
           if(gzgets(gz_fastq_file, seq, SINGLE_LINE_MAX)!=NULL){
            if( ((double)rand()/(double)RAND_MAX )<subsampling_rate){
              if(seq[strlen(seq)-1]=='\n') seq[strlen(seq)-1]='\0'; /*chomp*/
   
              last_i=strlen(seq)-patternlength;
              for(i=0;i<=last_i;i++){
                val=search_hash(seq+i,patternlength,h);
                if(val.i!=-1) { count_array[val.c][val.i]++; break; } // val.c==2 means decoy
              }
            }         
            gzgets(gz_fastq_file, line, SINGLE_LINE_MAX); // ignore 3rd line
            gzgets(gz_fastq_file, line, SINGLE_LINE_MAX); // ignore 4th line
           }
           n++;
         }
       }
       gzclose(gz_fastq_file);

    }else{

       // move to startpos
       n=0;
       while(n<startpos) {
         fgets(line, SINGLE_LINE_MAX, fastq_file);
         fgets(line, SINGLE_LINE_MAX, fastq_file);
         fgets(line, SINGLE_LINE_MAX, fastq_file);
         fgets(line, SINGLE_LINE_MAX, fastq_file);
         n++;
       }
   
       if(subsampling_rate==1.){
         while (fgets(line, SINGLE_LINE_MAX, fastq_file) != NULL && n<endpos) {
            if(fgets(seq, SINGLE_LINE_MAX, fastq_file)!=NULL){
              if(seq[strlen(seq)-1]=='\n') seq[strlen(seq)-1]='\0'; /*chomp*/
   
              last_i=strlen(seq)-patternlength;
              for(i=0;i<=last_i;i++){
                val=search_hash(seq+i,patternlength,h);
                if(val.i!=-1) { count_array[val.c][val.i]++; break; }
              }
              fgets(line, SINGLE_LINE_MAX, fastq_file); // ignore 3rd line
              fgets(line, SINGLE_LINE_MAX, fastq_file); // ignore 4th line
            }
            n++;
         }
       }else{
         while (fgets(line, SINGLE_LINE_MAX, fastq_file) != NULL && n<endpos) {
           if(fgets(seq, SINGLE_LINE_MAX, fastq_file)!=NULL){
            if( ((double)rand()/(double)RAND_MAX )<subsampling_rate){
              if(seq[strlen(seq)-1]=='\n') seq[strlen(seq)-1]='\0'; /*chomp*/
   
              last_i=strlen(seq)-patternlength;
              for(i=0;i<=last_i;i++){
                val=search_hash(seq+i,patternlength,h);
                if(val.i!=-1) { count_array[val.c][val.i]++; break; } // val.c==2 means decoy
              }
            }         
            fgets(line, SINGLE_LINE_MAX, fastq_file); // ignore 3rd line
            fgets(line, SINGLE_LINE_MAX, fastq_file); // ignore 4th line
           }
           n++;
         }
       }
       fclose(fastq_file);
    }
}




void read_fastq_PE (char* fastqfilename, char* fastqfilename2, hash* h, long** count_array, int patternlength, double subsampling_rate, long startpos, long endpos)
{
    FILE* fastq_file=0;
    FILE* fastq_file2=0;
    gzFile* gz_fastq_file=0;
    gzFile* gz_fastq_file2=0;
    char line[SINGLE_LINE_MAX],seq[SINGLE_LINE_MAX],line2[SINGLE_LINE_MAX],seq2[SINGLE_LINE_MAX];
    int i,last_i,last_i2;
    int2 val;
    long n,n2;
    char gzmode=0;

    if(strlen(fastqfilename)==0||strlen(fastqfilename2)==0){ fprintf(stderr,"Can't open fastq file.\n"); exit(1); }
    else {
      gz_fastq_file = gzopen((const char*)fastqfilename,"r");
      gz_fastq_file2 = gzopen((const char*)fastqfilename2,"r");
      if(gz_fastq_file !=0 && gz_fastq_file2 !=0){ gzmode=1; }
      else if(gz_fastq_file ==0 && gz_fastq_file2 ==0){
        fastq_file = fopen((const char*)fastqfilename,"r");
        fastq_file2 = fopen((const char*)fastqfilename2,"r");
        if(fastq_file ==0 || fastq_file2 == 0)
        {
           fprintf(stderr,"can't open fastq file.\n");
           exit(1);
        }
      }
      else { fprintf(stderr,"Both fastq files must be either gzipped or non-gzipped.\n"); exit(1); }
    }


    if(gzmode){

       // move to startpos
       n=0;
       while(n<startpos) {
         gzgets(gz_fastq_file,line, SINGLE_LINE_MAX);
         gzgets(gz_fastq_file,line, SINGLE_LINE_MAX);
         gzgets(gz_fastq_file,line, SINGLE_LINE_MAX);
         gzgets(gz_fastq_file,line, SINGLE_LINE_MAX);
         n++;
       }
       n2=0;
       while(n2<startpos) {
         gzgets(gz_fastq_file2, line, SINGLE_LINE_MAX);
         gzgets(gz_fastq_file2, line, SINGLE_LINE_MAX);
         gzgets(gz_fastq_file2, line, SINGLE_LINE_MAX);
         gzgets(gz_fastq_file2, line, SINGLE_LINE_MAX);
         n2++;
       }
  
   
       if(subsampling_rate==1.){
         while (gzgets(gz_fastq_file,line, SINGLE_LINE_MAX) != NULL && gzgets(gz_fastq_file2,line2, SINGLE_LINE_MAX) != NULL && n<endpos) {
           
           if(gzgets(gz_fastq_file, seq, SINGLE_LINE_MAX)!=NULL && gzgets(gz_fastq_file2, seq2, SINGLE_LINE_MAX)!=NULL){
              if(seq[strlen(seq)-1]=='\n') seq[strlen(seq)-1]='\0'; /*chomp*/
              if(seq2[strlen(seq2)-1]=='\n') seq2[strlen(seq2)-1]='\0'; /*chomp*/
   
              last_i=strlen(seq)-patternlength;
              for(i=0;i<=last_i;i++){
                val=search_hash(seq+i,patternlength,h);
                if(val.i!=-1) { count_array[val.c][val.i]++; break; } // val.c==2 means decoy
              }
              if(val.i==-1){ // If the pattern is not found in the first mate, search the second mate. If the pattern is found in the first mate, do not search the second mate, to avoid double-counting.
                last_i2=strlen(seq2)-patternlength;
                for(i=0;i<=last_i2;i++){
                  val=search_hash(seq2+i,patternlength,h);
                  if(val.i!=-1) { count_array[val.c][val.i]++; break; } // val.c==2 means decoy
                }
              }
   
              gzgets(gz_fastq_file,line, SINGLE_LINE_MAX); // ignore 3rd line
              gzgets(gz_fastq_file,line, SINGLE_LINE_MAX); // ignore 4th line
              gzgets(gz_fastq_file2,line2, SINGLE_LINE_MAX); // ignore 3rd line
              gzgets(gz_fastq_file2,line2, SINGLE_LINE_MAX); // ignore 4th line
            }
            n++;
         }
       }else{
         while (gzgets(gz_fastq_file,line, SINGLE_LINE_MAX) != NULL && gzgets(gz_fastq_file2,line2, SINGLE_LINE_MAX) != NULL && n<endpos) {
           
           if(gzgets(gz_fastq_file, seq, SINGLE_LINE_MAX)!=NULL && gzgets(gz_fastq_file2, seq2, SINGLE_LINE_MAX)!=NULL){
            if( ((double)rand()/(double)RAND_MAX )<subsampling_rate){
   
              if(seq[strlen(seq)-1]=='\n') seq[strlen(seq)-1]='\0'; /*chomp*/
              if(seq2[strlen(seq2)-1]=='\n') seq2[strlen(seq2)-1]='\0'; /*chomp*/
   
              last_i=strlen(seq)-patternlength;
              for(i=0;i<=last_i;i++){
                val=search_hash(seq+i,patternlength,h);
                if(val.i!=-1) { count_array[val.c][val.i]++; break; } // val.c==2 means decoy
              }
              if(val.i==-1){ // If the pattern is not found in the first mate, search the second mate. If the pattern is found in the first mate, do not search the second mate, to avoid double-counting.
                last_i2=strlen(seq2)-patternlength;
                for(i=0;i<=last_i2;i++){
                  val=search_hash(seq2+i,patternlength,h);
                  if(val.i!=-1) { count_array[val.c][val.i]++; break; } // val.c==2 means decoy
                }
              }
            }
            gzgets(gz_fastq_file,line, SINGLE_LINE_MAX); // ignore 3rd line
            gzgets(gz_fastq_file,line, SINGLE_LINE_MAX); // ignore 4th line
            gzgets(gz_fastq_file2,line2, SINGLE_LINE_MAX); // ignore 3rd line
            gzgets(gz_fastq_file2,line2, SINGLE_LINE_MAX); // ignore 4th line
           }
           n++;
         }
       }
       gzclose(gz_fastq_file);
       gzclose(gz_fastq_file2);

    }else{
       // move to startpos
       n=0;
       while(n<startpos) {
         fgets(line, SINGLE_LINE_MAX, fastq_file);
         fgets(line, SINGLE_LINE_MAX, fastq_file);
         fgets(line, SINGLE_LINE_MAX, fastq_file);
         fgets(line, SINGLE_LINE_MAX, fastq_file);
         n++;
       }
       n2=0;
       while(n2<startpos) {
         fgets(line, SINGLE_LINE_MAX, fastq_file2);
         fgets(line, SINGLE_LINE_MAX, fastq_file2);
         fgets(line, SINGLE_LINE_MAX, fastq_file2);
         fgets(line, SINGLE_LINE_MAX, fastq_file2);
         n2++;
       }
   
   
       if(subsampling_rate==1.){
         while (fgets(line, SINGLE_LINE_MAX, fastq_file) != NULL && fgets(line2, SINGLE_LINE_MAX, fastq_file2) != NULL && n<endpos) {
           
           if(fgets(seq, SINGLE_LINE_MAX, fastq_file)!=NULL && fgets(seq2, SINGLE_LINE_MAX, fastq_file2)!=NULL){
              if(seq[strlen(seq)-1]=='\n') seq[strlen(seq)-1]='\0'; /*chomp*/
              if(seq2[strlen(seq2)-1]=='\n') seq2[strlen(seq2)-1]='\0'; /*chomp*/
   
              last_i=strlen(seq)-patternlength;
              for(i=0;i<=last_i;i++){
                val=search_hash(seq+i,patternlength,h);
                if(val.i!=-1) { count_array[val.c][val.i]++; break; } // val.c==2 means decoy
              }
              if(val.i==-1){ // If the pattern is not found in the first mate, search the second mate. If the pattern is found in the first mate, do not search the second mate, to avoid double-counting.
                last_i2=strlen(seq2)-patternlength;
                for(i=0;i<=last_i2;i++){
                  val=search_hash(seq2+i,patternlength,h);
                  if(val.i!=-1) { count_array[val.c][val.i]++; break; } // val.c==2 means decoy
                }
              }
   
              fgets(line, SINGLE_LINE_MAX, fastq_file); // ignore 3rd line
              fgets(line, SINGLE_LINE_MAX, fastq_file); // ignore 4th line
              fgets(line2, SINGLE_LINE_MAX, fastq_file2); // ignore 3rd line
              fgets(line2, SINGLE_LINE_MAX, fastq_file2); // ignore 4th line
            }
            n++;
         }
       }else{
         while (fgets(line, SINGLE_LINE_MAX, fastq_file) != NULL && fgets(line2, SINGLE_LINE_MAX, fastq_file2) != NULL && n<endpos) {
           
           if(fgets(seq, SINGLE_LINE_MAX, fastq_file)!=NULL && fgets(seq2, SINGLE_LINE_MAX, fastq_file2)!=NULL){
            if( ((double)rand()/(double)RAND_MAX )<subsampling_rate){
   
              if(seq[strlen(seq)-1]=='\n') seq[strlen(seq)-1]='\0'; /*chomp*/
              if(seq2[strlen(seq2)-1]=='\n') seq2[strlen(seq2)-1]='\0'; /*chomp*/
   
              last_i=strlen(seq)-patternlength;
              for(i=0;i<=last_i;i++){
                val=search_hash(seq+i,patternlength,h);
                if(val.i!=-1) { count_array[val.c][val.i]++; break; } // val.c==2 means decoy
              }
              if(val.i==-1){ // If the pattern is not found in the first mate, search the second mate. If the pattern is found in the first mate, do not search the second mate, to avoid double-counting.
                last_i2=strlen(seq2)-patternlength;
                for(i=0;i<=last_i2;i++){
                  val=search_hash(seq2+i,patternlength,h);
                  if(val.i!=-1) { count_array[val.c][val.i]++; break; } // val.c==2 means decoy
                }
              }
            }
            fgets(line, SINGLE_LINE_MAX, fastq_file); // ignore 3rd line
            fgets(line, SINGLE_LINE_MAX, fastq_file); // ignore 4th line
            fgets(line2, SINGLE_LINE_MAX, fastq_file2); // ignore 3rd line
            fgets(line2, SINGLE_LINE_MAX, fastq_file2); // ignore 4th line
           }
           n++;
         }
       }
       fclose(fastq_file);
       fclose(fastq_file2);



    }
}


// the array marking which SNP index was present in the hash.
void create_index_array(char* patternfilename){
    FILE* pattern_file=0;
    int n;
    Patternline patternline; 

    index_array=calloc(max_index+1,sizeof(char));
    
    // 1. read pattern file to get the number of keys
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

    fread(&n,sizeof(int),1,pattern_file);
    while(fread(&patternline,sizeof(Patternline),1,pattern_file) != 0) index_array[patternline.index]=1;

}


void print_count_array(long*** count_arrays, char* patternfilename){
  int i,ti;
  long totalcount,refcount,altcount;
  double vaf;

  //mark in index array which SNP index was present in the hash
  create_index_array(patternfilename);

  printf("index\tref\talt\tvaf\n");
  for(i=0;i<=max_index;i++){
     if(index_array[i]==0) printf("%d\tNA\tNA\tNA\n",i);
     else {
        refcount=0;
        altcount=0;
        totalcount=0;
        for(ti=0;ti<nthread;ti++){
          refcount += count_arrays[ti][0][i];
          altcount += count_arrays[ti][1][i];
        }
        totalcount = refcount+altcount;
        if(totalcount>0) printf("%d\t%ld\t%ld\t%lf\n",i,refcount,altcount,(double)altcount / (double)totalcount);
        else printf("%d\t%ld\t%ld\tNA\n",i,refcount,altcount);
     }
  }

  //delete the index array
  free(index_array);
}


void printusage(char* program)
{
   printf("Usage : %s <options> -1 fastqfile1 [-2 fastqfile2]  patternfile(.pt)\n\n",program);

   printf("\tInput arguments (required)\n");
   printf("\t  patternfile : a binary file containing sequences flanking representative snv sites, along with markers indicating the snv index and whether the sequence represents reference or alternative allele.\n");
   printf("\t  -1, --fastq1 <fastq_file_1> : fastq file for SE data or first fastq file for a PE data. File can be gzipped (auto-detect).\n\n");

   printf("\tOptions\n");
   printf("\t  -2, --fastq2 <fastq_file_2> : second fastq file for a PE data. File can be gzipped (auto-detect)\n");
   printf("\t  -s, --ss <subsampling_rate> : subsampling rate (default 1.0)\n");
   printf("\t  -d, --depth <desired_depth> : as an alternative to a user-defined subsampling rate, let the program compute the subsampling rate given a user-defined desired_depth and the data.\n");
   printf("\t  -R, --reference_length <reference_length> : The reference length (default : 3E9) to be used for computing subsampling rate. If the data is NOT WGS from human, and if you're using the -d option, it is highly recommended to specify the reference length. For instance, if your data is human RNA-seq, the total reference length could be about 3%% of the human genome, which can be set as 1E8.\n");
   printf("\t  -L, --pattern_length <pattern_length> : The length of the flanking sequences being used to identify SNV sites. Default is 21bp. It is recommended not to change this value, unless you have created your own pattern file with a different pattern length.\n");
   printf("\t  -p, --maxthread <number_of_threads> : number of threads to use (default : 1 )\n\n");

}


