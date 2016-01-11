#include "ngscheckmate.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include "stringhash2.h"

pthread_mutex_t mutex1 = PTHREAD_MUTEX_INITIALIZER;






int parse_patternfileline (char* line, hash* h){
       char tmpstr[SINGLE_LINE_MAX],pattern[SINGLE_LINE_MAX];
       int2 hashval;
       int numtab,j,i,index;
       char ref_vs_alt;

       numtab=0;j=0;
       for(i=0;i<=strlen(line);i++){
          if(line[i]=='\t'||line[i]=='\0') {
            numtab++;
            tmpstr[j]='\0';
            if(numtab==1) strcpy(pattern,tmpstr);
            else if(numtab==2) index = atoi(tmpstr);
            else if(numtab==3) ref_vs_alt = (char)atoi(tmpstr);
            j=0;
          }
          else { tmpstr[j]=line[i]; j++; }
       }
       if(numtab<3) { fprintf(stderr,"Error: input pattern file must have three columns.\n"); exit(1); }

       hashval.c=ref_vs_alt;
       hashval.i=index;

       insert_into_hash(pattern,strlen(pattern),hashval,h);
       return(index);
}


void parse_decoyfileline (char* line, hash* h){
       int2 hashval;

       // c and i values for decoy.
       hashval.c=0;
       hashval.i=max_index+1; // This is the indicator that the sequence is decoy. It is still counted but not included in the final printed array. 

       insert_into_hash(line,strlen(line),hashval,h);
}




hash* read_patternfile_construct_hash (char* patternfilename, char* decoyfilename){
    FILE* pattern_file=0;
    FILE* decoy_file=0;
    char line[SINGLE_LINE_MAX],line2[SINGLE_LINE_MAX];
    int index;
    long n;
    hash* h;


    // 1. read pattern and decoy files to get the number of keys
    n=0; // number of patterns
    // patternfile
    if(strlen(patternfilename)==0){ pattern_file = stdin; }
    else {
      pattern_file = fopen((const char*)patternfilename,"r");
      /*Check for validity of the file.*/
      if(pattern_file == 0)
      {
         fprintf(stderr,"can't open pattern file.\n");
         exit(1);
      }
    }

    while (fgets(line, SINGLE_LINE_MAX, pattern_file) != NULL) n++;
    fclose(pattern_file);

    // decoy file
    if(strlen(decoyfilename)>0){ 
      decoy_file = fopen((const char*)decoyfilename,"r");
      /*Check for validity of the file.*/
      if(decoy_file == 0)
      {
         fprintf(stderr,"can't open pattern file.\n");
         exit(1);
      }

      while (fgets(line, SINGLE_LINE_MAX, decoy_file) != NULL) n++;
      fclose(decoy_file);
    }

    h=newhash(n);


    // now actually build hash
    // pattern file
    if(strlen(patternfilename)==0){ pattern_file = stdin; }
    else {
      pattern_file = fopen((const char*)patternfilename,"r");
      /*Check for validity of the file.*/
      if(pattern_file == 0)
      {
         fprintf(stderr,"can't open pattern file.\n");
         exit(1);
      }
    }

    max_index=0;
    while (fgets(line, SINGLE_LINE_MAX, pattern_file) != NULL) {
       if(line[strlen(line)-1]=='\n') line[strlen(line)-1]='\0'; /*chomp*/
       index=parse_patternfileline(line, h);
       if(index>max_index) max_index=index;
    }
    fclose(pattern_file);

    // decoy file
    if(strlen(decoyfilename)>0){ 
      decoy_file = fopen((const char*)decoyfilename,"r");
      /*Check for validity of the file.*/
      if(decoy_file == 0)
      {
         fprintf(stderr,"can't open decoy file.\n");
         exit(1);
      }

      while (fgets(line, SINGLE_LINE_MAX, decoy_file) != NULL) {
       if(line[strlen(line)-1]=='\n') line[strlen(line)-1]='\0'; /*chomp*/
       parse_decoyfileline(line, h);
      }
      // do not add to index (decoys don't have snv index)
      fclose(decoy_file);
    }

    return(h);
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
    char line[SINGLE_LINE_MAX],seq[SINGLE_LINE_MAX];

    if(strlen(fastqfilename)==0){ fastq_file = stdin; }
    else {
      fastq_file = fopen((const char*)fastqfilename,"r");
      /*Check for validity of the file.*/
      if(fastq_file == 0)
      {
         fprintf(stderr,"can't open fastq file.\n");
         exit(1);
      }
    }

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
    char line[SINGLE_LINE_MAX],seq[SINGLE_LINE_MAX];
    int i,last_i;
    int2 val;
    long n;

    if(strlen(fastqfilename)==0){ fastq_file = stdin; }
    else {
      fastq_file = fopen((const char*)fastqfilename,"r");
      /*Check for validity of the file.*/
      if(fastq_file == 0)
      {
         fprintf(stderr,"can't open fastq file.\n");
         exit(1);
      }
    }

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




void read_fastq_PE (char* fastqfilename, char* fastqfilename2, hash* h, long** count_array, int patternlength, double subsampling_rate, long startpos, long endpos)
{
    FILE* fastq_file=0;
    FILE* fastq_file2=0;
    char line[SINGLE_LINE_MAX],seq[SINGLE_LINE_MAX],line2[SINGLE_LINE_MAX],seq2[SINGLE_LINE_MAX];
    int i,last_i,last_i2;
    int2 val;
    long n,n2;

    //fprintf(stderr,"scanning fastq in a paired-end mode..\n"); // DEBUGGING

    if(strlen(fastqfilename)==0||strlen(fastqfilename2)==0){ 
       fprintf(stderr,"Can't open fastq file.\n");
    }
    else {
      fastq_file = fopen((const char*)fastqfilename,"r");
      fastq_file2 = fopen((const char*)fastqfilename2,"r");
      /*Check for validity of the file.*/
      if(fastq_file == 0 || fastq_file2 == 0)
      {
         fprintf(stderr,"can't open fastq file.\n");
         exit(1);
      }
    }

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



void print_count_array(long*** count_arrays){
  int i,ti;
  long totalcount,refcount,altcount;
  double vaf;
  printf("index\tref\talt\tvaf\n");
  for(i=0;i<=max_index;i++){
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


void printusage(char* program)
{
   printf("Usage : %s <options> -1 fastqfile1 [-2 fastqfile2]  patternfile\n\n",program);

   printf("\tInput arguments (required)\n");
   printf("\t  patternfile : a text file with sequences flanking representative snv sites, along with markers indicating the snv index and whether the sequence represents reference or alternative allele.\n");
   printf("\t  fastqfile1 : see below 'Options'.\n\n");

   printf("\tOptions\n");
   printf("\t  -1, --fastq1 <fastq_file_1> : fastq file for SE data or first fastq file for a PE data. (required)\n");
   printf("\t  -2, --fastq2 <fastq_file_2> : second fastq file for a PE data.\n");
   printf("\t  -b, --bait <bait_file> : An optional a text file with decoy bait sequences (one per line). The sequences must be of the same length as the pattern length (default 21bp, or specified by -L) and located far enough from the pattern sequences in the genome that the presence of a bait in a read indicates the absence of a pattern in the read. Using a decoy bait file can improve speed. Not recommended for spliced reads (e.g. RNA-seq) unless the bait file is splice-aware (currently unsupported).\n");
   printf("\t  -s, --ss <subsampling_rate> : subsampling rate (default 1.0)\n");
   printf("\t  -d, --depth <desired_depth> : as an alternative to a user-defined subsampling rate, let the program compute the subsampling rate given a user-defined desired_depth and the data.\n");
   printf("\t  -R, --reference_length <reference_length> : The reference length (default : 3E9) to be used for computing subsampling rate. If the data is NOT WGS from human, and if you're using the -d option, it is highly recommended to specify the reference length. For instance, if your data is human RNA-seq, the total reference length could be about 3%% of the human genome, which can be set as 1E8.\n");
   printf("\t  -L, --pattern_length <pattern_length> : The length of the flanking sequences being used to identify SNV sites. Default is 21bp. It is recommended not to change this value, unless you have created your own pattern file with a different pattern length.\n\n");

   printf("\t  -p, --maxthread <number_of_threads> : number of threads to use (default : 1 )\n");
   printf("\t  -j, --nodeptherror : in case estimated subsampling rate is larger than 1, do not stop but reset it to 1 and continue.\n");
}


