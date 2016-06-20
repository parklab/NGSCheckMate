#include "ngscheckmate.h"
#include <getopt.h>

int main(int argc, char* argv[])
{


   char *patternfilename, fastqfilename[FILENAMEMAX], fastqfilename2[FILENAMEMAX];
   int patternlength;
   double subsampling_rate;
   double reference_length;
   double desired_depth;
   int option_index;
   char c;
   char pe;
   int i,ti;
   hash* h;
   long n_per_thread;
   
   if(argc<=1) { printusage(argv[0]); return(0); }

   /* options */
   static struct option long_options[]= {
      {"fastq1", required_argument, 0, '1' }, /* SE fastq file or first fastq file for PE data */
      {"fastq2", required_argument, 0, '2' }, /* second fastq file for PE data */
      {"ss", required_argument, 0, 's' }, /* subsampling rate */
      {"depth", required_argument, 0, 'd'}, /* desired depth */
      {"reference_length", required_argument, 0, 'R'}, /* length of reference, (eg. whole genome for WGS or whole genic region for RNA-seq) */
      {"pattern_length", required_argument, 0, 'L'}, /* length of pattern */
      {"maxthread", required_argument, 0, 'p' }, /* number of threads to use */
      //{"nodeptherror", no_argument, 0, 'j' },  /* produces warning instead of stopping when estimated subsampling rate is larger than 1. */
      {0,0,0,0}
   };

   /* default values */
   subsampling_rate=1;
   patternlength=21;
   reference_length=3E9;
   desired_depth=0; // default. If depth is 0, use subsampling rate instead of desired_depth.
   nthread=1;
   strcpy(fastqfilename,"");
   strcpy(fastqfilename2,"");
   nodeptherror=1; // always 1. (produce warning instead of stopping when estimated subsmapling rate is larger than 1.

   /* parsing options */
   while(1){
       option_index=0;
       c= getopt_long (argc, argv, "1:2:s:d:R:L:p", long_options, &option_index);

       if(c==-1) break;
       switch(c){
          case 0:  /* ? */
             if (long_options[option_index].flag!=0) break;
             fprintf(stdout,"option %s", long_options[option_index].name);
             if(optarg) fprintf(stdout," with arg %s",optarg);
             fprintf(stdout,"\n");
             break;
          case '1':
             strcpy(fastqfilename,optarg);   
             break;
          case '2':
             strcpy(fastqfilename2,optarg);   
             break;
          case 's':
             subsampling_rate = atof(optarg);   
             break;
          case 'd':
             desired_depth = atof(optarg);   
             break;
          case 'R':
             reference_length = atof(optarg);   
             break;
          case 'L':
             patternlength = atoi(optarg);   
             break;
          case 'p':
             nthread= atoi(optarg);
             if(nthread<1) { fprintf(stderr, "error: number of threads must be at least 1 (option -p).\n"); exit(1); }
             if(nthread>1){
                pth = (pthread_t*)malloc(sizeof(pthread_t)*(nthread-1));
                working_thread = (int*)malloc(sizeof(int)*(nthread-1));
             }
             break;
          case '?':
             fprintf(stderr,"error: unknown option?\n"); /*?*/
          default:
             return(0);
       }
   }

   if(strlen(fastqfilename)==0) { printf("At least one fastq file must be provided.(option -1/--fastq1)\n"); exit(0); }

   /* non-option arguments */
   if(optind<argc) {
      patternfilename = (char*)argv[optind];
   }
   else{ printusage(argv[0]); return(0); }

   if(strlen(fastqfilename2)>0) pe=1; else pe=0;  // PE vs SE.
   fprintf(stderr,"pe=%d\n",pe); // DEBUGGING

   srand((unsigned)time(NULL));
   fprintf(stderr,"started reading patternfile\n");
   h = read_patternfile_construct_hash (patternfilename);
   fprintf(stderr,"finished reading patternfile\n");
   long*** count_arrays = build_count_array();
   fprintf(stderr,"finished building count array\n");


   // estimate subsampling rate based on desired_depth and reference_length. This requires reading the fastq file once and assumes the read length is constant throughout the file.
   get_readlength_nReads (fastqfilename);
   fprintf(stderr,"Number of reads in the file : %ld\n",nReads); // DEBUGGING
   if(desired_depth>0){
     fprintf(stderr,"Calculating subsampling rate..\n");
     subsampling_rate = guess_subsampling_rate(desired_depth, reference_length, pe);
   }

   // scanning fastq file with a subsampling rate.
   // this part is threaded.
   fprintf(stderr,"Scanning fastq file..\n");
   read_fastq_args args[nthread];
   n_per_thread = nReads/nthread+1; // each scanning checks EOF (NULL) as well so to make sure all reads are read, I added 1.
   for(ti=0;ti<nthread;ti++){
     args[ti].fastqfilename = fastqfilename;
     args[ti].fastqfilename2 = fastqfilename2;
     args[ti].h = h;
     args[ti].count_array = count_arrays[ti];
     args[ti].patternlength = patternlength;
     args[ti].subsampling_rate = subsampling_rate;
     args[ti].startpos= ti==0?0:args[ti-1].endpos;
     args[ti].endpos=args[ti].startpos+n_per_thread;
   }

   if(nthread>1)
     for(ti=0;ti<nthread-1;ti++){
       working_thread[ti]=1;
       pthread_create(&pth[ti],NULL,read_fastq_thread,&(args[ti]));
       sleep(1);  // sleep 1 sec, so that not every rand function generates the same random number.
     }
   read_fastq_thread(&(args[nthread-1])); // unthreaded part.

   // Wait until all the threads finish.
   if(nthread>1)
     for(ti=0;ti<nthread-1;ti++) {
       pthread_join(pth[ti],NULL);
       working_thread[ti]=0;
      }

   fprintf(stderr,"finished reading fastq file\n");

   delete_hash(h); h=NULL;
   fprintf(stderr,"finished deleting hash\n");

   print_count_array(count_arrays);
   fprintf(stderr,"finished printing count array\n");

   //freeing count_array
   for(ti=0;ti<nthread;ti++){
     for(i=0;i<2;i++) free(count_arrays[ti][i]); count_arrays[ti][i]=NULL;
     free(count_arrays[ti]); count_arrays[ti]=NULL;
   }
   free(count_arrays); count_arrays=NULL;


   //allocated memory for threading.
   if(nthread>1){
       free(pth); pth=NULL;
       free(working_thread); working_thread=NULL;
   }

   return(0);
}

