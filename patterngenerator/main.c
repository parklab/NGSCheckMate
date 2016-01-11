#include<stdio.h>
#include "functions.h"

int main(int argc, char* argv[])
{

   if(argc<=1) { printusage_converter(argv[0]); return(0); }

   char patternfilename[FILENAMEMAX], outpatternfilename[FILENAMEMAX];
   strcpy(patternfilename,argv[1]);
   strcpy(outpatternfilename,argv[2]);

   convert_pattern_file(patternfilename,outpatternfilename);
  
   return(0);
}

