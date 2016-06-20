#ifndef NGSCHECKMATE_PATTERNCONVERTER_H
#define NGSCHECKMATE_PATTERNCONVERTER_H

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#define FILENAMEMAX 1000
#define SINGLE_LINE_MAX 1000
#define PATTERNLEN 21
#define BPATTERNLEN 6   // length of binary string that represent pattern (=ceiling(PATTERNLEN/4))


typedef struct Patternline_str
{
   unsigned char pattern[BPATTERNLEN]; // This is not a string. It's an array of character (no '\0' at the end)
   int index;
} Patternline;


static char base[4]={'A','C','G','T'};
static char rcbase[4]={'T','G','C','A'};

unsigned char* convert2bi(char* str);
Patternline parse_patternfileline (char* line);
void convert_pattern_file(char* patternfilename, char* outpatternfilename);
void printusage(char*);

void convert2str(unsigned char*,char**,char**);


#endif

