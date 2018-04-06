#ifndef __COMMON_H__
#define __COMMON_H__

#define FILE_NAME_LENGTH 1000

extern bool pairedEnd;

extern int kmer;
extern int verboseMode;

extern char gtfFilename[FILE_NAME_LENGTH];
extern char referenceFilename[FILE_NAME_LENGTH];
extern char fastqFilename[FILE_NAME_LENGTH];
extern char outputFilename[FILE_NAME_LENGTH];
extern char outputDir[FILE_NAME_LENGTH];

extern char versionNumberMajor[10];
extern char versionNumberMinor[10];

#endif
