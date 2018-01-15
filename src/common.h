/*
 * Copyright (c) 2012 - 2013, Simon Fraser University
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without modification, 
 * are permitted provided that the following conditions are met:
 *   
 * Redistributions of source code must retain the above copyright notice, this list
 * of conditions and the following disclaimer.
 * - Redistributions in binary form must reproduce the above copyright notice, this
 *   list of conditions and the following disclaimer in the documentation and/or other
 *   materials provided with the distribution.
 * - Neither the name of the Simon Fraser University nor the names of its contributors may be
 *   used to endorse or promote products derived from this software without specific
 *   prior written permission.
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

/*
 * Author         : Yen-Yi Lin
 * Email          : yenyil AT sfu DOT ca
 * Last Update    : March 19, 2013.
 */

#ifndef __COMMON__
#define __COMMON__
#include <vector>
#include <map>
#include <string>
#include <inttypes.h>
#include <time.h>

// Maximum number of characters per line
#define MAX_LINE 200000
// Maximum length of file names
#define FILE_NAME_LENGTH 500
// Maximum number of characters in parsing
#define TOKEN_LENGTH 20000
// Prompting messages after reading fixed line
#define VERBOSE_LINE 100000

#define L(c,...) fprintf(stdout,c,##__VA_ARGS__)
#define E(c,...) fprintf(stderr,c,##__VA_ARGS__)

//extern char GtfFilename[FILE_NAME_LENGTH];
//extern char RepeatFilename[FILE_NAME_LENGTH];
//extern char FastqFilename[FILE_NAME_LENGTH];
//extern char SamFilename[FILE_NAME_LENGTH];
//extern char SeqFilename[FILE_NAME_LENGTH];
//extern char SamListname[FILE_NAME_LENGTH];
//extern char Predictname[FILE_NAME_LENGTH];
//extern char TargetFilename[FILE_NAME_LENGTH];
//extern char Outputname[FILE_NAME_LENGTH];
//extern char OutputFolder[FILE_NAME_LENGTH];
//extern char Libraryname[FILE_NAME_LENGTH];
//
//
//extern int Mode_TSVProtein;
//extern int Mode_CheckProtID;
//extern int Mode_CheckPepSEQ;
//extern int Mode_CheckDtaID;
//extern int Mode_CheckMultiHit;
//extern int Mode_CountEvent;
//// Probably splits here
//extern int Mode_GetFusionEvent_0;
//extern int Mode_GetFusionEvent;
//extern int Mode_GetSVEvent;
//extern int Mode_deFuseInfo;
//extern int Mode_GetFusionPep;
//extern int Mode_AdjustPvalue;
//extern int Mode_CheckVCF;
//extern int Mode_ConvertVCF;
//extern int Mode_GetSVPep;
//extern int Mode_MergePep;
//
////
//extern int Op_Mode;
//extern int Join_Mode;
//extern int Reads_Model;
//extern int Decoy_Mode;
//extern int Mix_Mode;
////extern int Verbose_Mode;
//extern int Library_Type;
//extern int Orman_Flag;
//extern int CPLEX_STR_LENGTH;
//extern int Abe_Mode;
//extern int Reads_Buffer;
//extern float Error_Bound;
//extern std::string DEST_FOLDER;
//
//
//extern char versionNumber[10];
//extern char versionNumberF[10];
//
//int overlap_int(const int s1, const int e1, const int s2, const int e2);
//
//FILE *fileOpen(char *fileName, char *mode);
//std::vector<std::string> splitStringOld(std::string str, char delim);
//std::vector<std::string> splitString(char *raw_str, char delim);
//
//void splitStringCopy(const char *raw_str, const char delim, std::vector<std::string> &token_array);
//std::string convert_illrgal_string(std::string raw_str);
//
// Utility Modules
void copyToString( char *src, std::string &new_str);
void attachToString( char *src, std::string &new_str);
void copyToStringStrip( char *src, std::string &new_str);
void stripString( std::string &s );
int get_nth_word( char *target, const char *src, int n);
int get_num_word( const char *src);
void make_rc( const std::string &raw_str, std::string &new_str);
void make_rev( const std::string &raw_str, std::string &new_str);
void ConfigureOutputFolder();
void ExtractSampleList( const char *listfile, std::vector<std::string> &sample_vector);
int over_l( const int s1, const int e1, const int s2, const int e2);
std::string add_extension( const char *src, const std::string &ext);
std::string get_gtf_str( const std::string &src, uint32_t start, uint32_t end );
void fold_output( FILE *out, const std::string &str);
#endif
