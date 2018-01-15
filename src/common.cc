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
 * Last Update    : June 30, 2013.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <dirent.h>
#include <errno.h>
#include <string>
#include <sstream>
#include <vector>
#include <map>

#include "common.h"

using namespace std;

//int READ_LENGTH = 75;
//int CPLEX_STR_LENGTH = 520;
//string DEST_FOLDER = "./temp/";
//
///**********************************************/
//int overlap_int( const int s1, const int e1, const int s2, const int e2)
//{
//	int length = 0;
//	if ( (s2 <= s1) && ( e1 <= e2) )
//	{	length = e1 - s1 + 1;	}
//	else if ( (s1 <= s2) && ( e2 <= e1) )
//	{	length = e2 - s2 + 1;	}
//	else if ( (s2 <= s1) && ( s1 <= e2))
//	{	length = e2 - s1 + 1;	}
//	else if ( (s2 <= e1) && ( e1 <= e2)) 
//	{	length = e1 - s2 + 1;	}
//	return length;
//}
//
//
///**********************************************/
//FILE *fileOpen(char *fileName, char *mode)
//{
//	FILE *fp;
//	fp = fopen(fileName, mode);
//	if (NULL == fp)
//	{
//		fprintf(stdout, "Error: Cannot Open the file %s\n", fileName);
//		fflush(stdout);
//		exit(0);
//	}
//	return fp;
//
//}
//
///**********************************************/
//vector<string> splitStringOld(string str, char delim)
//{
//	stringstream ss(str);
//	string item;
//	vector<string> token_array;
//
//	while(getline(ss, item, delim))
//	{
//		token_array.push_back(item);
//	}
//	return token_array;
//}
//
///**********************************************/
//vector<string> splitString(char *raw_str, char delim)
//{
//	vector<string> token_array;
//	token_array.reserve(20);
//	string tmp;
//	tmp.reserve(10000);
//
//	int l = strlen(raw_str);
//	int pos = 0;
//	while (pos < l)
//	{
//		if (raw_str[pos] == delim || raw_str[pos]=='\n')
//		{
//			token_array.push_back(tmp);
//			tmp = "";
//		}
//		else
//		{
//			tmp += raw_str[pos];
//		}
//		pos ++;
//
//	}
//	return token_array;
//}
//
//
///**********************************************/
//void splitStringCopy(const char *raw_str, const char delim, vector<string> &token_array)
//{
//	//vector<string> token_array;
//	//token_array.reserve(20);
//	string tmp;
//	tmp.reserve(10000);
//
//	int l = strlen(raw_str);
//	int pos = 0;
//	while (pos < l)
//	{
//		if (raw_str[pos] == delim || raw_str[pos]=='\n')
//		{
//			token_array.push_back(tmp);
//			tmp = "";
//		}
//		else
//		{
//			tmp += raw_str[pos];
//		}
//		pos ++;
//
//	}
//	//return token_array;
//}
//
/**********************************************/
void copyToString( char *src, string &new_str)
{
  int limit=strlen(src);
  new_str.clear();
  new_str.reserve(limit);
  for( int i = 0; i < limit; i ++)
  {
      new_str.push_back(src[i]);
  }
}
/**********************************************/
void attachToString( char *src, string &new_str)
{
  int limit=strlen(src);
  for( int i = 0; i < limit; i ++)
  {
      new_str += src[i];
  }
}

/**********************************************/
void copyToStringStrip( char *src, string &new_str)
{
  int limit=strlen(src);
  new_str.clear();
  new_str.reserve(limit);
  for( int i = 0; i < limit; i ++)
  {
	  if( '\n' != src[i])
	  {
		  new_str.push_back(src[i]);
	  }
  }
}

/**********************************************/
void stripString( string &s)
{
	if (!s.empty() && s[s.length()-1] == '\n') {
		s.erase(s.length()-1);
	}
}
/**********************************************/
int get_nth_word( char *target, const char *src, int n)
{
	int flag = 1;
	int of = 0;
	target[0]= '\0';
	while( n && (sscanf(src, "%s%n", target, &of)) )
	{
		src+=of;
		n--;
	}
	if ( n )
	{
		E("Too Many fields in %s\n", src);
		flag = 0;
	}
	return flag;
}

/**********************************************/
int get_num_word( const char *src )
{
	int flag = 0;
	int of = 0;

	const char * delim="\t";

	char *str1=(char*)malloc(MAX_LINE);
	char *token1=(char*)malloc(TOKEN_LENGTH);
	char *save;

	strncpy(str1, src, MAX_LINE);
	for (token1 = strtok_r(str1, delim, &save); token1; token1 = strtok_r(NULL, delim, &save))
	{	flag++; }
	
	free(token1);
	return flag;
}

/**********************************************/
void make_rc( const string &raw_str, string &new_str)
{
	new_str = "";
	int limit = (int)raw_str.size();
	for( int i =0; i < limit; i++ )
	{
		switch( raw_str[limit-1-i])
		{
			case 'A':
				new_str.push_back('T');
				break;
			case 'C':
				new_str.push_back('G');
				break;
			case 'G':
				new_str.push_back('C');
				break;
			case 'T':
				new_str.push_back('A');
				break;
			default:
				new_str.push_back('N');
		}
	}
}
/**********************************************/
void make_rev( const string &raw_str, string &new_str )
{
	new_str.clear();
	int limit = (int)raw_str.size();
	new_str.reserve(limit);
	for( int i =0; i < limit; i++ )
	{
		new_str.push_back( raw_str[limit-1-i] );
	}
}
///**********************************************/
//void ConfigureOutputFolder()
//{
//	struct stat st;
//	if ( 0 != stat(OutputFolder, &st) || !S_ISDIR(st.st_mode) )
//	{
//		if ( -1 == mkdir(OutputFolder, 0777) )
//		{
//			perror( "Error occurs in checking output folder" );
//			exit( EXIT_FAILURE );
//		}
//	}
//}
/**********************************************/
void ExtractSampleList( const char *list_file, vector<string> &sample_vector)
{
	int i = 0;

	FILE *fp = fopen( list_file, "r");
	char *readline = (char*)malloc(MAX_LINE);
	
	while( NULL != fgets( readline, MAX_LINE, fp) )
	{
		for( i = strlen(readline) - 1; i >= 0; i--)
		{
			if ( !isspace(readline[i]) ) {break;}
			E(">>%d %d %s", i+1, (int)strlen(readline), readline);
		}
		readline[i+1] = '\0';
		sample_vector.push_back(readline);
	}
	fclose( fp );
}
/**********************************************/
int over_l( const int s1, const int e1, const int s2, const int e2)
{
	int length = 0;
	if ( (s2 <= s1) && ( e1 <= e2) )
	{
		length = e1 - s1 + 1;
	}
	else if ( (s1 <= s2) && ( e2 <= e1) )
	{
		length = e2 - s2 + 1;
	}
	else if ( (s2 <= s1) && ( s1 <= e2))
	{
		length = e2 - s1 + 1;
	}
	else if ( (s2 <= e1) && ( e1 <= e2)) 
	{
		length = e1 - s2 + 1;
	}
	return length;
}

/**********************************************/
string add_extension( const char *src, const string &ext)
{
	string new_str;
	int limit=strlen(src);
	new_str.reserve(limit);
	
	for( int i = 0; i < limit; i ++)
	{
		new_str.push_back(src[i]);
	}
	new_str += ".";
	new_str += ext;
	return new_str;
}
/**********************************************/
string get_gtf_str( const string &src, uint32_t start, uint32_t end)
{
	uint32_t len = end - start + 1;
	return src.substr( start -1, len);
}

/**********************************************/
void fold_output( FILE *out, const string &src )
{
	int limit=(int)src.size();

	for( int i=0; i < limit; i++ )
	{
		fprintf( out, "%c", src[i]);
		if ( 0 == (i+1)%100)
		{fprintf(out, "\n");}
	}
	if ( 0 != limit%100 )
	{fprintf(out, "\n");}
}
