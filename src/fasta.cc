#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <string>
#include <map>

#include "fasta.h"
#include "common.h"

#define FASTA_LINE 50

using namespace std;

/**********************************************/
void LoadFasta( const char *fasta_file, map<string, string> &chr_map )
{
	FILE *fp = fopen( fasta_file, "r" );

	char *readline  = (char*)malloc(MAX_LINE);
	char *chr       = (char*)malloc(TOKEN_LENGTH);
	int str_len = 0, offset = 0;
	string tmp_str;
	string chr_name = "";
	string chr_seq = "";
	chr_seq.reserve(100000000);

	while( NULL != fgets(readline, MAX_LINE, fp) )
	{
		if ( '>' == readline[0] )
		{
			if ( 0 < (int)chr_name.size() )
			{
				chr_map[chr_name] = chr_seq;
				fprintf(stdout, "%s of length %d\n", chr_name.c_str(), (int)chr_seq.size()  );
				chr_seq= "";
				chr_seq.reserve(100000000);
			}

			//tmp_str = readline;
			//str_len = (int) tmp_str.length();
			sscanf(readline, ">%s %n\n", chr, &offset);
			copyToString( chr, chr_name);
		}
		else 
	    {
			tmp_str = readline;
			str_len = (int) tmp_str.length();

			chr_seq += tmp_str.substr(0, str_len-1);
			
	    }
	}

	chr_map[chr_name] =chr_seq;
	fprintf(stdout, "%s of length %d\n", chr_name.c_str(), (int)chr_seq.size()  );
	L("%lu sequences loaded.\n", chr_map.size() );
	fclose(fp);
	free(readline);
	free(chr);
}
///**********************************************/
//void LoadFastaShort( char *fasta_file, map<string, string> &chr_map, size_t buffer_size )
//{
//	
//	FILE *fp = fopen(fasta_file, "r");
//
//	char *readline  = (char*)malloc(MAX_LINE);
//	char *chr  = (char*)malloc(TOKEN_LENGTH);
//	int str_len = 0, offset = 0;
//	string  tmp_str;
//	string chr_name = "";
//	string chr_seq = "";
//	chr_seq.reserve(32);
//
//	while( NULL != fgets(readline, MAX_LINE, fp) )
//	{
//		if ( '>' == readline[0] )
//		{
//			if ( 0 < (int)chr_name.size() )
//			{
//				//fprintf(stdout, "%s of length %d\n", chr_name.c_str(), (int)chr_seq.size()  );
//				chr_map[chr_name] = chr_seq;
//				chr_seq.clear();
//				chr_seq.reserve(32);
//			}
//
//			tmp_str = readline;
//			str_len = (int) tmp_str.length();
//			sscanf(readline, ">%s %n\n", chr, &offset);
//			copyToString( chr, chr_name);
//		}
//		else 
//	    {
//			tmp_str = readline;
//			str_len = (int) tmp_str.length();
//
//			chr_seq += tmp_str.substr(0, str_len-1);
//			
//	    }
//	}
//
//	//fprintf(stdout, "%s of length %d\n", chr_name.c_str(), (int)chr_seq.size()  );
//	L("%lu sequences loaded.\n", (int)chr_map.size());
//	chr_map[chr_name] =chr_seq;
//	fclose(fp);
//	free(readline);
//	free(chr);
//}
//
///**********************************************/
//// Read Sequence Content for a FASTA file of single chromosome 
//void LoadFasta_Single( char *fasta_file, string &chr_name, string &chr_seq )
//{
//	
//	FILE *fp = fopen(fasta_file, "r");
//
//	char *readline  = (char*)malloc(MAX_LINE);
//	int str_len = 0;
//	string  tmp_str;
//	//string chr_name = "";
//	//string chr_seq = "";
//	chr_seq.reserve(100000000);
//
//	while( NULL != fgets(readline, MAX_LINE, fp) )
//	{
//		if ( '>' == readline[0] )
//		{
//			tmp_str = readline;
//			str_len = (int) tmp_str.length();
//			chr_name = tmp_str.substr(1, str_len - 2);
//		}
//		//{
//		//	if ( 0 < (int)chr_name.size() )
//		//	{
//		//		fprintf(stdout, "%s of length %d\n", chr_name.c_str(), (int)chr_seq.size()  );
//		//		chr_map[chr_name] = chr_seq;
//		//		chr_seq= "";
//		//		chr_seq.reserve(100000000);
//		//	}
//
//		//	tmp_str = readline;
//		//	str_len = (int) tmp_str.length();
//		//	
//		//	chr_name = tmp_str.substr(1, str_len - 2);
//		//}
//		else 
//		{
//			tmp_str = readline;
//			str_len = (int) tmp_str.length();
//			chr_seq += tmp_str.substr(0, str_len-1);
//		}	
//	    
//	}
//
//	fprintf(stdout, "%s of length %d\n", chr_name.c_str(), (int)chr_seq.size()  );
//	//chr_map[chr_name] =chr_seq;
//	fclose(fp);
//}
//// Reading Protein Fasta File for classifying a TSV results
///**********************************************/
//void readEnsemblProteinID( char *fasta_file, map<string, int> &prot_map )
//{
//	
//	FILE *fp = fopen(fasta_file, "r");
//
//	char *tag1=(char*)malloc(TOKEN_LENGTH);
//	char *tag2=(char*)malloc(TOKEN_LENGTH);
//	char *misc=(char*)malloc(TOKEN_LENGTH);
//
//	char *readline  = (char*)malloc(MAX_LINE);
//	int prot_type = 4; // If not specifiec, treat it as a known peptide
//	int str_len = 0;
//	int ret;
//	int offset;
//	int nov_pep = 0;
//
//	string chr_name = "";
//
//	while( NULL != fgets(readline, MAX_LINE, fp) )
//	{
//		// Header. In Ensembl, that will be something as 
//		// >ENSP00000354687 pep:known others
//		if ( 0 == strncmp(">", readline, 1) )
//		{
//			ret = sscanf( readline, ">%s %[^:]:%s %n", tag1, misc, tag2, &offset);
//			copyToString(tag1, chr_name);
//			
//			if ( 3 == ret ){	prot_type = getProteinType( tag2 ); }
//			else {prot_type = 4;}
//			
//			prot_map[chr_name] = prot_type;
//			//L("Adding %s with ret %d\n", chr_name.c_str(), ret);
//			//if ( 4 > prot_type){nov_pep++;}
//			if ( IsUnknownProtein(prot_type) ){nov_pep++;}
//		}
//	}
//
//	fprintf(stdout, "%ld records parsed with %d novel proteins\n", prot_map.size(), nov_pep  );
//	fclose(fp);
//	
//	free(tag1);
//	free(tag2);
//	free(misc);
//}
//// Each item is (peptide seq, Property )
///**********************************************/
//void readSeqEnsemblProteinID_0( char *fasta_file, map<string, int> &prot_map )
//{
//	
//	FILE *fp = fopen(fasta_file, "r");
//
//	char *tag1=(char*)malloc(TOKEN_LENGTH);
//	char *tag2=(char*)malloc(TOKEN_LENGTH);
//	char *misc=(char*)malloc(TOKEN_LENGTH);
//
//	char *readline  = (char*)malloc(MAX_LINE);
//	int prot_type = 3; // We treat it as a known peptide if the information NA
//	int str_len = 0;
//	int ret;
//	int offset;
//	int nov_pep = 0;
//
//	string chr_name = "";
//	string chr_seq = "";
//
//	while( NULL != fgets(readline, MAX_LINE, fp) )
//	{
//		// Header. In Ensembl, that will be something as 
//		// >ENSP00000354687 pep:known others
//		if ( 0 == strncmp(">", readline, 1) )
//		{
//			ret = sscanf( readline, ">%s %[^:]:%s %n", tag1, misc, tag2, &offset);
//			copyToString(tag1, chr_name);
//			
//			if ( 3 == ret ){	prot_type = getProteinType( tag2 ); }
//			else {prot_type = 3;}
//		
//			if ( 0 < (int)chr_seq.size() )
//			{	
//				prot_map[chr_seq] = prot_type; 
//				chr_seq.clear();
//			}
//			//L("Adding %s with ret %d\n", chr_name.c_str(), ret);
//			if ( 3 > prot_type){nov_pep++;}
//		}
//		else 
//		{
//			str_len = (int)strlen(readline);
//			for( int i = 0; i < str_len - 1; i++ )
//			{
//				chr_seq += readline[i];
//			}
//			
//		}
//	}
//			
//	if ( 0 < (int)chr_seq.size() )
//	{	
//		prot_map[chr_seq] = prot_type; 
//		chr_seq.clear();
//		if ( 3 > prot_type){nov_pep++;}
//	}
//
//	fprintf(stdout, "%ld records parsed with %d novel proteins\n", prot_map.size(), nov_pep  );
//	fclose(fp);
//	
//	free(tag1);
//	free(tag2);
//	free(misc);
//	free(readline);
//}
//// Each item is (Peptide seq, Property )
///**********************************************/
//void readSeqEnsemblProteinID( char *fasta_file, map<string, int> &k_map, map<string, int> &n_map, map<string, int> &p_map )
//{
//	
//	FILE *fp = fopen(fasta_file, "r");
//
//	char *tag1=(char*)malloc(TOKEN_LENGTH);
//	char *tag2=(char*)malloc(TOKEN_LENGTH);
//	char *misc=(char*)malloc(TOKEN_LENGTH);
//
//	char *readline  = (char*)malloc(MAX_LINE);
//	int prot_type = 4; // We treat it as a known peptide if the information NA
//	int str_len = 0;
//	int ret;
//	int offset;
//	int count = 0, nov_pep = 0;
//
//	string chr_name = "";
//	string chr_seq = "";
//
//	while( NULL != fgets(readline, MAX_LINE, fp) )
//	{
//		// Header. In Ensembl, that will be something as 
//		// >ENSP00000354687 pep:known others
//		if ( 0 == strncmp(">", readline, 1) )
//		{
//			if ( 0 < (int)chr_seq.size() )
//			{	
//				//if ( 4 == prot_type ){ k_map[chr_seq] = 1;}
//				//else if ( 2 == prot_type ){ n_map[chr_seq] = 1;}
//				//else if ( 1 == prot_type ){ p_map[chr_seq] = 1;}
//				if ( IsKnownProtein(prot_type) ){ k_map[chr_seq] = 1;}
//				else if ( IsNovelProtein( prot_type ) ){ n_map[chr_seq] = 1;}
//				else if ( IsPutativeProtein( prot_type ) ){ p_map[chr_seq] = 1;}
//				chr_seq.clear();
//				if ( IsUnknownProtein( prot_type ) ){nov_pep++;}
//				count++;
//			}// store the previous record
//
//			ret = sscanf( readline, ">%s %[^:]:%s %n", tag1, misc, tag2, &offset);
//			copyToString(tag1, chr_name);
//			
//			if ( 3 == ret ){	prot_type = getProteinType( tag2 ); }
//			else {prot_type = 4;}
//		}
//		else 
//		{
//			str_len = (int)strlen(readline);
//			for( int i = 0; i < str_len - 1; i++ )
//			{
//				chr_seq += readline[i];
//			}
//			
//		}
//	}
//	if ( 0 < (int)chr_seq.size() )
//	{	
//		//if ( 3 == prot_type ){ k_map[chr_seq] = 1;}
//		//else if ( 2 == prot_type ){ n_map[chr_seq] = 1;}
//		//else if ( 1 == prot_type ){ p_map[chr_seq] = 1;}
//		//if ( 3 > prot_type){nov_pep++;}
//		//count++;
//		if ( IsKnownProtein(prot_type) ){ k_map[chr_seq] = 1;}
//		else if ( IsNovelProtein( prot_type ) ){ n_map[chr_seq] = 1;}
//		else if ( IsPutativeProtein( prot_type ) ){ p_map[chr_seq] = 1;}
//		chr_seq.clear();
//		if ( IsUnknownProtein( prot_type ) ){nov_pep++;}
//		count++;
//	}// store the previous record
//			
//
//	fprintf(stdout, "%d records parsed with %d novel proteins\n", count++, nov_pep  );
//	fclose(fp);
//	
//	free(tag1);
//	free(tag2);
//	free(misc);
//	free(readline);
//}
///**********************************************/
//void MakeTranscript( map<string, GENEPosition> &PTMap, map<string, string> &chr_map, char *fasta_file, char *log_file, int run_mode, int read_length, int threshold  )
//{
//	//int fasta_line = 50;
//
//	map<string, GENEPosition>::iterator mit;
//	GENEPosition g_pos;
//	
//	int finished = 0;
//
//	int exon_num = 0;
//	
//	string trans_seq = "", exon_seq = "";
//	int trans_length;
//	int global_start = 1, global_end = 1;
//	int shift = 50;
//
//	string trans_name;
//	if ( 6 == run_mode )
//	{
//		FILE *fasta_fp = fopen(fasta_file, "a");
//		fprintf(fasta_fp, ">all_transcripts\n");
//		fclose(fasta_fp);
//	}
//	FILE *fp = fopen(log_file, "w");
//
//
//	for( mit = PTMap.begin(); mit != PTMap.end(); mit++ )
//	{
//		g_pos = (*mit).second;
//		exon_num = (int)g_pos.epos.size();
//
//		trans_seq = "";
//		trans_length = 0;
//		trans_seq.reserve(100000);
//
//		for( int i = 0; i < exon_num; i++ )
//		{
//			exon_seq = chr_map[g_pos.chr_name].substr(g_pos.epos[i].start-1, g_pos.epos[i].end - g_pos.epos[i].start  + 1 );
//			trans_seq += exon_seq;
//			trans_length += exon_seq.size();
//		}
//
//		//int line_num = 0;
//		//int last_num = 0;
//		trans_name = (*mit).first;
//		if ( 3 == run_mode )
//		{
//			if ( threshold <= trans_length )
//			{
//			write_fasta(fasta_file, trans_name, trans_seq, trans_length);
//			fprintf(fp, "%s\t%s\t%s\t%d\n", (*mit).second.chr_name.c_str(), (*mit).second.gene_name.c_str(), (*mit).first.c_str(), trans_length);
//			}
//			//fprintf(fp, "%s\n", trans_seq.c_str() );
//			//line_num = trans_length/fasta_line;
//			//last_num = trans_length%fasta_line;
//			//fprintf(fp, ">%d\t%d\n", trans_length, line_num );
//			//for( int i = 0; i < line_num; i++)
//			//{
//			//	fprintf(fp, "%d %s\n",i,  trans_seq.substr(50*i, 50).c_str() );
//			//}
//			//if ( 0 < last_num)
//			//{
//			//	fprintf(fp, "%d %s\n",line_num+1,  trans_seq.substr(50*line_num, last_num).c_str() );
//			//}
//		}
//		else if ( 4 == run_mode )
//		{
//			if ( threshold <= trans_length )
//			{
//				write_read_fasta(fasta_file, trans_name, trans_seq, trans_length, read_length, (*mit).second );
//				fprintf(fp, "%s\t%s\t%s\t%d\n", (*mit).second.chr_name.c_str(), (*mit).second.gene_name.c_str(), (*mit).first.c_str(), trans_length);
//				//write_read_fasta_separate(trans_name, trans_seq, trans_length, read_length, (*mit).second.chr_name );
//			}
//		}
//		else if ( 5 == run_mode )
//		{
//			if ( threshold <= trans_length )
//			{
//				write_annotate_fasta(fasta_file, trans_name, trans_seq, trans_length, (*mit).second);
//				fprintf(fp, "%s\t%s\t%s\t%c\t%d\n", (*mit).second.chr_name.c_str(), (*mit).second.gene_name.c_str(), (*mit).first.c_str(), (*mit).second.strand, trans_length);
//				//write_read_fasta_separate(trans_name, trans_seq, trans_length, read_length, (*mit).second.chr_name );
//			}
//		}
//		else if ( 6 == run_mode )
//		{
//			
//			if ( threshold <= trans_length )
//			{
//				write_concat_fasta(fasta_file, trans_name, trans_seq, trans_length);
//				global_end = global_start + trans_length - 1;
//				fprintf(fp, "%s\t%s\t%s\t%d\t%d\t%d\n", (*mit).second.chr_name.c_str(), (*mit).second.gene_name.c_str(), (*mit).first.c_str(), trans_length, global_start, global_end);
//				global_start = global_end + shift + 1;
//			}
//			//fprintf(fp, "%s\n", trans_seq.c_str() );
//			//line_num = trans_length/fasta_line;
//			//last_num = trans_length%fasta_line;
//			//fprintf(fp, ">%d\t%d\n", trans_length, line_num );
//			//for( int i = 0; i < line_num; i++)
//			//{
//			//	fprintf(fp, "%d %s\n",i,  trans_seq.substr(50*i, 50).c_str() );
//			//}
//			//if ( 0 < last_num)
//			//{
//			//	fprintf(fp, "%d %s\n",line_num+1,  trans_seq.substr(50*line_num, last_num).c_str() );
//			//}
//		}
//		//else if ( 5 == run_mode)
//		//{
//		//	write_stat( trans_name, exon_num, trans_length);
//		//}
//		finished++;
//		if ( 0 == finished%1000)
//		{
//			fprintf(stderr, ".");
//		}
//	}
//
//	if ( 6 == run_mode )
//	{
//		FILE *fasta_fp = fopen(fasta_file, "a");
//		fprintf(fasta_fp, "\n");
//		fclose(fasta_fp);
//	}
//	fprintf(stderr, "\n");
//	fclose( fp );
//}
///**********************************************/
//void write_fasta(char *fasta_file, string &pt_id, string &trans_seq, int trans_length)
//{
//	FILE *fp = fopen(fasta_file, "a");
//	fprintf(fp, ">%s\n%s\n", pt_id.c_str(), trans_seq.c_str());
//	//int line_num = trans_length/FASTA_LINE;
//	//int last_num = trans_length%FASTA_LINE;
//	//for( int i = 0; i < line_num; i++ )
//	//{
//	//	fprintf(fp, "%s\n", trans_seq.substr(50*i, 50).c_str() );
//	//}
//	//if ( 0 < last_num )
//	//{
//	//	fprintf(fp, "%s\n",  trans_seq.substr(50*line_num, last_num).c_str() );
//	//}
//	fclose(fp);
//}
//
///**********************************************/
//void write_annotate_fasta(char *fasta_file, string &pt_id, string &trans_seq, int trans_length, GENEPosition &gpos)
//{
//	FILE *fp = fopen(fasta_file, "a");
//	fprintf(fp, ">%s_%s_%s_%c\n", pt_id.c_str(), gpos.chr_name.c_str(), gpos.gene_name.c_str(), gpos.strand);
//	fprintf(fp, "%s\n", trans_seq.c_str() );
//	fclose(fp);
//}
///**********************************************/
//void write_concat_fasta(char *fasta_file, string &pt_id, string &trans_seq, int trans_length)
//{
//
//	FILE *fp = fopen(fasta_file, "a");
//	fprintf(fp, "%s", trans_seq.c_str());
//	fprintf(fp, "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN" );
//	fclose(fp);
//}
//
///**********************************************/
//void write_fasta_separate( string &pt_id, string &trans_seq, int trans_length, string &chr_name )
//{
//	string fasta_file = "ref_" + chr_name + ".fasta";
//	FILE *fp = fopen(fasta_file.c_str(), "a");
//
//	fprintf(fp, ">%s\n", pt_id.c_str());
//	int line_num = trans_length/FASTA_LINE;
//	int last_num = trans_length%FASTA_LINE;
//	for( int i = 0; i < line_num; i++ )
//	{
//		fprintf(fp, "%s\n",  trans_seq.substr(50*i, 50).c_str() );
//	}
//	if ( 0 < last_num )
//	{
//		fprintf(fp, "%s\n",  trans_seq.substr(50*line_num, last_num).c_str() );
//	}
//
//	fclose(fp);
//}
//
///**********************************************/
//void write_read_fasta(char *fasta_file, string &pt_id, string &trans_seq, int trans_length, int read_length, GENEPosition &gpos)
//{
//	FILE *fp = fopen(fasta_file, "a");
//
//	int max_start = trans_length - read_length;
//	int counter = 1; 
//	string read_seq = "";
//
//	//int line_num = read_length/FASTA_LINE;
//	//int last_num = read_length%FASTA_LINE;
//	
//	for(int x= 0; x < max_start; x++)
//	{
//		read_seq = trans_seq.substr(x, read_length);
//		fprintf(fp, ">%s_%s_%d\n", gpos.gene_name.c_str(), pt_id.c_str(), counter);
//		fprintf(fp, "%s\n", read_seq.c_str());
//		//for( int i = 0; i < line_num; i++ )
//		//{
//		//	fprintf(fp, "%s\n", read_seq.substr(50*i, 50).c_str() );
//		//}
//		//if ( 0 < last_num )
//		//{
//		//	fprintf(fp, "%s\n",  read_seq.substr(50*line_num, last_num).c_str() );
//		//}
//		counter++;
//	}
//
//	fclose(fp);
//}
//
///**********************************************/
//void write_read_fasta_separate(string &pt_id, string &trans_seq, int trans_length, int read_length, GENEPosition &gpos)
//{
//	string fasta_file = "read_" + gpos.chr_name + ".fasta";
//	FILE *fp = fopen(fasta_file.c_str(), "a");
//
//	int max_start = trans_length - read_length;
//	int counter = 1; 
//	string read_seq = "";
//
//	int line_num = read_length/FASTA_LINE;
//	int last_num = read_length%FASTA_LINE;
//	
//	for(int x= 0; x < max_start; x++)
//	{
//		read_seq = trans_seq.substr(x, read_length);
//		fprintf(fp, ">%s_%d\n", pt_id.c_str(), counter);
//		//fprintf(fp, "%s\n", read_seq.c_str());
//		for( int i = 0; i < line_num; i++ )
//		{
//			fprintf(fp, "%s\n", read_seq.substr(50*i, 50).c_str() );
//		}
//		if ( 0 < last_num )
//		{
//			fprintf(fp, "%s\n",  read_seq.substr(50*line_num, last_num).c_str() );
//		}
//		counter++;
//	}
//
//	fclose(fp);
//}
