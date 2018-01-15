#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <string.h>
#include <iostream>
#include <fstream>
#include <map>
#include <algorithm>
#include "geneannotation.h"

using namespace std;

/**********************************************/
void parse_attribute( const char *attr )
{
	E("Do %s", attr);	
	int offset;
	char *token 	= (char*)malloc(TOKEN_LENGTH);
	char *val 		= (char*)malloc(TOKEN_LENGTH);
	char c;
	while(attr)	
	{
		sscanf( attr, "%s \"%[^\"]\"%c %n", token, val, &c, &offset);
		L(">%s %s\n", token, val);
		//L("<%d %d\n", offset, attr+offset);
		//if ( strncmp( token, "gene_id", 7) )
		attr += offset;
	}
}


// Reading Gene modelGTF file with Ensembl format to get 
/**********************************************/
// void ensembl_gtf_reader( char *gtf_file, map<string, vector <isoform> > &iso_gene_map, map<string, vector<gene_data> > &gene_sorted_map)
void ensembl_gtf_reader( const char *gtf_file )//, map<string, vector <isoform> > &iso_gene_map, map<string, vector<gene_data> > &gene_sorted_map)
{
	char *readline 	= (char*)malloc(MAX_LINE);
	char *seq 		= (char*)malloc(TOKEN_LENGTH);
	char *src	 	= (char*)malloc(TOKEN_LENGTH);
	char *fea 		= (char*)malloc(TOKEN_LENGTH);
	char *misc 		= (char*)malloc(TOKEN_LENGTH);
	char *gid_str 	= (char*)malloc(TOKEN_LENGTH);
	char *gn_str 	= (char*)malloc(TOKEN_LENGTH);
	char *tid_str 	= (char*)malloc(TOKEN_LENGTH);
	char *tn_str 	= (char*)malloc(TOKEN_LENGTH);
	char *bio_str 	= (char*)malloc(TOKEN_LENGTH);
	
	uint32_t start = 0, end = 0; 
	
	char strand;
	char fr_cha;
	int frame = 0, offset = 0, of2 = 0;
	int count = 0;	
	
	int num_attr = 0;
	
	FILE *fp = fopen( gtf_file, "r" );
	while( NULL != fgets( readline, MAX_LINE, fp) )
	{
		if ( 0 == strncmp( "#", readline, 1 ) ) {	continue; }
		
		//sscanf( readline, "%s %s %s %u %u %s %c %c%n",
		sscanf( readline, "%s %s %s %u %u %s %c %c %n",
			seq, src, fea,  &start, &end, misc, &strand, &fr_cha, &offset );
		L("=%s %s\n", seq, src);
		L("=%s", readline+offset);
		parse_attribute( readline+offset);

		////L(">%s %s %s %u %u %s %c %c\n", seq, src, fea, start, end, misc, strand, fr_cha);
		////L("=%d\n", offset);
		////sscanf(readline +offset, "%[^\"]\"%[^\"]\"%[^\"]\"%[^\"]\" %n", // %s %s %u %u %s %c %d %n",
		////	misc, gid_str, misc, tid_str, &offset2//misc, tid_str, gn_st
		//////	seq, src, fea, &start, &end, score, &strand, &frame, &offset
		////);
		//
		//if ( !feature_of_interest( fea ) ) { continue; }
		//sscanf( readline + offset, 
		//	"%s %c%[^\"]%s %s %c%[^\"]%s\ 
		//	%s %s %s %c%[^\"]%s\
		//	%s %s %s %c%[^\"]%s %s %c%[^\"]%s %n",
		//	misc, &del, gid_str, misc, misc, &del, tid_str, misc,
		//	misc, misc, misc, &del, gn_str, misc, 
		//	misc, misc, misc, &del, bio_str, misc, misc, &del, tn_str, misc, &of2	);

		//// extra string copy slows down the following function
		//// num_attr = get_gtf_info( readline + offset, gtf_attr );
		//// L(">>%s\n", gid_str);L("$$%s\n", tid_str);L("^^%s\n", gn_str); //L("<<%s\n", tn_str); //L("--%s\n", bio_str); 
		//
		//// Start New Transcript
		//if ( 0 != strncmp( trans_id, tid_str, TOKEN_LENGTH ) )
		//{
		//	if ( '\0' != trans_id[0] )
		//	{	
		//		if ( 0 < exon_s + exon_e )	// last region is not a cds
		//		{	
		//			exon_token.e_s = exon_s;
		//			exon_token.e_e = exon_e;
		//			exon_token.c_s = cds_s;
		//			exon_token.c_e = cds_e;
		//			//// determine if a token is cds, utr, or intron
		//			exon_token.cds = 0;
		//			if ( 0 < cds_s + cds_e )
		//			{
		//				exon_token.cds = 2;
		//				if ( ( exon_s != cds_s) || ( exon_e != cds_e) )
		//				{	exon_token.cds = 1;	}
		//			}

		//			iso_token.exon.push_back(exon_token);	
		//		}
		//		//L("Finish Isoform %s with %d exons\n", trans_id, (int)iso_token.exon.size() ); 
		//		gstr = string(gene_id);				
		//		sort( iso_token.exon.begin(), iso_token.exon.end(), comp_exon);
		//		if ( 1 == cds_gene ){ iso_token.cds_iso = 1; }
		//		iso_gene_map[gstr].push_back( iso_token );

		//		if ( 0 == iso_token.exon.size()) { E("Warning: No Exon in %s\n", iso_token.id.c_str() ) ;}
		//		else{ adjust_gene( map_gene, gstr, gene_name, iso_token ); }
		//	}
		//	//L("Starting New Isoform %s\n", tid_str);
		//	strncpy( gene_id, gid_str, TOKEN_LENGTH);
		//	strncpy( gene_name, gn_str, TOKEN_LENGTH);
		//	strncpy( trans_id, tid_str, TOKEN_LENGTH);
		//	strncpy( trans_name, tn_str, TOKEN_LENGTH);
		//	iso_token.id = string(trans_id);
		//	iso_token.tname = string(trans_id);
		//	iso_token.gid = string(gene_id);
		//	iso_token.ref = string(seq);
		//	iso_token.src = string(src);
		//	iso_token.fea = string(fea);
		//	iso_token.strand = strand;
		//	iso_token.cds_iso = 0;
		//	iso_token.exon.clear();
		//	exon_s = 0; exon_e = 0; cds_gene = 0;
		//}

		//// CDS and Exon
		//// Note in Ensembl GTF CDS comes after Exon records
		////
		//if ( 0 == strncmp(fea, "exon", 4) )
		//{	// Previous records
		//	exon_token.e_s = exon_s;
		//	exon_token.e_e = exon_e;
		//	exon_token.c_s = cds_s;
		//	exon_token.c_e = cds_e;
		//	//// determine if a token is cds, utr, or intron
		//	exon_token.cds = 0;
		//	if ( 0 < cds_s + cds_e )
		//	{
		//		exon_token.cds = 2;
		//		if ( ( exon_s != cds_s) || ( exon_e != cds_e) )
		//		{	exon_token.cds = 1;	}
		//	}
		//	if ( 0 < exon_s + exon_e )	
		//	{	iso_token.exon.push_back(exon_token);	}

		//	
		//	// Reading Current Records
		//	exon_s = start;
		//	exon_e = end;
		//	cds_s = 0;
		//	cds_e = 0;
		//}
		//else if ( 0 == strncmp(fea, "CDS", 3) )
		//{
		//	cds_s = start;
		//	cds_e = end;
		//	cds_gene = 1; // reading a protein-coding transcript
		//}
		//
		count++;
		if ( 0 == count%1000000){E(".");}

	}
	////Last Record
	//exon_token.e_s = exon_s;
	//exon_token.e_e = exon_e;
	//exon_token.c_s = cds_s;
	//exon_token.c_e = cds_e;
	//exon_token.cds = 0;
	//if ( 0 < cds_s + cds_e )
	//{
	//	exon_token.cds = 2;
	//	if ( ( exon_s != cds_s) || ( exon_e != cds_e))
	//	{	exon_token.cds = 1;	}
	//}
	//if ( 0 < exon_s + exon_e )	
	//{	iso_token.exon.push_back(exon_token);	}
	////L("Finish Isoform %s with %d exon\n", trans_id, (int)iso_token.exon.size() );
	////strncpy( gene_name, gn_str, TOKEN_LENGTH);
	
	//gstr = string(gene_id);
	//sort( iso_token.exon.begin(), iso_token.exon.end(), comp_exon);
	//if ( 1 == cds_gene ){ iso_token.cds_iso = 1;}
	//iso_gene_map[gstr].push_back(iso_token);
	//if ( 0 == iso_token.exon.size()) { E("Warning: No Exon in %s\n", iso_token.id.c_str() ) ;}
	//else{adjust_gene( map_gene, gstr, gene_name, iso_token );}
	
	E("\n");
	free(readline);
	free(seq);
	free(src);
	free(fea);
	free(misc);
	free(gid_str);
	free(gn_str); 
	free(tid_str);
	free(tn_str);
	free(bio_str);
	//free(gene_id);
	//free(gene_name);
	//free(trans_id);
	//free(trans_name);

	//map<string, vector<isoform> >::iterator it;
	//for( it = iso_gene_map.begin(); it != iso_gene_map.end(); it++)
	//{
	//	int limit=(int)it->second.size();
	//	for( int i = 0; i < limit; i++)
	//	{
	//		int exon_size = it->second[i].exon.size();
	//		for( int j = 0; j < exon_size; j++)
	//		{	L("%s\t%s\t%lu\t%lu\t%d\t%d\n", it->first.c_str(), it->second[i].id.c_str(), it->second[i].exon[j].e_s, it->second[i].exon[j].e_e, it->second[i].exon[j].cds, it->second[i].cds_iso);	}
	//	}
	//}
	
	//map<string, gene_data >::iterator it;
	//for( it = map_gene.begin(); it != map_gene.end(); it++)
	//{
	//	//L("Add_Gene\t%s\t%s\t%s\t%u\t%u\t%d\n", it->second.chr.c_str(), it->second.gene_id.c_str(), it->second.gene_name.c_str(), it->second.start, it->second.end, it->second.cds_gene );
	//	gene_sorted_map[it->second.chr].push_back(it->second);
	//}
	//E("Scanning total %d genes\n", (int)map_gene.size() );
}

