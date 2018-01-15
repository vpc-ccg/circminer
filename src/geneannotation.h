#ifndef __GENEANNOTATION__
#define __GENEANNOTATION__
#include <map>
#include <vector>
#include <string>
#include "common.h"
using namespace std;

void ensembl_gtf_reader( const char *gtf_file );

//typedef struct
//{
//	string chr_id;
//	int length;
//}tr_length;
//
//typedef struct
//{
//	string id;
//	string chr1;
//	string chr2;
//	char strand1;
//	char strand2;
//	vector<int> coor1;
//	vector<int> coor2;
//	int group;// for cluster indicator
//}fusion;
//
//
//typedef struct
//{
//	int index;
//	int end;
//	string chr;
//	int left;
//	int right;
//}f_segment;
//
//typedef struct
//{
//	int start;
//	int end;
//}Region;
//
//typedef struct
//{
//	uint32_t s;
//	uint32_t e;
//	int f;
//}CDS;
//
//typedef struct
//{
//	uint32_t e_s;
//	uint32_t e_e;
//	uint32_t c_s;
//	uint32_t c_e;
//	int f;
//	int cds; // 0 for UTR, 2 for CDS, 1 for state-changing exon
//}Exon;
//
//typedef struct
//{
//	string id;
//	string gid;
//	string ref;
//	string tname;
//	string gname;
//	string src;
//	string fea;
//	char strand;
//	vector<Exon> exon;
//	int cds_iso;
//}isoform;
//
//typedef struct
//{
//	uint32_t start;
//	uint32_t end;
//	char strand;
//	string chr;
//	string gene_id;
//	string gene_name;
//	int cds_gene;
//}gene_data;
//
//typedef struct
//{
//	string gene_id;
//	string trans_id;
//	string gene_name;
//	string trans_name;
//	string biotype;
//	int exon_num;
//} GTFAttr;
//
//typedef struct 
//{
//	string seqname;
//	string source;
//	string feature;
//	int start;
//	int end;
//	//float score;
//	char strand;
//	int frame;
//	GTFAttr attr;
//
//} GTFRecord;
//
//typedef struct
//{
//	int start;
//	int end;
//	int id;
//	char status;
//	int group_id;
//} ExonPosition;
//
//typedef struct
//{
//	int start;
//	int end;
//	char strand;
//	vector<ExonPosition> epos;
//	string chr_name;
//	string gene_name;
//	// We need to use gene NAME for some predictor. That's the only reason gn_name is here
//	string gn_name;
//} GENEPosition;
//
//typedef struct
//{
//	string gene_name;
//	int start;
//	int end;
//	char strand;
//	char *status;//char status;//int status;
//	FILE *fp;
//	int *count;
//	int finish_sample; // the number of sample finished
//	int flag;	// 1 indicates no reads in any sample; 0 o.w.
//	//vector<string> buffer;
//	vector< vector<string> > reads;
//} GENEBoundary;
//
//
//typedef struct
//{
//	string chr_name;
//	string gene_name;
//	int length;
//	int start;
//	int end;
//} coor_record;
//
//typedef struct
//{
//	string ref;
//	string fea;
//	uint32_t start;
//	uint32_t end;
//} t_mask;
//
//bool comp_isoform_start(const isoform &is_1, const isoform &is_2);
//bool comp_mask(const t_mask &m1, const t_mask &m2);
//bool comp_gene(const GENEBoundary &gb1, const GENEBoundary &gb2);
//bool comp_gtf( const GTFRecord &gtf1, const GTFRecord &gtf2);
//bool comp_exon( const Exon &exon1, const Exon &exon2); 
//bool comp_gene_in_contig( const gene_data &g1, const gene_data &g2); 
//bool comp_gene_in_contig_end( const gene_data &g1, const gene_data &g2); 
//
//void adjust_gene( map<string, gene_data> &map_gene, const string g_id, char *gene_name, const isoform &new_iso);
//
//void ensembl_Reader( char *gtf_file, map<string, vector<isoform> > &iso_gene_map, map<string, vector<gene_data> > &gene_sorted_map);
//void output_Isoform_Boundary( map<string, vector<isoform> > &iso_gene_map, char *outfile );
//void masker_reader( char *gtf_file, char *out_file );
//
//void output_Isoform_Boundary_Specific( map<string, vector<isoform> > &iso_gene_map, char *outfile, map<string, int> &gene_dict );
//void output_Gene_Boundary_Specific( map<string, vector<gene_data> > &gene_sorted_map, char *outfile, map<string, int> &gene_dict );
//void genename_reader( char *gtf_file, map<string, int> &gene_dict );
//
//void SortGeneInContig( map<string, vector<gene_data> > &gene_sorted_map);
//void SortGeneInContigByEnd( map<string, vector<gene_data> > &gene_sorted_map);
//
//
//void LengthCalculator( char *gtf_file, map<string, vector<tr_length> > &length_map );
////void LengthCalculator( char *gtf_file, map<string, int> &length_vec );
//
//void regionReader( char *pred_file, vector<fusion> &fusion_vec );
//
//void read_uniq_region( std::map<std::string, std::vector<Region> > &UniqMap, char *uniq_file);
//
//void LoadCoordinate(char *info_file, map<string, coor_record> &coor_map);
//void mapping_transcript(string id_str, string name_string, map<string, string> &enst_table);
//
//map<string, GENEPosition> gtfReader_PT(char *gtfFileName, map<string, vector<string> > &gene_table, map<string, string> &enst_table, map<string, string> &ensg_table );
//map<string, map<string, GENEPosition> > gtfReader(char *gtfFileName );
//GENEPosition generate_genepos(vector<GTFRecord> exon_vector);
//GTFRecord MakeGtfItem(vector<string> token_array);
//GTFAttr parse_gtf_attr(string attr_string);
//string convert_genename(string raw_str);
//string extract_attr(string raw_string);
//
//map<string, vector<GENEBoundary> > convertGeneBoundary(map<string, map<string, GENEPosition> > GeneMap, int s_num);
//
//void AnalyzeAllMapping(char *map_file, map<string, coor_record> &coor_map, int read_length, char *log_file);
//void check_mapping(string &trans_name, vector<int> &mapping_count, int max_count, FILE *log_fp, string &chr_name, string &gene_name);
//
//
//void AnalyzeUniqMapping(char *map_file, map<string, coor_record> &coor_map, int read_length, char *log_file);
//void check_uniq_mapping(string &trans_name, vector<int> &mapping_count, int max_count, FILE *log_fp, string &chr_name, string &gene_name);
//
//void printGTF(vector<GENEBoundary> gb_vector);
//void printGeneMap(map<string, map<string, GENEPosition> > GeneMap);
//void printPTGeneMap( map<string, GENEPosition>  PTMap);
//void printExonPosition( vector<ExonPosition> epos);
//void printExonPositionRange( vector<ExonPosition> epos, int start, int end);
//void OutputGeneStat(map<string, GENEPosition> &PTMap);
//void OutputFusionPrediction( vector<fusion> &fusion_vec);
//void OutputLengthSummary( map<string, vector<tr_length> > &length_map);
#endif
