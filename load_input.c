
#include <unistd.h>
#include <malloc.h>
#include <string.h>
#include <dirent.h>
#include <inttypes.h>
#include <inttypes.h>

#include "load_input.h"

//debug
#include "bit_operation.h"

//KSEQ_INIT(gzFile, gzread)
kseq_t *seq1 = NULL;
kseq_t *seq2 = NULL;

uint64_t* buffer_ref_seq = NULL;
#ifdef UNI_SEQ64
uint64_t* buffer_seq = NULL;
#else
uint8_t* buffer_seq = NULL;
#endif

#ifdef UNPIPATH_OFF_K20
uint64_t* buffer_seqf = NULL;
uint64_t* buffer_off_g = NULL;
uint64_t* buffer_p = NULL;
uint64_t* buffer_pp = NULL;
uint64_t* buffer_hash_g = NULL;
uint64_t* chr_end_n = NULL;
uint64_t* chr_end_n_rec = NULL;
uint64_t first_cnt_cl = 0;
#else
uint32_t* buffer_seqf = NULL;
uint32_t* buffer_off_g = NULL;
uint32_t* buffer_p = NULL;
uint32_t* buffer_pp = NULL;
uint32_t* buffer_hash_g = NULL;
uint32_t* chr_end_n = NULL;
uint32_t* chr_end_n_rec = NULL;
uint32_t first_cnt_cl = 0;
#endif

//uint8_t* buffer_edge = NULL;

uint32_t tol_s_n = 0;
uint32_t tol_u = 0;

uint32_t* buffer_pupos = NULL;
uint32_t* buffer_puid = NULL;
uint32_t* buffer_kmer_g = NULL;
	
uint64_t result_ref_seq = 0;
uint64_t result_seq = 0;
uint64_t result_seqf = 0;
uint64_t result_edge = 0;
uint64_t result_p = 0;
uint64_t result_pp = 0;
uint64_t result_pu = 0;
uint64_t result_hash_g = 0;
uint64_t result_kmer_g = 0;
uint64_t result_off_g = 0;
uint64_t result_ref_g = 0;
uint64_t ref_length_p = 0;

uint8_t* aseq_array = NULL;
uint8_t* type_array = NULL;
uint32_t* rpos_array = NULL;
uint32_t* aindex_array = NULL;
uint32_t* apos_array = NULL;
uint32_t* rpos_end_array = NULL;
uint64_t* ref_seqv_array = NULL;
uint32_t* chr_rpos_a = NULL;

uint8_t* aseq_array_rec = NULL;
uint8_t* type_array_rec = NULL;
uint32_t* rpos_array_rec = NULL;
uint32_t* aindex_array_rec = NULL;
uint32_t* apos_array_rec = NULL;
uint32_t* rpos_end_array_rec = NULL;
uint64_t* ref_seqv_array_rec = NULL;
uint32_t* chr_rpos_a_rec = NULL;

const char* refseqs = "ref.seq";
const char* uniseqs = "unipath.seq";
const char* uniseq_bs = "unipath.seqb";
const char* uniseqf_bs = "unipath.seqfb";
const char* uniedges = "unipath.edge";
const char* unipus = "unipath.pu";
const char* uniposs = "unipath.pos";
const char* uniposps = "unipath.posp";
//const char* unistas = "unipath.sta";
const char* unihash_gs = "unipath_g.hash";
const char* unikmer_gs = "unipath_g.kmer";
const char* unioff_gs = "unipath_g.offset";
const char* unichrs = "unipath.chr";
const char* f_n = "N.sta";
const char* graph = "graph.sta";
const char* divs = "div/";
const char* f_size = "unipath.size";
const char* sys_c_mkdir = "mkdir ";
const char* sys_c_rm = "rm -rf ";
const char* sys_rm = "rm ";

const char* vcf_rposs = "vcf_rpos";
const char* vcf_aindexs = "vcf_aindex";
const char* vcf_aposs = "vcf_apos";
const char* vcf_aseqs = "vcf_aseq";
const char* vcf_rpos_ends = "vcf_rpos_end";
const char* ref_vcfs = "ref_vcf.seq";
const char* vcf_sizes = "vcf_size";
const char* vcf_types = "vcf_type";
const char* vcf_chr_rposs = "vcf_chr_rpos";

const char* vcf_rposs_rev = "vcf_rpos_rev";
const char* vcf_aindexs_rev = "vcf_aindex_rev";
const char* vcf_aposs_rev = "vcf_apos_rev";
const char* vcf_aseqs_rev = "vcf_aseq_rev";
const char* vcf_rpos_ends_rev = "vcf_rpos_end_rev";
const char* ref_vcfs_rev = "ref_vcf_rev.seq";
const char* vcf_sizes_rev = "vcf_size_rev";
const char* vcf_types_rev = "vcf_type_rev";

char vcf_rpos[ROUTE_LENGTH_MAX];
char vcf_aindex[ROUTE_LENGTH_MAX];
char vcf_apos[ROUTE_LENGTH_MAX];
char vcf_aseq[ROUTE_LENGTH_MAX];
char vcf_rpos_end[ROUTE_LENGTH_MAX];
char ref_vcf[ROUTE_LENGTH_MAX];
char vcf_size[ROUTE_LENGTH_MAX];
char vcf_type[ROUTE_LENGTH_MAX];
char vcf_chr_rpos[ROUTE_LENGTH_MAX];

char vcf_rpos_rev[ROUTE_LENGTH_MAX];
char vcf_aindex_rev[ROUTE_LENGTH_MAX];
char vcf_apos_rev[ROUTE_LENGTH_MAX];
char vcf_aseq_rev[ROUTE_LENGTH_MAX];
char vcf_rpos_end_rev[ROUTE_LENGTH_MAX];
char ref_vcf_rev[ROUTE_LENGTH_MAX];
char vcf_size_rev[ROUTE_LENGTH_MAX];
char vcf_type_rev[ROUTE_LENGTH_MAX];


char sam_result[ROUTE_LENGTH_MAX];
char read_fastq1[ROUTE_LENGTH_MAX];
char read_fastq2[ROUTE_LENGTH_MAX];
char index_route[ROUTE_LENGTH_MAX];
char filename_ref[ROUTE_LENGTH_MAX];
char filename_div[ROUTE_LENGTH_MAX];
char filename_sta[ROUTE_LENGTH_MAX];
char ref_seq[ROUTE_LENGTH_MAX];
char uniseq[ROUTE_LENGTH_MAX];
char uniseq_b[ROUTE_LENGTH_MAX];
char uniseqf_b[ROUTE_LENGTH_MAX];
char uniedge[ROUTE_LENGTH_MAX];
char unipu[ROUTE_LENGTH_MAX];
char unipos[ROUTE_LENGTH_MAX];
char uniposp[ROUTE_LENGTH_MAX];
//char unista[ROUTE_LENGTH_MAX];
char unihash_g[ROUTE_LENGTH_MAX];
char unikmer_g[ROUTE_LENGTH_MAX];
char unioff_g[ROUTE_LENGTH_MAX];
char unichr[ROUTE_LENGTH_MAX];
char N_route[ROUTE_LENGTH_MAX];
char unisize[ROUTE_LENGTH_MAX];
char vcf_file_name[ROUTE_LENGTH_MAX];
char chr_names[MAX_CHR_NUM][MAX_CHR_NAME_LENGTH];
char chr_line_content[MAX_CHR_NAME_LENGTH];

//files
FILE* fp_ref_seq = NULL;
FILE* fp_us_b = NULL;
FILE* fp_usf_b = NULL;
FILE* fp_ub = NULL;
//FILE* fp_ue = NULL;
FILE* fp_up = NULL;
FILE* fp_upp = NULL;
FILE* fp_pu = NULL;
FILE* fp_sta = NULL;
FILE* fp_chr = NULL;
FILE* fp_uh = NULL;
FILE* fp_uf = NULL;
FILE* fp_hash = NULL;
FILE* fp_kmer = NULL;
FILE* fp_off = NULL;
FILE* fp_n = NULL;
FILE* fp_us = NULL;
FILE* fp_num = NULL;
FILE* unipath_debug = NULL;

int** chr_res_buffer = NULL;
char** xa_d1s_buffer = NULL;
char** xa_d2s_buffer = NULL;
uint32_t** sam_pos1s_buffer = NULL;
uint32_t** sam_pos2s_buffer = NULL;
int** lv_re1s_buffer = NULL;
int** lv_re2s_buffer = NULL;

const uint8_t f = 4;
const uint8_t k = 14;
uint8_t k_t = 22;	
uint32_t chr_file_n = 1;
uint8_t thread_n = 1;
uint32_t upper_ins = 700;
uint32_t floor_ins = 300;
uint16_t readlen_max = 512;//252
uint16_t seed_l_max = 0;
uint8_t seed_l_l = 5;
uint16_t seed_step = 5; 
uint8_t cir_fix_n = 4;
uint16_t pos_n_max = 300;
uint16_t length_reduce = 22;
float va_ra = 0.02;

//important parameter 
float mis_match_r = 0.04;
float lv_rate = 0.06;//0.1
float last_circle_rate = 0;

//
float max_pair_score_r = 0.05;//0.06 0.2 
float mis_match_r_single = 0.05;//0.04 0.4
float lv_rate_anchor = 0.05;

uint8_t anchor_ext_flag = 0;
uint8_t flag_std = 0;
uint8_t local_ksw = 0;
uint8_t mgn_flag = 1;
uint8_t tree_flag = 0;

uint16_t cus_ali_n = 20; //150
uint16_t cus_max_output_ali = 150;//300 anchor
float max_single_score_r = 0.06;
uint16_t pr_single_outputn = 1000;
uint16_t seed_filter_pos_numn = 100;
uint16_t seed_filter_pos_num_singlen = 100;//

int8_t mat_score = 1, mis_score = 4, gapo_score = 6, gape_score = 1;

static int index_build_usage()
{
	fprintf(stderr, "\n");
	fprintf(stderr, "Program: de Brijn Graph-based mapping system index building\n");
	fprintf(stderr, "Version: %s\n", PACKAGE_VERSION);
	fprintf(stderr, "Contact: Hongzhe Guo <hzguo@hit.edu>\n\n");
	fprintf(stderr, "Usage:   deBGA index [options] reference.fasta <index_route> \n\n");
	fprintf(stderr, "Options: -k	INT      the k-mer length of the vertices of RdBG [20-28]\n");
	fprintf(stderr, "	 --ext-alt	STR (default: not set) vcf-info-based alignment\n");
	fprintf(stderr, "\n");
	return 1;
}

static int load_input_usage()
{
	fprintf(stderr, "\n");
	fprintf(stderr, "Program:	de Brijn Graph-based mapping system seed reduction and alignment\n");
	fprintf(stderr, "Version:	%s\n", PACKAGE_VERSION);
	fprintf(stderr, "Contact:	Hongzhe Guo <hzguo@hit.edu>\n\n");
	fprintf(stderr, "Usage:  	deBGA aln [options] <index_route> <read pair-end1.fq> [read pair-end2.fq] <result_file.sam>\n\nOptions:\n");
	fprintf(stderr, "	-k INT	the minimum length of a valid Uni-MEM seed [21-28]\n");	
	fprintf(stderr, "	-s INT	the number of iterations of re-seeding [%u]\n", cir_fix_n);
	fprintf(stderr, "	-i INT	the minimum interval of seeding [%u]\n", seed_step);
	fprintf(stderr, "	-n INT	the maximum allowed number of hits per seed [%u]\n", pos_n_max);
	fprintf(stderr, "	-c NUM	the threshold on the edit distance for early stop [%.2f]\n", max_pair_score_r);
	fprintf(stderr, "	--cl NUM the adjusted threshold on the edit distance for early stop [%.2f]\n", last_circle_rate);
	fprintf(stderr, "	--local  the local alignment option for confident alignment\n");
	fprintf(stderr, "	--local-match NUM the score for a matched base in the local alignment [%d]\n", mat_score);
	fprintf(stderr, "	--local-mismatch NUM the penalty for a mismatched base in the local alignment [%d]\n", mis_score);
	fprintf(stderr, "	--local-gap-open NUM the penalty for a gap open in the local alignment [%d]\n", gapo_score);
	fprintf(stderr, "	--local-gap-extension NUM the penalty for gap extension in the local alignment [%d]\n", gape_score);
	fprintf(stderr, "	--ext-alt-aln  (default: not set) vcf-info-based alignment\n");
	fprintf(stderr, "	--stdout   (default: not set) output alignments by stdout\n");
	fprintf(stderr, "	-u INT	the upper limit of insert size (only for pair-end reads) [%u] \n", upper_ins);
	fprintf(stderr, "	-f INT	the lower limit of insert size (only for pair-end reads) [%u] \n", floor_ins);
	fprintf(stderr, "	-o INT	the maximum number of alignment output [%u]\n", cus_ali_n);
	fprintf(stderr, "	-x INT	the maximum number of alignment output for anchoring alignment [%u]\n", cus_max_output_ali);
	fprintf(stderr, "	-l INT	the maximum allowed read length [%u]\n", readlen_max);
	fprintf(stderr, "	-e INT	the budget for single-end alignment [%u]\n", seed_filter_pos_numn);
	fprintf(stderr, "	-p INT	the number of threads [%u]\n", thread_n);
	//fprintf(stderr, "	--mg 	use the mode of multi-genomes\n");

	fprintf(stderr, "	Please refer to the following link for more detailed information about the options: https://github.com/HIT-Bioinformatics/deBGA\n");
	fprintf(stderr, "\n");

	return 1;
}

int help_usage()
{
	fprintf(stderr, "\n");
	fprintf(stderr, "Program:	deBGA (De bruijn graph nucleotide alignment)\n");
	fprintf(stderr, "Usage:  	deBGA <command> [options]\n\n");
	fprintf(stderr, "Command:	index		index sequences in the FASTA format\n");
	fprintf(stderr, "		aln      	pair-end and single-end reads seed reduction and alignment based on De bruijn graph\n");
	fprintf(stderr, "\n");
	fprintf(stderr, "Usage:   deBGA index [options] reference.fasta <index_route> \n\n");
	fprintf(stderr, "Options: -k INT      the k-mer length of the vertices of RdBG [20-28]\n");
	fprintf(stderr, "\n");
	fprintf(stderr, "Usage:  	deBGA aln [options] <index_route> <read pair-end1.fq> [read pair-end2.fq] <result_file.sam>\n\nOptions:\n");
	fprintf(stderr, "	-k INT	the minimum length of a valid Uni-MEM seed [21-28]\n");	
	fprintf(stderr, "	-s INT	the number of iterations of re-seeding [%u]\n", cir_fix_n);
	fprintf(stderr, "	-i INT	the minimum interval of seeding [%u]\n", seed_step);
	fprintf(stderr, "	-n INT	the maximum allowed number of hits per seed [%u]\n", pos_n_max);
	fprintf(stderr, "	-c NUM	the threshold on the edit distance for early stop [%.2f]\n", max_pair_score_r);
	fprintf(stderr, "	--cl NUM the adjusted threshold on the edit distance for early stop [%.2f]\n", last_circle_rate);
	fprintf(stderr, "	--local  the local alignment option for confident alignment\n");
	fprintf(stderr, "	--local-match NUM the score for a matched base in the local alignment [%d]\n", mat_score);
	fprintf(stderr, "	--local-mismatch NUM the penalty for a mismatched base in the local alignment [%d]\n", mis_score);
	fprintf(stderr, "	--local-gap-open NUM the penalty for a gap open in the local alignment [%d]\n", gapo_score);
	fprintf(stderr, "	--local-gap-extension NUM the penalty for gap extension in the local alignment [%d]\n", gape_score);
	fprintf(stderr, "	--anchor-ext (default: not set) deBGA do anchor alignment if the tree-extension do not output an alignment record\n");
	fprintf(stderr, "	--stdout   (default: not set) output alignments by stdout\n");
	fprintf(stderr, "	-u INT	the upper limit of insert size (only for pair-end reads) [%u] \n", upper_ins);
	fprintf(stderr, "	-f INT	the lower limit of insert size (only for pair-end reads) [%u] \n", floor_ins);
	fprintf(stderr, "	-o INT	the maximum number of alignment output [%u]\n", cus_ali_n);
	fprintf(stderr, "	-x INT	the maximum number of alignment output for anchoring alignment [%u]\n", cus_max_output_ali);
	fprintf(stderr, "	-l INT	the maximum allowed read length [%u]\n", readlen_max);
	fprintf(stderr, "	-e INT	the budget for single-end alignment [%u]\n", seed_filter_pos_numn);
	fprintf(stderr, "	-p INT	the number of threads [%u]\n", thread_n);
	//fprintf(stderr, "	--mg 	use the mode of multi-genomes\n");

	fprintf(stderr, "	Please refer to the following link for more detailed information about the options: https://github.com/hitbc/deBGA2\n");
	fprintf(stderr, "\n");

	
	return 1;
}

enum {
	PAR_LOCAL_KSW,
	PAR_LAST_CIRCLE,
	PAR_MULTI_GENOMES,
	PAR_MAT_SCORE,
	PAR_MIS_SCORE,
	PAR_GAPO_SCORE,
	PAR_GAPE_SCORE,
	PAR_A_ANCHOR,
	PAR_LV_ANCHOR,
	PAR_EXT_ALT,
	PAR_EXT_ALT_ALN,
	PAR_HELP,
	PAR_ANCHOR_EXT,
	PAR_STD_OUT
};

static const char *short_option = "k:p:u:f:l:r:i:s:n:o:x:a:c:g:e:v:";
static struct option long_option[] = {
	{(char*)"local", no_argument,  0, PAR_LOCAL_KSW},
	{(char*)"cl", required_argument, 0, PAR_LAST_CIRCLE},
	{(char*)"mg", no_argument, 0, PAR_MULTI_GENOMES},
	{(char*)"local-match", required_argument, 0, PAR_MAT_SCORE},
	{(char*)"local-mismatch", required_argument, 0, PAR_MIS_SCORE},
	{(char*)"local-gap-open", required_argument, 0, PAR_GAPO_SCORE},
	{(char*)"local-gap-extension", required_argument, 0, PAR_GAPE_SCORE},
	{(char*)"aanchor", required_argument, 0, PAR_A_ANCHOR},
	{(char*)"vanchor", required_argument, 0, PAR_LV_ANCHOR},
	{(char*)"ext-alt", required_argument, 0, PAR_EXT_ALT},
	{(char*)"help", no_argument, 0, PAR_HELP},
	{(char*)"anchor-ext", no_argument, 0, PAR_ANCHOR_EXT},
	{(char*)"stdout", no_argument, 0, PAR_STD_OUT},
	{(char*)"ext-alt-aln", no_argument, 0, PAR_EXT_ALT_ALN},
	{(char*)0, 0, 0, 0}
};

int load_input_index(int argc, char *argv[])
{
	int c = 0;

	//while ((c = getopt(argc, argv, "k:")) >= 0) {
	while((c = getopt_long(argc, argv, short_option, long_option, NULL)) != -1){
		switch (c) {
		case 'k': k_t = atoi(optarg); break;
		case PAR_EXT_ALT: strcpy(vcf_file_name, optarg); tree_flag = 1; break;
		case PAR_HELP: return index_build_usage();
		//default: return index_build_usage();
		}
	}
	
	if (argc - optind < 2) return index_build_usage();
	
	memset(filename_ref, 0, ROUTE_LENGTH_MAX);
	memset(index_route, 0, ROUTE_LENGTH_MAX);
	
	strcpy(filename_ref, argv[optind]);
	strcpy(index_route, argv[optind + 1]);

	if(index_route[strlen(index_route) - 1] != '/')	strcat(index_route, "/");
	
	fprintf(stderr, "reference file: %s\nindex route: %s\n", filename_ref, index_route);
	
	memset(ref_seq, 0, ROUTE_LENGTH_MAX);
	strcpy(ref_seq, index_route);
	strcat(ref_seq, refseqs);
	
	memset(uniseq, 0, ROUTE_LENGTH_MAX);
	strcpy(uniseq, index_route);
	strcat(uniseq, uniseqs);
	
	memset(uniseq_b, 0, ROUTE_LENGTH_MAX);
	strcpy(uniseq_b, index_route);
	strcat(uniseq_b, uniseq_bs);
	
	memset(uniseqf_b, 0, ROUTE_LENGTH_MAX);
	strcpy(uniseqf_b, index_route);
	strcat(uniseqf_b, uniseqf_bs);
	
	memset(uniedge, 0, ROUTE_LENGTH_MAX);
	strcpy(uniedge, index_route);
	strcat(uniedge, uniedges);
	
	memset(unipu, 0, ROUTE_LENGTH_MAX);
	strcpy(unipu, index_route);
	strcat(unipu, unipus);
	
	memset(unipos, 0, ROUTE_LENGTH_MAX);
	strcpy(unipos, index_route);
	strcat(unipos, uniposs);
	
	memset(uniposp, 0, ROUTE_LENGTH_MAX);
	strcpy(uniposp, index_route);
	strcat(uniposp, uniposps);
	/*
	memset(unista, 0, ROUTE_LENGTH_MAX);
	strcpy(unista, index_route);
	strcat(unista, unistas);
	*/
	memset(unihash_g, 0, ROUTE_LENGTH_MAX);
	strcpy(unihash_g, index_route);
	strcat(unihash_g, unihash_gs);
	
	memset(unikmer_g, 0, ROUTE_LENGTH_MAX);
	strcpy(unikmer_g, index_route);
	strcat(unikmer_g, unikmer_gs);
	
	memset(unioff_g, 0, ROUTE_LENGTH_MAX);
	strcpy(unioff_g, index_route);
	strcat(unioff_g, unioff_gs);
	
	memset(unichr, 0, ROUTE_LENGTH_MAX);
	strcpy(unichr, index_route);
	strcat(unichr, unichrs);
	
	memset(N_route, 0, ROUTE_LENGTH_MAX);
	strcpy(N_route, index_route);
	strcat(N_route, f_n);
	
	memset(filename_div, 0, ROUTE_LENGTH_MAX);
	strcpy(filename_div, index_route);
	strcat(filename_div, divs);
	
	memset(filename_sta, 0, ROUTE_LENGTH_MAX);
	strcpy(filename_sta, filename_div);
	strcat(filename_sta, graph);
	
	memset(unisize, 0, ROUTE_LENGTH_MAX);
	strcpy(unisize, index_route);
	strcat(unisize, f_size);

	if(tree_flag)
	{
		fprintf(stderr, "use the VCF-info-based extension and VCF file is %s\n", vcf_file_name);
		
		memset(vcf_rpos, 0, ROUTE_LENGTH_MAX);
		strcpy(vcf_rpos, index_route);
		strcat(vcf_rpos, vcf_rposs);
		
		memset(vcf_aindex, 0, ROUTE_LENGTH_MAX);
		strcpy(vcf_aindex, index_route);
		strcat(vcf_aindex, vcf_aindexs);
		
		memset(vcf_apos, 0, ROUTE_LENGTH_MAX);
		strcpy(vcf_apos, index_route);
		strcat(vcf_apos, vcf_aposs);
		
		memset(vcf_aseq, 0, ROUTE_LENGTH_MAX);
		strcpy(vcf_aseq, index_route);
		strcat(vcf_aseq, vcf_aseqs);
		
		memset(vcf_rpos_end, 0, ROUTE_LENGTH_MAX);
		strcpy(vcf_rpos_end, index_route);
		strcat(vcf_rpos_end, vcf_rpos_ends);
		
		memset(ref_vcf, 0, ROUTE_LENGTH_MAX);
		strcpy(ref_vcf, index_route);
		strcat(ref_vcf, ref_vcfs);
		
		memset(vcf_size, 0, ROUTE_LENGTH_MAX);
		strcpy(vcf_size, index_route);
		strcat(vcf_size, vcf_sizes);
		
		memset(vcf_type, 0, ROUTE_LENGTH_MAX);
		strcpy(vcf_type, index_route);
		strcat(vcf_type, vcf_types);
		
		memset(vcf_rpos_rev, 0, ROUTE_LENGTH_MAX);
		strcpy(vcf_rpos_rev, index_route);
		strcat(vcf_rpos_rev, vcf_rposs_rev);
		
		memset(vcf_aindex_rev, 0, ROUTE_LENGTH_MAX);
		strcpy(vcf_aindex_rev, index_route);
		strcat(vcf_aindex_rev, vcf_aindexs_rev);
		
		memset(vcf_apos_rev, 0, ROUTE_LENGTH_MAX);
		strcpy(vcf_apos_rev, index_route);
		strcat(vcf_apos_rev, vcf_aposs_rev);
		
		memset(vcf_aseq_rev, 0, ROUTE_LENGTH_MAX);
		strcpy(vcf_aseq_rev, index_route);
		strcat(vcf_aseq_rev, vcf_aseqs_rev);
		
		memset(vcf_rpos_end_rev, 0, ROUTE_LENGTH_MAX);
		strcpy(vcf_rpos_end_rev, index_route);
		strcat(vcf_rpos_end_rev, vcf_rpos_ends_rev);
		
		memset(ref_vcf_rev, 0, ROUTE_LENGTH_MAX);
		strcpy(ref_vcf_rev, index_route);
		strcat(ref_vcf_rev, ref_vcfs_rev);
		
		memset(vcf_size_rev, 0, ROUTE_LENGTH_MAX);
		strcpy(vcf_size_rev, index_route);
		strcat(vcf_size_rev, vcf_sizes_rev);
		
		memset(vcf_type_rev, 0, ROUTE_LENGTH_MAX);
		strcpy(vcf_type_rev, index_route);
		strcat(vcf_type_rev, vcf_types_rev);
		
		memset(vcf_chr_rpos, 0, ROUTE_LENGTH_MAX);
		strcpy(vcf_chr_rpos, index_route);
		strcat(vcf_chr_rpos, vcf_chr_rposs);
	}

	chr_end_n = (uint64_t* )calloc(MAX_CHR_NUM, 8);
#ifdef	HANDLE_DIR
	char create_route[ROUTE_LENGTH_MAX];
	char rm_route[ROUTE_LENGTH_MAX];
	DIR* directory_pointer = NULL;
#ifndef	DEBUG_VCF
	if((directory_pointer = opendir(index_route)) != NULL)//filename_div
	{
		memset(rm_route, 0, ROUTE_LENGTH_MAX);
		strcpy(rm_route, sys_c_rm);
		strcat(rm_route, index_route);

		system(rm_route);
	}
#endif	
	memset(create_route, 0, ROUTE_LENGTH_MAX);
	strcpy(create_route, sys_c_mkdir);
	strcat(create_route, index_route);

	if((directory_pointer = opendir(index_route))==NULL)
		system(create_route);
	
	memset(create_route, 0, ROUTE_LENGTH_MAX);
	strcpy(create_route, sys_c_mkdir);
	strcat(create_route, filename_div);
	
	if((directory_pointer = opendir(filename_div))==NULL)
		system(create_route);
#endif
		
	return 0;
}

int load_input_map(int argc, char *argv[])
{
	int c = 0;
	while((c = getopt_long(argc, argv, short_option, long_option, NULL)) != -1){
		switch (c) {
			case 'k': k_t = atoi(optarg); break;
			case 'p': thread_n = atoi(optarg); break;
			case 'u': upper_ins = atoi(optarg); break;
			case 'f': floor_ins = atoi(optarg); break;
			case 'l': readlen_max = atoi(optarg); break;
			case 'r': length_reduce = atoi(optarg); break;
			case 'i': seed_step = atoi(optarg); break;
			case 's': cir_fix_n = atoi(optarg); break;
			case 'n': pos_n_max = atoi(optarg); break;
			case 'o': cus_ali_n = atoi(optarg); break;
			case 'x': cus_max_output_ali = atoi(optarg); break;
			case 'c': max_pair_score_r = atof(optarg); lv_rate_anchor = max_pair_score_r; mis_match_r_single = max_pair_score_r; break;
			case 'e': seed_filter_pos_numn = atoi(optarg); break;
			case 'v': lv_rate = atof(optarg); break;
			case PAR_A_ANCHOR: mis_match_r_single = atof(optarg); break;
			case PAR_LV_ANCHOR: lv_rate_anchor = atof(optarg); break;
			case PAR_LAST_CIRCLE: last_circle_rate = atof(optarg); mis_match_r_single = last_circle_rate; break;
			case PAR_LOCAL_KSW: local_ksw = 1; last_circle_rate = 0; mis_match_r_single = 0.4; break;
			case PAR_MULTI_GENOMES: mgn_flag = 0;break;
			case PAR_MAT_SCORE: mat_score = atoi(optarg); break;
			case PAR_MIS_SCORE: mis_score = atoi(optarg); break;
			case PAR_GAPO_SCORE: gapo_score = atoi(optarg); break;
			case PAR_GAPE_SCORE: gape_score = atoi(optarg); break;
			case PAR_ANCHOR_EXT: anchor_ext_flag = 1; break;
			case PAR_STD_OUT: flag_std = 1; break;
			case PAR_EXT_ALT_ALN: tree_flag = 1; break;
			case PAR_HELP: return load_input_usage();
			default: return load_input_usage();
		} 
	}   

	if (argc - optind < 3) return load_input_usage();
	if((k_t < 21) || (k_t > 28))
	{
		fprintf(stderr, "Input error: -k cannot be less than 21 or more than 28\n");
		exit(1);
	}
	if((thread_n < 1) || (thread_n > 32))
	{
		fprintf(stderr, "Input error: -p cannot be less than 1 or more than 32\n");
		exit(1);
	}
	if(upper_ins <= floor_ins)
	{
		fprintf(stderr, "Input error: -u should be more than floor limit\n");
		exit(1);
	}
	if(readlen_max > 2048)
	{
		fprintf(stderr, "Input error: -l cannot be more than 2048\n");
		exit(1);
	}
	if(length_reduce > 50)
	{
		fprintf(stderr, "Input error: -r cannot be more than 50\n");
		exit(1);
	}
	if((seed_step < 1) || (seed_step > 20))//5 20
	{
		fprintf(stderr, "Input error: -i cannot be more than 20 or less than 5\n");
		exit(1);
	}
	if(cir_fix_n > 10)
	{
		fprintf(stderr, "Input error: -s cannot be more than 10\n");
		exit(1);
	}
	if(cus_ali_n > 1000)//300
	{
		fprintf(stderr, "Input error: -o cannot be more than 1000\n");//300
		exit(1);
	}
	if(cus_max_output_ali > 1000)//500
	{
		fprintf(stderr, "Input error: -x cannot be more than 1000\n");//500
		exit(1);
	}
	if(max_pair_score_r > 0.35)//0.5
	{
		fprintf(stderr, "Input error: -c cannot be more than 0.5\n");
		exit(1);
	}
	if(seed_filter_pos_numn > 2000)
	{
		fprintf(stderr, "Input error: -e cannot be more than 2000\n");
		exit(1);
	}
	if(lv_rate > 0.5)
	{
		fprintf(stderr, "Input error: -v cannot be more than 0.5\n");
		exit(1);
	}
	if(last_circle_rate > 0.5)
	{
		fprintf(stderr, "Input error: --cl cannot be more than 0.5\n");
		exit(1);
	}
	
	memset(index_route, 0, ROUTE_LENGTH_MAX);
	strcpy(index_route, argv[optind]);
	
	if(index_route[strlen(index_route) - 1] != '/')	strcat(index_route, "/");
	
	memset(read_fastq1, 0, ROUTE_LENGTH_MAX);
	strcpy(read_fastq1, argv[optind + 1]);

	memset(ref_seq, 0, ROUTE_LENGTH_MAX);
	strcpy(ref_seq, index_route);
	strcat(ref_seq, refseqs);
	
	memset(uniseq_b, 0, ROUTE_LENGTH_MAX);
	strcpy(uniseq_b, index_route);
	strcat(uniseq_b, uniseq_bs);
	
	memset(uniseqf_b, 0, ROUTE_LENGTH_MAX);
	strcpy(uniseqf_b, index_route);
	strcat(uniseqf_b, uniseqf_bs);
	
	memset(uniedge, 0, ROUTE_LENGTH_MAX);
	strcpy(uniedge, index_route);
	strcat(uniedge, uniedges);
	
	memset(unipu, 0, ROUTE_LENGTH_MAX);
	strcpy(unipu, index_route);
	strcat(unipu, unipus);
	
	memset(unipos, 0, ROUTE_LENGTH_MAX);
	strcpy(unipos, index_route);
	strcat(unipos, uniposs);
	
	memset(uniposp, 0, ROUTE_LENGTH_MAX);
	strcpy(uniposp, index_route);
	strcat(uniposp, uniposps);
	/*
	memset(unista, 0, ROUTE_LENGTH_MAX);
	strcpy(unista, index_route);
	strcat(unista, unistas);
	*/
	memset(unihash_g, 0, ROUTE_LENGTH_MAX);
	strcpy(unihash_g, index_route);
	strcat(unihash_g, unihash_gs);
	
	memset(unikmer_g, 0, ROUTE_LENGTH_MAX);
	strcpy(unikmer_g, index_route);
	strcat(unikmer_g, unikmer_gs);
	
	memset(unioff_g, 0, ROUTE_LENGTH_MAX);
	strcpy(unioff_g, index_route);
	strcat(unioff_g, unioff_gs);
	
	memset(unichr, 0, ROUTE_LENGTH_MAX);
	strcpy(unichr, index_route);
	strcat(unichr, unichrs);
	
	memset(N_route, 0, ROUTE_LENGTH_MAX);
	strcpy(N_route, index_route);
	strcat(N_route, f_n);
	
	memset(filename_div, 0, ROUTE_LENGTH_MAX);
	strcpy(filename_div, index_route);
	strcat(filename_div, divs);
	
	memset(filename_sta, 0, ROUTE_LENGTH_MAX);
	strcpy(filename_sta, filename_div);
	strcat(filename_sta, graph);
	
	memset(unisize, 0, ROUTE_LENGTH_MAX);
	strcpy(unisize, index_route);
	strcat(unisize, f_size);
	
	if(tree_flag)
	{
		fprintf(stderr, "use the VCF-info-based extension\n");
		
		memset(vcf_rpos, 0, ROUTE_LENGTH_MAX);
		strcpy(vcf_rpos, index_route);
		strcat(vcf_rpos, vcf_rposs);
		
		memset(vcf_aindex, 0, ROUTE_LENGTH_MAX);
		strcpy(vcf_aindex, index_route);
		strcat(vcf_aindex, vcf_aindexs);
		
		memset(vcf_apos, 0, ROUTE_LENGTH_MAX);
		strcpy(vcf_apos, index_route);
		strcat(vcf_apos, vcf_aposs);
		
		memset(vcf_aseq, 0, ROUTE_LENGTH_MAX);
		strcpy(vcf_aseq, index_route);
		strcat(vcf_aseq, vcf_aseqs);
		
		memset(vcf_rpos_end, 0, ROUTE_LENGTH_MAX);
		strcpy(vcf_rpos_end, index_route);
		strcat(vcf_rpos_end, vcf_rpos_ends);
		
		memset(ref_vcf, 0, ROUTE_LENGTH_MAX);
		strcpy(ref_vcf, index_route);
		strcat(ref_vcf, ref_vcfs);
		
		memset(vcf_size, 0, ROUTE_LENGTH_MAX);
		strcpy(vcf_size, index_route);
		strcat(vcf_size, vcf_sizes);
		
		memset(vcf_type, 0, ROUTE_LENGTH_MAX);
		strcpy(vcf_type, index_route);
		strcat(vcf_type, vcf_types);
		
		memset(vcf_rpos_rev, 0, ROUTE_LENGTH_MAX);
		strcpy(vcf_rpos_rev, index_route);
		strcat(vcf_rpos_rev, vcf_rposs_rev);
		
		memset(vcf_aindex_rev, 0, ROUTE_LENGTH_MAX);
		strcpy(vcf_aindex_rev, index_route);
		strcat(vcf_aindex_rev, vcf_aindexs_rev);
		
		memset(vcf_apos_rev, 0, ROUTE_LENGTH_MAX);
		strcpy(vcf_apos_rev, index_route);
		strcat(vcf_apos_rev, vcf_aposs_rev);
		
		memset(vcf_aseq_rev, 0, ROUTE_LENGTH_MAX);
		strcpy(vcf_aseq_rev, index_route);
		strcat(vcf_aseq_rev, vcf_aseqs_rev);
		
		memset(vcf_rpos_end_rev, 0, ROUTE_LENGTH_MAX);
		strcpy(vcf_rpos_end_rev, index_route);
		strcat(vcf_rpos_end_rev, vcf_rpos_ends_rev);
		
		memset(ref_vcf_rev, 0, ROUTE_LENGTH_MAX);
		strcpy(ref_vcf_rev, index_route);
		strcat(ref_vcf_rev, ref_vcfs_rev);
		
		memset(vcf_size_rev, 0, ROUTE_LENGTH_MAX);
		strcpy(vcf_size_rev, index_route);
		strcat(vcf_size_rev, vcf_sizes_rev);
		
		memset(vcf_type_rev, 0, ROUTE_LENGTH_MAX);
		strcpy(vcf_type_rev, index_route);
		strcat(vcf_type_rev, vcf_types_rev);
	
		memset(vcf_chr_rpos, 0, ROUTE_LENGTH_MAX);
		strcpy(vcf_chr_rpos, index_route);
		strcat(vcf_chr_rpos, vcf_chr_rposs);
	}
	
	if ((argc - optind < 4) && (flag_std == 0))
	{
		max_single_score_r = max_pair_score_r;
		seed_filter_pos_num_singlen = seed_filter_pos_numn;

		if(strcmp(argv[optind + 2] + strlen(argv[optind + 2]) - 3, ".fq") || strcmp(argv[optind + 2] + strlen(argv[optind + 2]) - 6, ".fastq"))
		{
			fprintf(stderr, "result file cannot be .fq or .fastq\n");
			exit(1);
		}
		strcpy(sam_result, argv[optind + 2]);
		seed_ali_single_end();
	}else{
		memset(read_fastq2, 0, ROUTE_LENGTH_MAX);
		strcpy(read_fastq2, argv[optind + 2]);
		if(flag_std == 0)
			strcpy(sam_result, argv[optind + 3]);
		seed_ali();
	}	

	return 0;
}

void load_index_file()
{
	uint64_t a_size = 0;
	uint64_t us_n = 0;
	uint64_t usf_n = 0;
	uint64_t ue_n = 0;
	uint64_t up_n = 0;
	uint64_t upp_n = 0;
	uint64_t hash_n = 0;
	uint64_t kmer_n = 0;
	uint64_t off_n = 0;
	uint64_t pu_n = 0;
	uint64_t ref_seq_n = 0;
	uint64_t file_size = 0;
	
	//read index file size
	fp_num = fopen(unisize, "rb");
	if (fp_num == NULL)
    {
        fputs ("File error opening the size of index file: wrong index route or wrong index file name\n",stderr);
        exit (1);
    }
	fread(&us_n, 8, 1, fp_num);
	fread(&usf_n, 8, 1, fp_num);
	fread(&ue_n, 8, 1, fp_num);
	fread(&up_n, 8, 1, fp_num);
	fread(&upp_n, 8, 1, fp_num);
	fread(&hash_n, 8, 1, fp_num);
	fread(&kmer_n, 8, 1, fp_num);
	fread(&off_n, 8, 1, fp_num);
	fread(&pu_n, 8, 1, fp_num);
	fread(&ref_seq_n, 8, 1, fp_num);

	if(tree_flag)	fread(&ref_length_p, 8, 1, fp_num);
	
	fclose(fp_num);
	
	//read ref seq file
	fprintf(stderr, "Load ref seq\n");
	
	fp_ref_seq = fopen(ref_seq, "rb");
	if (fp_ref_seq == NULL)
    {
        fputs ("File error opening the  seq file\n",stderr);
        exit (1);
    }
	
	fseek(fp_ref_seq, 0, SEEK_END);// non-portable
    file_size = ftell(fp_ref_seq);
	rewind(fp_ref_seq);

	ref_seq_n = file_size;
	buffer_ref_seq = (uint64_t* )calloc(ref_seq_n + 536, 1);//536 = (2048 >> 5 + 3) << 3
	
	a_size = ref_seq_n >> 3;
	result_ref_seq = fread(buffer_ref_seq, 8, a_size, fp_ref_seq);
	
	if (result_ref_seq != a_size)
    {
        fprintf(stderr, "Reading error");
        exit (3);
    }

    fclose(fp_ref_seq);

    //read input unipath seq file
    fprintf(stderr, "Load unipath seq\n");

    fp_us_b = fopen (uniseq_b, "rb" );
    if (fp_us_b == NULL)
    {
        fprintf(stderr, "File error opening the unipath seq file\n");
        exit (1);
    }

	fseek(fp_us_b, 0, SEEK_END);// non-portable
    file_size = ftell(fp_us_b);
	rewind(fp_us_b);

	us_n = file_size;
	
    // allocate memory to contain the whole file:
#ifdef UNI_SEQ64
	buffer_seq = (uint64_t* ) malloc (us_n);
#else	
    buffer_seq = (uint8_t* ) malloc (us_n);
#endif
    if (buffer_seq == NULL)
    {
        fputs ("Memory error buffer_seq\n",stderr);
        exit (2);
    }

    // copy the file into the buffer:
	
#ifdef UNI_SEQ64	
	a_size = us_n >> 3;
	result_seq = fread (buffer_seq, 8, a_size, fp_us_b);
#else	
	a_size = us_n;
    result_seq = fread (buffer_seq, 1, a_size, fp_us_b);
#endif

    if (result_seq != a_size)
    {
        fputs ("Reading error",stderr);
        exit (3);
    }

    fclose(fp_us_b);

    //read input unipath offset file
    fprintf(stderr, "Load unipath offset\n");

    fp_usf_b = fopen(uniseqf_b, "rb");
    if(fp_usf_b == NULL)
    {
        fputs ("File error opening the unipath offset file\n",stderr);
        exit (1);
    }

	fseek(fp_usf_b, 0, SEEK_END);// non-portable
    file_size = ftell(fp_usf_b);
	rewind(fp_usf_b);

	usf_n = file_size;
	
    // allocate memory to contain the whole file:
	// copy the file into the buffer:
#ifdef UNPIPATH_OFF_K20
	buffer_seqf = (uint64_t* ) malloc (usf_n);
	if (buffer_seqf == NULL)
    {
        fprintf(stderr, "Memory error buffer_seqf");
        exit (2);
    }
	a_size = (usf_n >> 3);
    result_seqf = fread (buffer_seqf, 8, a_size, fp_usf_b);
#else	
    buffer_seqf = (uint32_t* ) malloc (usf_n);
	if (buffer_seqf == NULL)
    {
        fputs ("Memory error buffer_seqf",stderr);
        exit (2);
    }
	a_size = (usf_n >> 2);
    result_seqf = fread (buffer_seqf, 4, a_size, fp_usf_b);
#endif	
	
    if (result_seqf != a_size)
    {
        fputs ("Reading error",stderr);
        exit (3);
    }

    fclose(fp_usf_b);

    //read input unipath position file
    fprintf(stderr, "Load unipath position\n");

    fp_up = fopen(unipos, "rb");
    if(fp_up == NULL)
    {
        fprintf(stderr, "File error opening the unipath position file\n");
        exit (1);
    }

	fseek(fp_up, 0, SEEK_END);// non-portable
    file_size = ftell(fp_up);
	rewind(fp_up);
	
	up_n = file_size;
	
    // allocate memory to contain the whole file:

	// copy the file into the buffer:
#ifdef	UNPIPATH_OFF_K20
	buffer_p = (uint64_t* ) malloc (up_n);
    if (buffer_p == NULL)
    {
        fputs ("Memory error buffer_p",stderr);
        exit (2);
    }
	a_size = (up_n >> 3);
	result_p = fread (buffer_p, 8, a_size, fp_up);
#else
	buffer_p = (uint32_t* ) malloc (up_n);
    if (buffer_p == NULL)
    {
        fputs ("Memory error buffer_p",stderr);
        exit (2);
    }
	a_size = (up_n >> 2);
	result_p = fread (buffer_p, 4, a_size, fp_up);
#endif
    
    if (result_p != a_size)
    {
        fputs ("Reading error",stderr);
        exit (3);
    }
	
    fclose(fp_up);

    //read input unipath position point file
    fprintf(stderr, "Load unipath position point\n");

    fp_upp = fopen(uniposp, "rb");
    if(fp_upp == NULL)
    {
        fputs ("File error opening the unipath position point file\n",stderr);
        exit (1);
    }
	
	fseek(fp_upp, 0, SEEK_END);// non-portable
    file_size = ftell(fp_upp);
	rewind(fp_upp);

	upp_n = file_size;
	
#ifdef	UNPIPATH_OFF_K20
	// allocate memory to contain the whole file:
    a_size = (upp_n >> 3);
    buffer_pp = (uint64_t* ) malloc (upp_n);
    if (buffer_pp == NULL)
    {
        fputs ("Memory error buffer_pp",stderr);
        exit (2);
    }
    // copy the file into the buffer:
    result_pp = fread (buffer_pp, 8, a_size, fp_upp);
#else
    // allocate memory to contain the whole file:
    a_size = (upp_n >> 2);
    buffer_pp = (uint32_t* ) malloc (upp_n);
    if (buffer_pp == NULL)
    {
        fputs ("Memory error buffer_pp",stderr);
        exit (2);
    }
    // copy the file into the buffer:
    result_pp = fread (buffer_pp, 4, a_size, fp_upp);
#endif
	if (result_pp != a_size)
    {
        fputs ("Reading error",stderr);
        exit (3);
    }
    fclose(fp_upp);
	
    //read input unipath hash file
    fprintf(stderr, "Load unipath hash\n");

    fp_hash = fopen(unihash_g, "rb");
    if(fp_hash == NULL)
    {
        fputs ("File error opening the graph hash file\n",stderr);
        exit (1);
    }
	
	fseek(fp_hash, 0, SEEK_END);// non-portable
    file_size = ftell(fp_hash);
	rewind(fp_hash);

	hash_n = file_size;
	
#ifdef	UNPIPATH_OFF_K20
	a_size = (hash_n >> 3);
    buffer_hash_g = (uint64_t* ) malloc (hash_n);
    if (buffer_hash_g == NULL)
    {
        fputs("Memory error",stderr);
        exit(2);
    }

    // copy the file into the buffer:
    result_hash_g = fread (buffer_hash_g, 8, a_size, fp_hash);
    if (result_hash_g != a_size)
    {
        fputs("Reading error",stderr);
        exit(3);
    }
#else
    a_size = (hash_n >> 2);
    buffer_hash_g = (uint32_t* ) malloc (hash_n);
    if (buffer_hash_g == NULL)
    {
        fputs("Memory error",stderr);
        exit(2);
    }

    // copy the file into the buffer:
    result_hash_g = fread (buffer_hash_g, 4, a_size, fp_hash);
    if (result_hash_g != a_size)
    {
        fputs("Reading error",stderr);
        exit(3);
    }
#endif

    fclose(fp_hash);

    //read input graph kmer file
    fprintf(stderr, "Load unipath kmer\n");

    fp_kmer = fopen(unikmer_g, "rb");
    if(fp_kmer == NULL)
    {
        fputs ("File error opening the graph hash file\n",stderr);
        exit (1);
    }
	
	fseek(fp_kmer, 0, SEEK_END);// non-portable
    file_size = ftell(fp_kmer);
	rewind(fp_kmer);

	kmer_n = file_size;
	
    a_size = (kmer_n >> 2);
    
    //a_size = kmer_num;

    buffer_kmer_g = (uint32_t* ) malloc (kmer_n);
    if (buffer_kmer_g == NULL)
    {
        fputs ("Memory error buffer_kmer_g",stderr);
        exit (2);
    }

    // copy the file into the buffer:
    result_kmer_g = fread (buffer_kmer_g, 4, a_size, fp_kmer);
    if (result_kmer_g != a_size)
    {
        fputs ("Reading error buffer_kmer_g",stderr);
        exit (3);
    }

    fclose(fp_kmer);

    //read input graph off file
    fprintf(stderr, "Load unipath off\n");

    fp_off = fopen(unioff_g, "rb");
    if(fp_off == NULL)
    {
        fprintf(stderr, "File error opening the graph hash file\n");
        exit (1);
    }
	
	fseek(fp_off, 0, SEEK_END);// non-portable
    file_size = ftell(fp_off);
	rewind(fp_off);

	off_n = file_size;

    // copy the file into the buffer:
#ifdef UNPIPATH_OFF_K20
	buffer_off_g = (uint64_t* ) malloc (off_n);
	if (buffer_off_g == NULL)
    {
        fputs ("Memory error",stderr);
        exit (2);
    }
	a_size = (off_n >> 3);
    result_off_g = fread (buffer_off_g, 8, a_size, fp_off);
#else	
    buffer_off_g = (uint32_t* ) malloc (off_n);
	if (buffer_off_g == NULL)
    {
        fputs ("Memory error",stderr);
        exit (2);
    }
	a_size = (off_n >> 2);
    result_off_g = fread (buffer_off_g, 4, a_size, fp_off);
#endif	

    if (result_off_g != a_size)
    {
        fputs ("Reading error",stderr);
        exit (3);
    }

    fclose(fp_off);

	fp_chr = fopen (unichr, "r" );
    if (fp_chr == NULL)
    {
        fputs ("File error opening the chr file\n",stderr);
        exit (1);
    }
	
	chr_end_n = (uint64_t* )calloc(MAX_CHR_NUM, 8);
	chr_end_n_rec = (uint64_t* )calloc(MAX_CHR_NUM, 8);

	uint32_t chr_line_n = 0;
    while(!feof(fp_chr))
    {
		fscanf(fp_chr, "%s", chr_line_content);
		//fgets(chr_line_content, MAX_CHR_NAME_LENGTH, fp_chr);

		if((chr_line_n & 0X1) == 0)
		{
			strcpy(chr_names[chr_file_n], chr_line_content);
			chr_names[chr_file_n][strlen(chr_names[chr_file_n])] = '\0';

		}else{
#ifdef UNPIPATH_OFF_K20
			sscanf(chr_line_content, "%"PRId64"", &chr_end_n[chr_file_n]);
#else
			sscanf(chr_line_content, "%u", &chr_end_n[chr_file_n]);
#endif
			chr_file_n++;
		}	
		
		fflush(stdout);
		
		chr_line_n++;
    }

	chr_end_n[0] = START_POS_REF + 1;
	strcpy(chr_names[chr_file_n], "*");

	if(tree_flag)
	{
		int chr_i = 0;
		chr_end_n_rec[0] = 0;
		for(chr_i = chr_file_n - 1, chr_line_n = 1; chr_i > 0 ; chr_i--, chr_line_n++)
			chr_end_n_rec[chr_line_n] = chr_end_n_rec[chr_line_n - 1] + chr_end_n[chr_i] - chr_end_n[chr_i - 1];

		chr_rpos_a = (uint32_t* ) calloc (chr_line_n + 1, 4);
		chr_rpos_a_rec = (uint32_t* ) calloc (chr_line_n + 1, 4);
	
		FILE* fp_chr_rpos = fopen(vcf_chr_rpos, "r");
		if (fp_chr_rpos == NULL)
		{
			fputs ("File error opening the vcf chr rpos file\n",stderr);
			exit (1);
		}

		fread (chr_rpos_a, 4, chr_line_n, fp_chr_rpos);
		fread (chr_rpos_a_rec, 4, chr_line_n, fp_chr_rpos);
	
		if(fp_chr_rpos)	fclose(fp_chr_rpos);	
	}
}


void load_input_sorted_vcf()
{
	uint32_t rpos_cnt = 0;
	uint32_t apos_cnt = 0;
	uint32_t ref_cv_n = 0;
	uint64_t aseq_cnt = 0;
	uint64_t rpos_rsize = 0;
	uint64_t apos_rsize = 0;
	uint64_t aindex_rsize = 0;
	uint64_t rpos_end_rsize = 0;
	uint64_t ref_seqv_rsize = 0;
	uint64_t aseq_rsize = 0;
	uint64_t type_rsize = 0;
	
	fprintf(stderr, "begin reading vcf info file\n");

	FILE* fp_in_rpos = fopen(vcf_rpos,"rb");
	if(fp_in_rpos == NULL)
	{
		fprintf(stderr, "Error of opening vcf rpos file\n");
		exit(1);
	}

	FILE* fp_in_aindex = fopen(vcf_aindex,"rb");
	if(fp_in_aindex == NULL)
	{
		fprintf(stderr, "Error of opening vcf ALT index file\n");
		exit(1);
	}

	FILE* fp_in_apos = fopen(vcf_apos,"rb");
	if(fp_in_apos == NULL)
	{
		fprintf(stderr, "Error of opening vcf ALT pos file\n");
		exit(1);
	}

	FILE* fp_in_aseq = fopen(vcf_aseq,"rb");
	if(fp_in_aseq == NULL)
	{
		fprintf(stderr, "Error of opening vcf aseq file\n");
		exit(1);
	}

	FILE* fp_in_rpos_end = fopen(vcf_rpos_end,"rb");
	if(fp_in_rpos_end == NULL)
	{
		fprintf(stderr, "Error of opening vcf rpos end file\n");
		exit(1);
	}

	FILE* fp_in_ref_seq_v = fopen(ref_vcf,"rb");
	if(fp_in_ref_seq_v == NULL)
	{
		fprintf(stderr, "Error of opening vcf ref seq file\n");
		exit(1);
	}

	FILE* fp_in_size = fopen(vcf_size,"rb");
	if(fp_in_size == NULL)
	{
		fprintf(stderr, "Error of opening vcf size file\n");
		exit(1);
	}

	FILE* fp_in_type = fopen(vcf_type,"rb");
	if(fp_in_type == NULL)
	{
		fprintf(stderr, "Error of opening vcf size file\n");
		exit(1);
	}

	fread(&rpos_cnt, 4, 1, fp_in_size);
	fread(&apos_cnt, 4, 1, fp_in_size);
	fread(&aseq_cnt, 4, 1, fp_in_size);
	fread(&ref_cv_n, 4, 1, fp_in_size);

	fprintf(stderr, "begin allocating memory for vcf info\n");

	rpos_array = (uint32_t* )calloc(rpos_cnt, 4);
	if(rpos_array == NULL)
	{
		fprintf(stderr, "Error of allocating memory rpos_array\n");
		exit(3);
	}
	aindex_array = (uint32_t* )calloc(rpos_cnt + 1, 4);
	if(aindex_array == NULL)
	{
		fprintf(stderr, "Error of allocating memory aindex_array\n");
		exit(3);
	}
	apos_array = (uint32_t* )calloc(apos_cnt + 1, 4);
	if(apos_array == NULL)
	{
		fprintf(stderr, "Error of allocating memory apos_array\n");
		exit(3);
	}
	rpos_end_array = (uint32_t* )calloc(apos_cnt + 1, 4);
	if(rpos_end_array == NULL)
	{
		fprintf(stderr, "Error of allocating memory rpos_end_array\n");
		exit(3);
	}
	aseq_array = (uint8_t* )calloc(aseq_cnt, 1);
	if(aseq_array == NULL)
	{
		fprintf(stderr, "Error of allocating memory aseq_array\n");
		exit(3);
	}
	ref_seqv_array = (uint64_t* )calloc((ref_cv_n >> 1) + 1000, 1);//(ref_cv_n >> 4 + 1) + 1000, 4
	if(ref_seqv_array == NULL)
	{
		fprintf(stderr, "Error of allocating memory ref_seqv_array\n");
		exit(3);
	}
	type_array = (uint8_t* )calloc(apos_cnt + 1, 1);
	if(type_array == NULL)
	{
		fprintf(stderr, "Error of allocating memory ref_seqv_array\n");
		exit(3);
	}

	if((rpos_rsize = fread(rpos_array, 4, rpos_cnt, fp_in_rpos)) != rpos_cnt)
	{
		fprintf (stderr,"Reading error of rpos file %u %u",rpos_cnt, rpos_rsize);
		exit(2);
	}

	if((aindex_rsize = fread(aindex_array, 4, rpos_cnt + 1, fp_in_aindex)) != (rpos_cnt + 1))
	{
		fputs ("Reading error of aindex file",stderr);
		exit(2);
	}

	if((apos_rsize = fread(apos_array, 4, apos_cnt + 1, fp_in_apos)) != (apos_cnt + 1))
	{
		fputs ("Reading error of apos file",stderr);
		exit(2);
	}

	if((rpos_end_rsize = fread(rpos_end_array, 4, apos_cnt, fp_in_rpos_end)) != apos_cnt)
	{
		fputs ("Reading error of rpos end file",stderr);
		exit(2);
	}

	if((ref_seqv_rsize = fread(ref_seqv_array, 8, (ref_cv_n >> 4) + 1, fp_in_ref_seq_v)) != ((ref_cv_n >> 4) + 1))
	{
		fputs ("Reading error of ref seq vcf file",stderr);
		exit(2);
	}

	if((aseq_rsize = fread(aseq_array, 1, aseq_cnt, fp_in_aseq)) != aseq_cnt)
	{
		fputs ("Reading error of aseq end file",stderr);
		exit(2);
	}

	if((type_rsize = fread(type_array, 1, apos_cnt, fp_in_type)) != apos_cnt)
	{
		fputs ("Reading error of vcf type file",stderr);
		exit(2);
	}

	//check whether ref and alt char is same to test fp_out_ref_seq_v

	if(fp_in_rpos)	fclose(fp_in_rpos);
	if(fp_in_aindex)	fclose(fp_in_aindex);
	if(fp_in_apos)	fclose(fp_in_apos);
	if(fp_in_rpos_end)	fclose(fp_in_rpos_end);
	if(fp_in_aseq)	fclose(fp_in_aseq);
	if(fp_in_ref_seq_v)	fclose(fp_in_ref_seq_v);
	if(fp_in_size)	fclose(fp_in_size);
	if(fp_in_type)	fclose(fp_in_type);

	fp_in_rpos = fopen(vcf_rpos_rev,"rb");
	if(fp_in_rpos == NULL)
	{
		printf("Error of opening vcf rpos rec file\n"); 
		exit(1);
	}		
	
	fp_in_aindex = fopen(vcf_aindex_rev,"rb");
	if(fp_in_aindex == NULL)
	{
		printf("Error of opening vcf ALT index rec file\n"); 
		exit(1);
	}
	
	fp_in_apos = fopen(vcf_apos_rev,"rb");
	if(fp_in_apos == NULL)
	{
		printf("Error of opening vcf ALT pos rec file\n"); 
		exit(1);
	}

	fp_in_aseq = fopen(vcf_aseq_rev,"rb");
	if(fp_in_aseq == NULL)
	{
		printf("Error of opening vcf aseq rec file\n"); 
		exit(1);
	}
	
	fp_in_rpos_end = fopen(vcf_rpos_end_rev,"rb");
	if(fp_in_rpos_end == NULL)
	{
		printf("Error of opening vcf rpos end rec file\n"); 
		exit(1);
	}
	
	fp_in_ref_seq_v = fopen(ref_vcf_rev,"rb");
	if(fp_in_ref_seq_v == NULL)
	{
		printf("Error of opening vcf ref seq rec file\n"); 
		exit(1);
	}
	
	fp_in_size = fopen(vcf_size_rev,"rb");
	if(fp_in_size == NULL)
	{
		printf("Error of opening vcf size rec file\n"); 
		exit(1);
	}
	
	fp_in_type = fopen(vcf_type_rev,"rb");
	if(fp_in_type == NULL)
	{
		printf("Error of opening vcf size rec file\n"); 
		exit(1);
	}
	
	fread(&rpos_cnt, 4, 1, fp_in_size);
	fread(&apos_cnt, 4, 1, fp_in_size);
	fread(&aseq_cnt, 4, 1, fp_in_size);
	fread(&ref_cv_n, 4, 1, fp_in_size);
	
	rpos_array_rec = (uint32_t* )calloc(rpos_cnt, 4);
	if(rpos_array_rec == NULL){
		printf("Error of allocating memory rpos_array_rec\n");
		exit(3);
	}
	aindex_array_rec = (uint32_t* )calloc(rpos_cnt + 1, 4); 
	if(aindex_array_rec == NULL){
		printf("Error of allocating memory aindex_array_rec\n");
		exit(3);
	}
	apos_array_rec = (uint32_t* )calloc(apos_cnt + 1, 4);
	if(apos_array_rec == NULL){
		printf("Error of allocating memory apos_array_rec\n");
		exit(3);
	}
	rpos_end_array_rec = (uint32_t* )calloc(apos_cnt + 1, 4);
	if(rpos_end_array_rec == NULL){
		printf("Error of allocating memory rpos_end_array_rec\n");
		exit(3);
	}
	aseq_array_rec = (uint8_t* )calloc(aseq_cnt, 1);
	if(aseq_array_rec == NULL){
		printf("Error of allocating memory aseq_array_rec\n");
		exit(3);
	}
	ref_seqv_array_rec = (uint64_t* )calloc((ref_cv_n >> 1) + 1000, 1);//(ref_cv_n >> 4 + 1) + 1000, 4
	if(ref_seqv_array_rec == NULL){
		printf("Error of allocating memory ref_seqv_array_rec\n");
		exit(3);
	}
	type_array_rec = (uint8_t* )calloc(apos_cnt + 1, 1);
	if(type_array_rec == NULL){
		printf("Error of allocating memory type_array_rec\n");
		exit(3);
	}
	
	if((rpos_rsize = fread(rpos_array_rec, 4, rpos_cnt, fp_in_rpos)) != rpos_cnt)
	{
		fprintf (stderr,"Reading error of rpos rev file %u %u",rpos_cnt,rpos_rsize);
		exit(2);
	}

	if((aindex_rsize = fread(aindex_array_rec, 4, rpos_cnt + 1, fp_in_aindex)) != (rpos_cnt + 1))
	{
		fputs ("Reading error of aindex rev file",stderr);
		exit(2);
	}

	if((apos_rsize = fread(apos_array_rec, 4, apos_cnt + 1, fp_in_apos)) != (apos_cnt + 1))
	{
		fputs ("Reading error of apos rev file",stderr);
		exit(2);
	}

	if((rpos_end_rsize = fread(rpos_end_array_rec, 4, apos_cnt, fp_in_rpos_end)) != apos_cnt)
	{
		fputs ("Reading error of rpos end rev file",stderr);
		exit(2);
	}

	if((ref_seqv_rsize = fread(ref_seqv_array_rec, 8, (ref_cv_n >> 4) + 1, fp_in_ref_seq_v)) != ((ref_cv_n >> 4) + 1))
	{

		fputs ("Reading error of ref seq vcf rev file",stderr);
		exit(2);
	}

	if((aseq_rsize = fread(aseq_array_rec, 1, aseq_cnt, fp_in_aseq)) != aseq_cnt)
	{
		fputs ("Reading error of aseq end rev file",stderr);
		exit(2);
	}

	if((type_rsize = fread(type_array_rec, 1, apos_cnt, fp_in_type)) != apos_cnt)
	{
		fputs ("Reading error of vcf type rev file",stderr);
		exit(2);
	}
	
	if(fp_in_rpos)	fclose(fp_in_rpos);
	if(fp_in_aindex)	fclose(fp_in_aindex);
	if(fp_in_apos)	fclose(fp_in_apos);
	if(fp_in_rpos_end)	fclose(fp_in_rpos_end);
	if(fp_in_aseq)	fclose(fp_in_aseq);
	if(fp_in_ref_seq_v)	fclose(fp_in_ref_seq_v);
	if(fp_in_size)	fclose(fp_in_size);
	if(fp_in_type)	fclose(fp_in_type);	
}

