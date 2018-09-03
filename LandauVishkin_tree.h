
#ifndef LANDAUVISHKIN_TREE_H_
#define LANDAUVISHKIN_TREE_H_

// Computes the edit distance between two strings and returns a CIGAR string for the edits.
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

#include "bit_operation.h"

#define	LV_MP

#define bool int
#define _ASSERT assert
#define _uint64 uint64_t
#define CountTrailingZeroes(x, ans) do{ans = __builtin_ctzll(x);} while(0)
#define __min(a, b) (a<b?a:b)
#define MAX_K_EXT 100
//250
//720
//500
//410
//820
//65
//129
#define true 1
#define false 0

#define	BUCKET_SIZE	16
#define	EXTENSION_REF_LENGTH	2048
//1000
#define	NODES_QUEUE_MAX	5900000
//10000
//100

#define	MAX_SEQ_N	1000
#define	MAX_ALT_CNT	1000	
#define	MAX_CIGAR_LEN	100

#define	TEST_RAW	10
#define	TEST_COL	2048

#define	END_FREE
#define	BEGIN_ALLOC
#define	ALT_TR

//debug
//#define	DEBUG_DETAIL
//#define	DEBUG_LRESULT
//#define	DEBUG_ROUTE

//#define	DEBUG_TREEALTROUTE

//sv extension 0 14265 57061 TAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAAC

//#define	LVE_DETAIL
//#define	LVE_RESULT

#define	FIX_SV

#define	MAX_READLEN	2049
#define MAX_COV_A	(((MAX_READLEN - 1) >> 6) + 1)
#define MAX_OP_SCORE	MAX_READLEN
#define	MAX_EXT_ALT_POS	2000
#define	CUS_SEED_SET	2000

//#define	CUS_MAX_OUTPUT_ALI	150
//#define CUS_MAX_OUTPUT_ALI2 (CUS_MAX_OUTPUT_ALI << 1)
#define CUS_MAX_OUTPUT_ALI2 2000

#define ANCHOR_HASH_ALI
#define	REDUCE_ANCHOR

#ifdef	ANCHOR_HASH_ALI

typedef struct anchor_seeds {
	int read_left_off;
	int read_right_off;
	uint32_t ref_left_off;
	uint32_t ref_right_off;
	uint16_t seed_length;
} anchor_seed;

extern uint8_t k_anchor;
extern uint8_t anchor_seed_d;
extern uint8_t k_anchor_back;
extern uint8_t anchor_back_mask;
extern uint8_t k_anchor_re;
extern uint8_t k_anchor_b;

extern uint16_t*** anchor_hash;
extern uint8_t*** anchor_array;
extern uint16_t*** anchor_point;
extern uint16_t*** anchor_pos;

extern anchor_seed** anchor_seed_buffer;

#endif

extern uint8_t anchor_seed_length_thr;//11: 1993234; 12 1992946; 15 1990***; 20 1989***; 30 1988***

extern uint32_t insert_dis;
extern int64_t devi;

typedef struct read_hs {
	uint32_t des;
	uint16_t off_set;
} read_h;

extern read_h** rh;

typedef struct seed_PARRAY
{
	uint64_t cov[MAX_COV_A];

	uint32_t pos_start;
	uint32_t pos_n;
	uint32_t ref_pos_off;
	uint32_t ref_pos_off_r;

	int s_r_o_l;
	int s_r_o_r;

	uint16_t length;
	uint8_t ui;
}seed_pa;

typedef struct SEQIOIN
{
	char* read_seq1;
	char* read_seq2;

	char* qual1;
	char* qual2;

	char* name;
	char* name_other;
	
	uint16_t read_length1;
	uint16_t read_length2;

	uint16_t length_h1;
	uint16_t length_h2;
	uint16_t nm1;
	uint16_t nm2;

	uint16_t v_cnt;
	uint8_t flag1;
	uint8_t flag2;
	int chr_re;
	int chr_re1;
	int chr_re2;

	int64_t pos1;
	int64_t pos2;
	int qualc1;
	int qualc2;
	int64_t cross;
	char* seq1;
	char* seq2;
	char* cigar1;
	char* cigar2;
	uint16_t xa_n;
	uint16_t xa_n_p1;
	uint16_t xa_n_p2;	
	uint16_t xa_n1;
	uint16_t xa_n2;
	
#ifdef	FIX_SV	
	uint16_t xa_n_x1;
	uint16_t xa_n_x2;
	int* chr_res_s1;
	int* chr_res_s2;
#endif
	int* chr_res;
	char* xa_d1s;
	uint32_t* sam_pos1s;
	char* cigar_p1s[CUS_MAX_OUTPUT_ALI2];
	char* cigar_p2s[CUS_MAX_OUTPUT_ALI2];

	int* lv_re1s;
	char* xa_d2s;
	uint32_t* sam_pos2s;

	int* lv_re2s;

	uint8_t* chr_res1;
	char* xa_ds1;
	uint32_t* sam_poss1;
	uint8_t* lv_res1;

	uint8_t* chr_res2;
	char* xa_ds2;
	uint32_t* sam_poss2;
	uint8_t* lv_res2;

}seq_io;

typedef struct Tree_Node //11
{
	uint8_t re_flag;
	uint16_t seqs_start;
	uint16_t seqs_length;
	uint16_t seq_dis;
	uint16_t seq_start;
	uint16_t next_nodes_n;
	uint16_t seq_nums_n;
	uint16_t node_offset;
	uint16_t alt_nums;
	uint16_t alt_info_cnt;
	uint16_t alt_info_n;
	uint16_t next_nodes[BUCKET_SIZE];	//how many nodes are there at most?
	uint16_t seq_nums[MAX_SEQ_N];
} treenode;

typedef struct route_node_info{
	uint8_t alt_type;
	uint8_t alt_flag;
	uint8_t cigar;
	
	int16_t pdis;
	uint16_t nodeid;
	uint16_t offset;
	uint16_t alt_index;
	uint32_t route_id;
	uint16_t del_p;
	//uint16_t nodern;
	uint32_t pre_index;
	
	//int16_t tdis;
	//uint16_t noder_begin;
	//uint16_t noder_len;
}route_node;

typedef struct ext_lv_queue{
	
	uint8_t alt_flag;
	uint8_t alt_type;
	uint16_t nodeid;
	uint16_t pdis;
	//uint16_t tdis;
	uint16_t offset;
	uint16_t alt_info_cnt;
	uint16_t del_p;
	uint32_t route_id;
}ext_queue;

typedef struct ext_lv_result{
	uint8_t edit;
	uint8_t cigar;
	uint16_t nodeid;
	uint16_t offset;
	uint32_t rnode_index;
	uint32_t routeid;
}ext_result;

typedef struct ext_alt_info{
	uint16_t nodeid;
	uint16_t node_offset;
	uint16_t alt_length;
	uint16_t ending_node_id;
	uint16_t ending_node_offset;
	uint16_t ori_index;
	uint16_t alt_start;
	uint16_t alt_end;
	uint16_t seq_index;
	uint32_t alt_seq_begin;
	uint32_t vcf_index;//debug
	uint8_t alt_type;//debug
}alt_info_a;

typedef struct alt_seq_info{
	uint32_t ref_index;
	uint32_t aindex;
}alt_seq_a;

typedef struct ori_ext_results
{
	//uint16_t route_id;
	//uint16_t node_id;
	//uint16_t seqpos;
	//uint16_t seqpos_n;
	uint16_t edit;
	uint16_t seq_dis;
	uint32_t chrno;
	char cigars[MAX_CIGAR_LEN];
	uint32_t altseq[MAX_CIGAR_LEN << 1];
} ext_results;


typedef struct ori_com_result
{
	uint16_t seqid;
	uint16_t resultid;
} com_result;

typedef struct ext_outputs
{
	uint8_t dir;
	uint16_t edit;
	uint32_t chrno;
	uint32_t pos;//int64_t
	char cigar[MAX_CIGAR_LEN];
} ext_output;

uint16_t max_route_ids;

void load_vcf_chr(char* );
void chr_vcf_index();
void load_input_sorted_vcf();
int compare_vcf_sort(const void * , const void * );
int compare_alt_sort(const void * , const void * );
int compare_com_result(const void * , const void * );
int compare_read_hash(const void * , const void * );
int comepare_anchor_seed(const void * , const void * );
int binsearch_com(uint16_t , com_result* , int );

uint16_t ext_single_tree_alt(uint8_t** , uint8_t** , uint8_t* , uint8_t** , ext_results** , com_result** , uint16_t* , uint8_t , uint32_t* , uint16_t , uint16_t , uint16_t , uint8_t , uint16_t, uint64_t*, uint8_t );
void create_route_array(uint16_t , uint16_t , uint16_t** , uint32_t** , uint32_t* , uint32_t** , uint32_t* , alt_info_a* );
int compare_ext_out(const void * , const void * );
int binsearch_pair_pos_ext(uint32_t , ext_output** , uint16_t, int );
char* combine_cigars(char* , char* , int16_t , int16_t , uint16_t , uint8_t , uint8_t , uint8_t );


int computeEditDistance_tree(
    treenode** ,
    uint8_t* ,
    int ,
    int ,
    ext_queue** ,
    ext_result*** ,
    uint16_t* ,
    alt_info_a* ,
    uint8_t* ,
    uint8_t* ,
    uint32_t,//
	uint32_t*, 
	uint32_t*,
	route_node** ,
	uint32_t** ,
	uint16_t* ,
	uint32_t* ,
	uint32_t* ,
	uint16_t,
	uint32_t*,
	uint16_t*);

int compare_route_mask(const void * , const void * );

#endif /* LANDAUVISHKIN_TREE_H_ */