
#ifndef ALI_EXT_H_
#define ALI_EXT_H_


#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

//#include "load_input.h"
#include "LandauVishkin_tree.h"

#define	MAX_K_VAR	16
#define	MAX_K_VAR_K	15

#define	EXT_ANCHOR

//tree extension
extern uint8_t** seqss1;
extern uint8_t** seqss2;
extern treenode***	node_arrays;
extern uint16_t*** seq_nodes;
extern uint16_t** 	seq_node_ns;
extern uint16_t***	seq_node_lens;
extern uint16_t** 	seq_node_len_ns;
extern alt_info_a** alt_info_arrays;
extern alt_seq_a** alt_seqss;
extern uint16_t** route_seps;
extern ext_queue*** lv_queues;
extern ext_result*** lv_results;
extern route_node*** rnode_array;
extern uint32_t*** routeid_arrays;
extern uint32_t** route_id_marks;
extern uint32_t** route_id_marks_pre;
extern ext_results** ex_results1;
extern ext_results** ex_results2;
extern com_result** c_results1;
extern com_result** c_results2;
extern uint8_t*** ref_seq_ori;
extern uint8_t*** ref_seq_var;
extern char** cigar_ps;
extern char** f_cigars;
extern char** b_cigars;
extern char** str_os;

extern ext_output*** ext_out1s;
extern ext_output*** ext_out2s;
extern uint16_t** mat_index1s;
extern uint16_t** mat_index2s;

extern const uint8_t max_extension_length;

extern uint16_t* end_dis1;
extern uint16_t* end_dis2;

extern int** op_dm_l1;
extern int** op_dm_r1;
extern int** ops_dm_l1;
extern int** ops_dm_r1;
extern int** op_dm_ex1;
extern int** op_dm_ex2;
extern int** ops_dm_ex1;
extern int** ops_dm_ex2;

extern int** chr_res;
extern uint32_t** sam_pos1s;
extern uint32_t** sam_pos2s;
extern char*** cigar_p1s;
extern char*** cigar_p2s;
extern char** xa_d1s;
extern char** xa_d2s;
extern int** lv_re1s;
extern int** lv_re2s;

extern seq_io* seqio;

extern char** qual1_buffer;
extern char** qual2_buffer;
extern char** read_rev_buffer;
extern char** read_rev_buffer_1;
extern char** pr_cigar1_buffer;
extern char** pr_cigar2_buffer;

int compare_ext_out(const void * , const void * );
int binsearch_pair_pos_ext(uint32_t , ext_output** , uint16_t , int );
	
#endif /* ALI_EXT_H_ */
