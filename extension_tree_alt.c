
//problem: unknown files: not enough array size to store vcf long alt record

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

#include <memory.h>
#include <math.h>
#include <inttypes.h>

#include "load_input.h"
#include "LandauVishkin_tree.h"
#include "extension_tree_alt.h"
//#include "binarys_qsort.h"

#define	EXT_TWICE

#define	MAX_EDIT	32

uint32_t arr_tol = 0;

//tree extension
uint8_t** seqss1 = NULL;
uint8_t** seqss2 = NULL;
treenode***	node_arrays = NULL;
uint16_t*** seq_nodes = NULL;
uint16_t** 	seq_node_ns = NULL;
uint16_t***	seq_node_lens = NULL;
uint16_t** 	seq_node_len_ns = NULL;
alt_info_a** alt_info_arrays = NULL;
alt_seq_a** alt_seqss = NULL;
uint16_t** route_seps = NULL;
ext_queue*** lv_queues = NULL;
ext_result*** lv_results = NULL;
route_node*** rnode_array = NULL;
uint32_t*** routeid_arrays = NULL;
uint32_t** route_id_marks = NULL;
uint32_t** route_id_marks_pre = NULL;
ext_results** ex_results1 = NULL;
ext_results** ex_results2 = NULL;
com_result** c_results1 = NULL;
com_result** c_results2 = NULL;
uint8_t*** ref_seq_ori = NULL;
uint8_t*** ref_seq_var = NULL;
char** cigar_ps = NULL;
char** f_cigars = NULL;
char** b_cigars = NULL;
char** str_os = NULL;

ext_output*** ext_out1s = NULL;
ext_output*** ext_out2s = NULL;
uint16_t** mat_index1s = NULL;
uint16_t** mat_index2s = NULL;



//no consideration of repetitive function



uint8_t extension_tree_alt(seed_pa* seed_pr[2][2],
                           uint32_t spa_i[2][2],
                           uint8_t unmatch[2][2],
                           uint64_t** seed_set_pos[2][2],
                           uint64_t* buffer_ref_seq,
                           uint64_t* buffer_ref_seq_var,
                           uint8_t** ref_seq,
                           uint8_t** ref_seq_var,
                           uint8_t* read_char,
                           uint64_t read_bit1[2][((MAX_READLEN - 1) >> 5) + 1],
                           uint64_t read_bit2[2][((MAX_READLEN - 1) >> 5) + 1],
                           uint16_t read_len1,
                           uint16_t read_len2,
                           uint8_t tid,
                           uint16_t lv_k1,
                           uint16_t lv_k2,
                           uint32_t seqi,
                           int devi,
                           read_h* rh,
                           uint16_t** anchor_hash,
                           uint8_t** anchor_array,
                           uint16_t** anchor_point,
                           uint16_t** anchor_pos
                          )
{

	uint8_t rc_i = 0, un_ii = 0, l_j = 0, r_j = 0;
	uint8_t sam_flag1 = 0, sam_flag2 = 0;
	uint8_t re_flag = 1;
	uint8_t ori_i = 0;
	uint8_t base_re = 0;
	uint8_t anchor_hash_i = 0;
	uint8_t dir_tmp = 0;

	int16_t s_r_o_l = 0;
	uint16_t pos_n = 0, s_r_o_r = 0, read_b_i = 0, psp_i = 0, tmp_i = 0, xa_i = 0;
	uint16_t c_result1_n = 0, c_result2_n = 0;
	uint16_t ext_out_n1 = 0, ext_out_n2 = 0, ext_out_n1_tmp = 0, ext_out_n2_tmp = 0,ext_out_n = 0;
	uint16_t mat_n = 0;
	uint16_t tmp_edit = 0;
	uint16_t dm = 0;
	uint16_t read_len = 0;
	uint16_t k_tmp = 0;
	uint16_t mat_tmp1 = 0, mat_tmp2 = 0;
	uint16_t edit1 = 0, edit2 = 0;
	uint16_t v_cnt = 0, vs_cnt = 0, v_cnt1 = 0, vs_cnt1 = 0, v_cnt2 = 0, vs_cnt2 = 0;
	uint16_t dm_op = MAX_OP_SCORE, dm_ops = MAX_OP_SCORE;//, dm_op1 = 0, dm_ops1 = 0, dm_op2 = 0, dm_ops2 = 0;
	uint16_t v_cnt_i = 0;
	uint16_t tra_i;
	uint16_t des_i = 0;
	uint16_t print_i;
	uint16_t array_index;
	uint16_t char_i = 0;
	uint16_t buffer_i = 0;
	uint16_t seed_length = 0;
	uint16_t max_seed_length = 0;
	uint16_t lv_k = 0;
	uint16_t lv_k_1 = 0;
	uint16_t lv_k_2 = 0;
	uint16_t alt_num_tmp = 0;

	uint32_t read_k_p = 0;
	uint32_t read_k_t = 0;
	uint32_t anchor_mask = 0;
	uint32_t array_i = 0;
	uint32_t array_i_p = 0;
	uint32_t anchor_ref = 0;
	uint32_t r_b_v = 0;
	uint32_t max_right = 0;
	uint32_t pos_tmp = 0;
	uint32_t pos_s = 0;
	uint32_t seed1_i = 0;
	uint32_t posi = 0;
	uint32_t spa = 0;
	uint32_t chr_re1 = 0, chr_re2 = 0;
	uint32_t index1 = 0;
	uint32_t index2 = 0;

	int left_i = 0;
	int right_i = 0;
	int tmp_bin = 0;
	int bit_char_i = 0;
	int end_dis = 0;
	int end_dis_a[2];

	int64_t sam_cross = 0, sam_pos1 = 0, sam_pos2 = 0;

	uint64_t ksw_s = 0;
	uint64_t ksw_e = 0;
	uint64_t base_i = 0;
	uint64_t base_i_off_l = 0;
	uint64_t base_i_off_r = 0;

	char* cigar_tmp = NULL;
	char* cigar_tmp_o = NULL;
	char* tmp_cigar_l = NULL;
	char* tmp_cigar_r = NULL;
	char* output_cigar[2][2];
	char sam_seq[MAX_READLEN] = {};

	uint8_t* anchor_array_p = NULL;

	uint16_t* mat_index1 = NULL;
	uint16_t* mat_index2 = NULL;
	uint16_t* anchor_hash_p = NULL;
	uint16_t* anchor_point_p = NULL;
	uint16_t* anchor_pos_p = NULL;

	uint32_t posesp[MAX_EXT_ALT_POS];
	uint32_t posesp_l[MAX_EXT_ALT_POS];
	uint32_t posesp_r[MAX_EXT_ALT_POS];

	int low = 0, high = 0, mid = 0, chr_re = 0;
	
	int* op_rc_p1 = NULL;
	int* ops_rc_p1 = NULL;
	int* op_rc_p2 = NULL;
	int* ops_rc_p2 = NULL;
	int* op_dm_ex1_p = NULL;
	int* ops_dm_ex1_p = NULL;
	int* op_dm_ex2_p = NULL;
	int* ops_dm_ex2_p = NULL;

	uint64_t* read_bit = NULL;
	uint64_t* read_bit_2 = NULL;
	uint64_t* ori = NULL;

	seed_pa* seed_pr1 = NULL;

	ext_output* ext_p1 = NULL;
	ext_output* ext_p2 = NULL;
	ext_output** ext_out1 = NULL;
	ext_output** ext_out2 = NULL;
	ext_output** ext_out = NULL;

	com_result* c_result1 = NULL;
	com_result* c_result2 = NULL;
	ext_results* ex_result1 = NULL;
	ext_results* ex_result2 = NULL;


	ext_out1 = ext_out1s[tid];
	ext_out2 = ext_out2s[tid];
	mat_index1 = mat_index1s[tid];
	mat_index2 = mat_index2s[tid];

	max_route_ids = (MAX_K_VAR << 1) + 2;

	anchor_mask = ((uint32_t )1 << k_anchor_b) - 1;

	for(rc_i = 0; rc_i < 2; rc_i++)
	{
		for(un_ii = 0; un_ii < 2; un_ii++)
		{
			if(unmatch[rc_i][un_ii] == 1)	continue;

			seed_pr1 = seed_pr[rc_i][un_ii];
			spa = spa_i[rc_i][un_ii];

			if(spa)
			{
				if(rc_i ^ un_ii)  	//end2
				{
					read_bit = read_bit2[un_ii];
					read_len = read_len2;
					//k_tmp = lv_k2;
					ext_out = ext_out2;
					ext_out_n = ext_out_n2;
				}
				else  			//end1
				{
					read_bit = read_bit1[un_ii];
					read_len = read_len1;
					//k_tmp = lv_k1;
					ext_out = ext_out1;
					ext_out_n = ext_out_n1;
				}
			}

			for(seed1_i = 0; seed1_i < spa; seed1_i++)
			{
				s_r_o_l = seed_pr1[seed1_i].s_r_o_l;
				s_r_o_r = seed_pr1[seed1_i].s_r_o_r;

				for(psp_i = 0, pos_n = 0; psp_i < seed_pr1[seed1_i].pos_n; psp_i++)
				{
					posi = seed_set_pos[rc_i][un_ii][tid][seed_pr1[seed1_i].pos_start + psp_i] - 1;// - max_extension_length

					if(pos_n < MAX_EXT_ALT_POS)
					{
						posesp[pos_n] = posi;
						pos_n++;

					}
				}

				if((s_r_o_l == -1) && (s_r_o_r == read_len))
				{
					for(psp_i = 0; psp_i < pos_n; psp_i++)
					{
						pos_tmp = posesp[psp_i];

						low = 0;
						high = chr_file_n - 1;
						while ( low <= high )
						{
							mid = (low + high) >> 1;
							if(pos_tmp < (chr_end_n[mid]))
							{
								high = mid - 1;
							}
							else if(pos_tmp > (chr_end_n[mid]))
							{
								low = mid + 1;
							}
							else
							{
								chr_re =  mid;
								break;
							}
							chr_re = low;
						}

						/**/
						if(ext_out_n < TREE_EXT_OUTPUT)
						{
							ext_out[ext_out_n]->chrno = chr_re;
							ext_out[ext_out_n]->edit = 0;
							ext_out[ext_out_n]->pos = pos_tmp;
							ext_out[ext_out_n]->dir = un_ii;

							sprintf(ext_out[ext_out_n]->cigar, "%uM", read_len);

							ext_out_n++;
						}

					}

					if(rc_i ^ un_ii)	ext_out_n2 = ext_out_n;
					else	ext_out_n1 = ext_out_n;

					continue;
				}

				c_result1_n = 0;
				c_result2_n = 0;
				//ext_out_n = 0;

				//left alignment
				if(s_r_o_l != -1)
				{
					for(bit_char_i = s_r_o_l, read_b_i = 0; bit_char_i >= 0; bit_char_i--, read_b_i++)
						read_char[read_b_i] = ((read_bit[bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

					for(psp_i = 0; psp_i < pos_n; psp_i++)
					{
						for(bit_char_i = s_r_o_l + MAX_EDIT, read_b_i = 0; bit_char_i > - 1; bit_char_i--, read_b_i++)
						{
							pos_tmp = posesp[psp_i] + bit_char_i - MAX_EDIT;
							ref_seq[psp_i][read_b_i] = ((buffer_ref_seq[pos_tmp >> 5] >> ((31 - (pos_tmp & 0X1f)) << 1)) & 0X3);
							ref_seq_var[psp_i][read_b_i] = ((buffer_ref_seq_var[pos_tmp >> 4] >> ((15 - (pos_tmp & 0Xf)) << 2)) & 0Xf);

						}
						posesp_l[psp_i] = ref_length_p - (posesp[psp_i] + s_r_o_l);

					}
					k_tmp = ((s_r_o_l + 1) / 10);
					k_tmp = (k_tmp < 2?2:k_tmp);
					k_tmp = (k_tmp > 5?5:k_tmp);

					alt_num_tmp = ext_single_tree_alt(ref_seq, ref_seq_var, read_char, seqss1, ex_results1, c_results1, &c_result1_n, tid, posesp_l, pos_n, s_r_o_l + 1 + MAX_EDIT, s_r_o_l + 1, 0, k_tmp, buffer_ref_seq_var, 1);

					if(!c_result1_n)
					{
						//seed1_i++;
						continue;
					}
				}


				//right alignment
				if(s_r_o_r < read_len)
				{
					for(bit_char_i = s_r_o_r, read_b_i = 0; bit_char_i < read_len; bit_char_i++, read_b_i++)
						read_char[read_b_i] = ((read_bit[bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

					for(psp_i = 0; psp_i < pos_n; psp_i++)
					{
						for(bit_char_i = s_r_o_r, read_b_i = 0; bit_char_i < read_len + MAX_EDIT; bit_char_i++, read_b_i++)
						{
							pos_tmp = posesp[psp_i] + bit_char_i;
							ref_seq[psp_i][read_b_i] = ((buffer_ref_seq[pos_tmp >> 5] >> ((31 - (pos_tmp & 0X1f)) << 1)) & 0X3);
							ref_seq_var[psp_i][read_b_i] = ((buffer_ref_seq_var[pos_tmp >> 4] >> ((15 - (pos_tmp & 0Xf)) << 2)) & 0Xf);

						}
						posesp_r[psp_i] =  posesp[psp_i] + s_r_o_r;

					}
					k_tmp = ((read_len - s_r_o_r) / 10);
					k_tmp = (k_tmp < 2?2:k_tmp);
					k_tmp = (k_tmp > 5?5:k_tmp);

					alt_num_tmp = ext_single_tree_alt(ref_seq, ref_seq_var, read_char, seqss2, ex_results2, c_results2, &c_result2_n, tid, posesp_r, pos_n, read_len - s_r_o_r + MAX_EDIT, read_len - s_r_o_r, 1, k_tmp, buffer_ref_seq_var, 1);

					if(!c_result2_n)
					{
						continue;
					}

				}

				c_result1 = c_results1[tid];
				c_result2 = c_results2[tid];

				ex_result1 = ex_results1[tid];
				ex_result2 = ex_results2[tid];

				if(c_result1_n && c_result2_n)
				{
					l_j = 1;
					r_j = 1;
					if(c_result1_n < c_result2_n)
					{
						qsort(c_result2, c_result2_n, sizeof(com_result), compare_com_result);
						for(tmp_i = 0; tmp_i < c_result1_n; tmp_i++)
						{
							if((tmp_bin = binsearch_com(c_result1[tmp_i].seqid, c_result2, c_result2_n)) != -1)   //found
							{
								tmp_edit = ex_result1[c_result1[tmp_i].resultid].edit + ex_result2[c_result2[tmp_bin].resultid].edit;

								//combine cigars
								cigar_tmp = combine_cigars(ex_result1[c_result1[tmp_i].resultid].cigars, ex_result2[c_result2[tmp_bin].resultid].cigars, s_r_o_l, s_r_o_r, read_len, tid, l_j, r_j);

								if(ext_out_n < TREE_EXT_OUTPUT)
								{
									ext_out[ext_out_n]->edit = tmp_edit;
									ext_out[ext_out_n]->chrno = ex_result2[c_result2[tmp_bin].resultid].chrno;
									ext_out[ext_out_n]->pos = posesp[c_result1[tmp_i].seqid] + s_r_o_l - ex_result1[c_result1[tmp_i].resultid].seq_dis;
									ext_out[ext_out_n]->dir = un_ii;
									strcpy(ext_out[ext_out_n]->cigar, cigar_tmp);

									ext_out_n++;
								}
							}
						}
					}
					else
					{
						qsort(c_result1, c_result1_n, sizeof(com_result), compare_com_result);
						for(tmp_i = 0; tmp_i < c_result2_n; tmp_i++)
						{
							if((tmp_bin = binsearch_com(c_result2[tmp_i].seqid, c_result1, c_result1_n)) != -1)   //found
							{
								tmp_edit = ex_result2[c_result2[tmp_i].resultid].edit + ex_result1[c_result1[tmp_bin].resultid].edit;

								//combine cigars
								cigar_tmp = combine_cigars(ex_result1[c_result1[tmp_bin].resultid].cigars, ex_result2[c_result2[tmp_i].resultid].cigars, s_r_o_l, s_r_o_r, read_len, tid, l_j, r_j);

								if(ext_out_n < TREE_EXT_OUTPUT)
								{
									ext_out[ext_out_n]->edit = tmp_edit;
									ext_out[ext_out_n]->chrno = ex_result2[c_result2[tmp_i].resultid].chrno;
									ext_out[ext_out_n]->pos = posesp[c_result1[tmp_bin].seqid] + s_r_o_l - ex_result1[c_result1[tmp_bin].resultid].seq_dis;
									ext_out[ext_out_n]->dir = un_ii;
									strcpy(ext_out[ext_out_n]->cigar, cigar_tmp);

									ext_out_n++;
								}
							}
						}
					}
				}
				else if(c_result1_n)
				{
					l_j = 1;
					r_j = 0;
					for(tmp_i = 0; tmp_i < c_result1_n; tmp_i++)
					{
						tmp_edit = ex_result1[c_result1[tmp_i].resultid].edit;

						//combine cigars
						cigar_tmp = combine_cigars(ex_result1[c_result1[tmp_i].resultid].cigars, NULL, s_r_o_l, s_r_o_r, read_len, tid, l_j, r_j);

						if(ext_out_n < TREE_EXT_OUTPUT)
						{
							ext_out[ext_out_n]->edit = tmp_edit;
							ext_out[ext_out_n]->chrno = ex_result1[c_result1[tmp_i].resultid].chrno;
							ext_out[ext_out_n]->pos = posesp[c_result1[tmp_i].seqid] + s_r_o_l - ex_result1[c_result1[tmp_i].resultid].seq_dis;
							ext_out[ext_out_n]->dir = un_ii;
							strcpy(ext_out[ext_out_n]->cigar, cigar_tmp);

							ext_out_n++;
						}
					}
				}
				else if(c_result2_n)
				{
					l_j = 0;
					r_j = 1;
					for(tmp_i = 0; tmp_i < c_result2_n; tmp_i++)
					{
						tmp_edit = ex_result2[c_result2[tmp_i].resultid].edit;

						//combine cigars
						cigar_tmp = combine_cigars(NULL, ex_result2[c_result2[tmp_i].resultid].cigars, s_r_o_l, s_r_o_r, read_len, tid, l_j, r_j);

						if(ext_out_n < TREE_EXT_OUTPUT)
						{
							ext_out[ext_out_n]->edit = tmp_edit;
							ext_out[ext_out_n]->chrno = ex_result2[c_result2[tmp_i].resultid].chrno;
							ext_out[ext_out_n]->pos = posesp[c_result2[tmp_i].seqid];
							ext_out[ext_out_n]->dir = un_ii;
							strcpy(ext_out[ext_out_n]->cigar, cigar_tmp);

							ext_out_n++;
						}
					}
				}

				if(rc_i ^ un_ii)	ext_out_n2 = ext_out_n;
				else	ext_out_n1 = ext_out_n;
			}
		}
	}

	if(ext_out_n1 < ext_out_n2)
	{
		qsort(ext_out2, ext_out_n2, sizeof(ext_output* ), compare_ext_out);

		end_dis_a[0] = end_dis1[tid];
		end_dis_a[1] = -end_dis2[tid];

		for(tmp_i = 0; tmp_i < ext_out_n1; tmp_i++)
		{
			end_dis = end_dis_a[ext_out1[tmp_i]->dir];

			posi = ext_out1[tmp_i]->pos + end_dis;

			if((tmp_bin = binsearch_pair_pos_ext(posi, ext_out2, ext_out_n2, devi)) != -1)
			{
				if(mat_n < CUS_SEED_SET)
				{
					mat_index1[mat_n] = tmp_i;
					mat_index2[mat_n] = tmp_bin;

					mat_n++;
				}

				for(bit_char_i = tmp_bin + 1; (ext_out2[bit_char_i]->pos < posi + devi) && (bit_char_i < ext_out_n2); bit_char_i++)
				{
					if(mat_n < CUS_SEED_SET)
					{
						mat_index1[mat_n] = tmp_i;
						mat_index2[mat_n] = bit_char_i;

						mat_n++;
					}

				}

				for(bit_char_i = tmp_bin - 1; (bit_char_i > -1) && (ext_out2[bit_char_i]->pos > posi - devi); bit_char_i--)
				{
					if(mat_n < CUS_SEED_SET)
					{
						mat_index1[mat_n] = tmp_i;
						mat_index2[mat_n] = bit_char_i;

						mat_n++;
					}
				}
			}
		}
	}
	else
	{
		qsort(ext_out1, ext_out_n1, sizeof(ext_output* ), compare_ext_out);

		end_dis_a[0] = end_dis2[tid];
		end_dis_a[1] = -end_dis1[tid];

		for(tmp_i = 0; tmp_i < ext_out_n2; tmp_i++)
		{
			end_dis = end_dis_a[ext_out1[tmp_i]->dir];

			posi = ext_out2[tmp_i]->pos + end_dis;
			if((tmp_bin = binsearch_pair_pos_ext(posi, ext_out1, ext_out_n1, devi)) != -1)
			{
				if(mat_n < CUS_SEED_SET)
				{
					mat_index1[mat_n] = tmp_bin;
					mat_index2[mat_n] = tmp_i;

					mat_n++;
				}

#ifdef	PAIR_TRAVERSE
				for(bit_char_i = tmp_bin + 1; (ext_out1[bit_char_i]->pos < posi + devi) && (bit_char_i < ext_out_n1); bit_char_i++)
				{
					if(mat_n < CUS_SEED_SET)
					{
						mat_index1[mat_n] = bit_char_i;
						mat_index2[mat_n] = tmp_i;

						mat_n++;
					}

				}

				for(bit_char_i = tmp_bin - 1; (ext_out1[bit_char_i]->pos > posi - devi) && (bit_char_i > -1); bit_char_i--)
				{
					if(mat_n < CUS_SEED_SET)
					{
						mat_index1[mat_n] = bit_char_i;
						mat_index2[mat_n] = tmp_i;

						mat_n++;
					}
				}
#endif
			}
		}

	}

	op_rc_p1 = op_dm_l1[tid];
	op_rc_p2 = op_dm_r1[tid];

	ops_rc_p1 = ops_dm_l1[tid];
	ops_rc_p2 = ops_dm_r1[tid];

	op_dm_ex1_p = op_dm_ex1[tid];
	op_dm_ex2_p = op_dm_ex2[tid];

	ops_dm_ex1_p = ops_dm_ex1[tid];
	ops_dm_ex2_p = ops_dm_ex2[tid];

	ext_output* ext_out_p = NULL;

	uint8_t find_anchor_flag = 0;
	
	alt_num_tmp = 0;
	
	if(mat_n == 0)
	{
		for(anchor_hash_i = 0; anchor_hash_i < 4; anchor_hash_i++)
		{
			if(anchor_hash_i == 0)
			{
				read_len = read_len2;
				ori = read_bit2[1];
			}
			else if(anchor_hash_i == 1)
			{
				read_len = read_len1;
				ori = read_bit1[0];
			}
			else if(anchor_hash_i == 2)
			{
				read_len = read_len1;
				ori = read_bit1[1];
			}
			else
			{
				read_len = read_len2;
				ori = read_bit2[0];
			}

			des_i = 0;
			for(ori_i = 0; ori_i < read_len >> 5; ori_i++)
			{
				for(char_i = 0; (char_i <= 64 - k_anchor_b) && (des_i <= read_len - k_anchor); char_i += 2, des_i++)
				{
					rh[des_i].des = ((ori[ori_i] >> (64 - char_i - k_anchor_b)) & anchor_mask);
					rh[des_i].off_set = des_i;
				}
				//note that boundary value
				for(char_i = 2; (char_i < k_anchor_b) && (des_i <= read_len - k_anchor); char_i += 2, des_i++)
				{
					rh[des_i].des = ((ori[ori_i] & anchor_mask_boundary_re[char_i >> 1]) << char_i) | (ori[ori_i + 1] >> (64 - char_i));// & anchor_mask_boundary[char_i >> 1]
					rh[des_i].off_set = des_i;
				}
			}

			qsort(rh, des_i, sizeof(read_h), compare_read_hash);

			memset(anchor_hash[anchor_hash_i], 0, 1 << 17);

			read_k_p = rh[0].des;
			anchor_hash[anchor_hash_i][read_k_p >> k_anchor_back] = 0;
			anchor_array[anchor_hash_i][0] = (read_k_p & anchor_back_mask);
			anchor_point[anchor_hash_i][0] = 0;
			anchor_pos[anchor_hash_i][0] = rh[0].off_set;
			array_i = 1;
			array_i_p = 0;

			for(char_i = 1; char_i < des_i; char_i++)
			{
				anchor_pos[anchor_hash_i][char_i] = rh[char_i].off_set;

				read_k_t = rh[char_i].des;

				if((read_k_t >> k_anchor_back) == (read_k_p >> k_anchor_back))
				{
					if((read_k_t & anchor_back_mask) != (read_k_p & anchor_back_mask))
					{
						anchor_array[anchor_hash_i][array_i] = (read_k_t & anchor_back_mask);
						anchor_point[anchor_hash_i][array_i] = char_i;
						array_i++;
					}
				}
				else
				{
					anchor_hash[anchor_hash_i][read_k_t >> k_anchor_back] = (array_i << 5);
					anchor_hash[anchor_hash_i][read_k_p >> k_anchor_back] |= (array_i - array_i_p);

					array_i_p = array_i;

					anchor_array[anchor_hash_i][array_i] = (read_k_t & anchor_back_mask);
					anchor_point[anchor_hash_i][array_i] = char_i;
					array_i++;
				}
				read_k_p = read_k_t;
			}
			anchor_hash[anchor_hash_i][read_k_p >> k_anchor_back] |= (array_i - array_i_p);
			anchor_point[anchor_hash_i][array_i] = char_i;
		}

		ext_out_n1_tmp = ext_out_n1;
		ext_out_n2_tmp = ext_out_n2;

		//
		find_anchor_flag = 0;
		for(tmp_i = 0; tmp_i < ext_out_n1 + ext_out_n2; tmp_i++)
		{
			if(tmp_i < ext_out_n1)
			{
				ext_out_p = ext_out1[tmp_i];

				anchor_hash_i = ext_out_p->dir;
				dir_tmp = (!anchor_hash_i);
				anchor_hash_i = (anchor_hash_i << 1) + anchor_hash_i;

				ext_out = ext_out2;
				ext_out_n = ext_out_n2;
				read_len = read_len2;
			}
			else
			{
				ext_out_p = ext_out2[tmp_i - ext_out_n1];

				anchor_hash_i = ext_out_p->dir;
				dir_tmp = (!anchor_hash_i);
				anchor_hash_i = 2 - anchor_hash_i;

				ext_out = ext_out1;
				ext_out_n = ext_out_n1;
				read_len = read_len1;
			}

			if(!((ext_out_n < TREE_EXT_OUTPUT) && (mat_n < CUS_SEED_SET)))
				continue;

			chr_re = ext_out_p->chrno;
			posi = ext_out_p->pos;

			if(anchor_hash_i == 0)
			{
				ksw_s = posi + end_dis1[tid] - devi;
				ksw_e = posi + insert_dis + devi;
				read_bit_2 = read_bit2[1];
			}
			else if(anchor_hash_i == 1)
			{
				ksw_s = posi - (end_dis1[tid] + devi);
				ksw_e = posi - (end_dis1[tid] - devi) + read_len1;
				read_bit_2 = read_bit1[0];
			}
			else if(anchor_hash_i == 2)
			{
				ksw_s = posi + (end_dis2[tid] - devi);
				ksw_e = posi + insert_dis + devi;
				read_bit_2 = read_bit1[1];
			}
			else
			{
				ksw_s = posi - (end_dis2[tid] + devi);
				ksw_e = posi - (end_dis2[tid] - devi) + read_len2;
				read_bit_2 = read_bit2[0];
			}

			anchor_hash_p = anchor_hash[anchor_hash_i];
			anchor_array_p = anchor_array[anchor_hash_i];
			anchor_point_p = anchor_point[anchor_hash_i];
			anchor_pos_p = anchor_pos[anchor_hash_i];

			buffer_i = 0;
			r_b_v = 0;
			for(base_i = ksw_s - 1; base_i < ksw_e - k_anchor; base_i += anchor_seed_d)
			{
				if(base_i + k_anchor - 1 < r_b_v)	continue;

				base_re = (base_i & 0X1f);
				if(base_re <= k_anchor_re)
				{
					anchor_ref = ((buffer_ref_seq[base_i >> 5] >> ((k_anchor_re - base_re) << 1)) & anchor_mask);
				}
				else
				{
					anchor_ref = (((buffer_ref_seq[base_i >> 5] & anchor_mask_boundary[32 - base_re]) << ((base_re - k_anchor_re) << 1)) | (buffer_ref_seq[(base_i >> 5) + 1] >> ((32 + k_anchor_re - base_re) << 1)));
				}

				max_right = 0;
				for(tra_i = 0; tra_i < (anchor_hash_p[anchor_ref >> 4] & 0X1f); tra_i++)
				{
					array_index = (anchor_hash_p[anchor_ref >> 4] >> 5) + tra_i;
					if(anchor_array_p[array_index] == (anchor_ref & 0Xf))
					{
						max_seed_length = 0;
						for(print_i = anchor_point_p[array_index]; print_i < anchor_point_p[array_index + 1]; print_i++)
						{
							//extension on both sides
							for(left_i = anchor_pos_p[print_i] - 1, base_i_off_l = base_i - 1; (left_i >= 0) && (base_i_off_l >= ksw_s - 1); left_i--, base_i_off_l--)
							{
								if(((read_bit_2[left_i >> 5] >> ((31 - (left_i  & 0X1f)) << 1)) & 0X3) != ((buffer_ref_seq[base_i_off_l >> 5] >> ((31 - (base_i_off_l & 0X1f)) << 1)) & 0X3))
									break;
							}

							for(right_i = anchor_pos_p[print_i] + k_anchor, base_i_off_r = base_i + k_anchor; (right_i < read_len) && (base_i_off_r < ksw_e); right_i++, base_i_off_r++)
							{
								if(((read_bit_2[right_i >> 5] >> ((31 - (right_i  & 0X1f)) << 1)) & 0X3) != ((buffer_ref_seq[base_i_off_r >> 5] >> ((31 - (base_i_off_r & 0X1f)) << 1)) & 0X3))
									break;
							}

							seed_length = right_i - left_i - 1;
							if(seed_length > max_seed_length)
							{
								max_seed_length = seed_length;
								max_right = base_i_off_r;
							}

							anchor_seed_buffer[tid][buffer_i].read_left_off = left_i;
							anchor_seed_buffer[tid][buffer_i].read_right_off = right_i;
							anchor_seed_buffer[tid][buffer_i].ref_left_off = base_i_off_l;
							anchor_seed_buffer[tid][buffer_i].ref_right_off = base_i_off_r;
							anchor_seed_buffer[tid][buffer_i].seed_length = seed_length;
							buffer_i++;
						}
						break;
					}
				}
				r_b_v = max_right;
			}

			//if(buffer_i == 0)	tol_s_n++;
			
			if((buffer_i > 0) && (max_seed_length > anchor_seed_length_thr))
			{
				find_anchor_flag = 1;
				
				qsort(anchor_seed_buffer[tid], buffer_i, sizeof(anchor_seed), comepare_anchor_seed);

				s_r_o_l = anchor_seed_buffer[tid][0].read_left_off;
				s_r_o_r = anchor_seed_buffer[tid][0].read_right_off;
				base_i_off_l = anchor_seed_buffer[tid][0].ref_left_off;
				base_i_off_r = anchor_seed_buffer[tid][0].ref_right_off;

				pos_s = base_i_off_l - s_r_o_l;//base_i_off_r - s_r_o_r

				if((s_r_o_l == -1) && (s_r_o_r == read_len))
				{
					sprintf(sam_seq, "%uM\0", read_len);
					cigar_tmp = sam_seq;
					pos_tmp = pos_s;
					tmp_edit = 0;
				}
				else
				{

					c_result1_n = 0;
					c_result2_n = 0;
					l_j = 0;
					r_j = 0;
					//left alignment
					if(s_r_o_l != -1)
					{
						for(bit_char_i = s_r_o_l, read_b_i = 0; bit_char_i >= 0; bit_char_i--, read_b_i++)
							read_char[read_b_i] = ((read_bit_2[bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);

						l_j = 1;
						psp_i = 0;
						pos_n = 1;
						// + 50
						for(bit_char_i = s_r_o_l + MAX_EDIT, read_b_i = 0; bit_char_i > - 1; bit_char_i--, read_b_i++)
						{
							pos_tmp = pos_s + bit_char_i - MAX_EDIT;
							ref_seq[psp_i][read_b_i] = ((buffer_ref_seq[pos_tmp >> 5] >> ((31 - (pos_tmp & 0X1f)) << 1)) & 0X3);
							ref_seq_var[psp_i][read_b_i] = ((buffer_ref_seq_var[pos_tmp >> 4] >> ((15 - (pos_tmp & 0Xf)) << 2)) & 0Xf);

						}
						posesp_l[psp_i] = ref_length_p - (pos_s + s_r_o_l);//posesp[psp_i]

#ifdef	EXT_TWICE
						k_tmp = ((s_r_o_l + 1) / 7);
						k_tmp = (k_tmp < 2?2:k_tmp);
						k_tmp = (k_tmp > 7?7:k_tmp);
#else
						k_tmp = ((s_r_o_l + 1) / 7);
						k_tmp = (k_tmp < 2?2:k_tmp);
						k_tmp = (k_tmp > 8?8:k_tmp);
#endif

						alt_num_tmp = ext_single_tree_alt(ref_seq, ref_seq_var, read_char, seqss1, ex_results1, c_results1, &c_result1_n, tid, posesp_l, pos_n, s_r_o_l + 1 + MAX_EDIT, s_r_o_l + 1, 0, k_tmp, buffer_ref_seq_var, 1);					
			
#ifdef	EXT_TWICE
						if((!alt_num_tmp) && (!c_result1_n))
						{
							k_tmp = ((s_r_o_l + 1) / 7) + 7;
							k_tmp = (k_tmp > MAX_K_VAR_K ? MAX_K_VAR_K : k_tmp);
							k_tmp = (k_tmp > (s_r_o_l + 1)?s_r_o_l:k_tmp);
				
							ext_single_tree_alt(ref_seq, ref_seq_var, read_char, seqss1, ex_results1, c_results1, &c_result1_n, tid, posesp_l, pos_n, s_r_o_l + 1 + MAX_EDIT, s_r_o_l + 1, 0, k_tmp, buffer_ref_seq_var, 0);
						}
#endif				

						if(!c_result1_n)
						{
							continue;
						}
					}


					//right alignment
					if(s_r_o_r < read_len)
					{
						for(bit_char_i = s_r_o_r, read_b_i = 0; bit_char_i < read_len; bit_char_i++, read_b_i++)
							read_char[read_b_i] = ((read_bit_2[bit_char_i >> 5] >> ((31 - (bit_char_i & 0X1f)) << 1)) & 0X3);
						
						r_j = 1;
						psp_i = 0;
						pos_n = 1;

						for(bit_char_i = s_r_o_r, read_b_i = 0; bit_char_i < read_len + MAX_EDIT; bit_char_i++, read_b_i++)
						{
							pos_tmp = pos_s + bit_char_i;
							ref_seq[psp_i][read_b_i] = ((buffer_ref_seq[pos_tmp >> 5] >> ((31 - (pos_tmp & 0X1f)) << 1)) & 0X3);
							ref_seq_var[psp_i][read_b_i] = ((buffer_ref_seq_var[pos_tmp >> 4] >> ((15 - (pos_tmp & 0Xf)) << 2)) & 0Xf);

						}
						posesp_r[psp_i] = pos_s + s_r_o_r;//posesp[psp_i]

#ifdef	EXT_TWICE
						k_tmp = ((read_len - s_r_o_r) / 7);
						k_tmp = (k_tmp < 2?2:k_tmp);
						k_tmp = (k_tmp > 7?7:k_tmp);
#else
						k_tmp = ((read_len - s_r_o_r) / 7);//5
						k_tmp = (k_tmp < 2?2:k_tmp);
						k_tmp = (k_tmp > 8?8:k_tmp);//9
#endif							

						alt_num_tmp = ext_single_tree_alt(ref_seq, ref_seq_var, read_char, seqss2, ex_results2, c_results2, &c_result2_n, tid, posesp_r, pos_n, read_len - s_r_o_r + MAX_EDIT, read_len - s_r_o_r, 1, k_tmp, buffer_ref_seq_var, 1);
				
#ifdef	EXT_TWICE
						if((!alt_num_tmp) && (!c_result2_n))
						{
							k_tmp = ((read_len - s_r_o_r) / 7) + 7;
							k_tmp = (k_tmp > MAX_K_VAR_K ? MAX_K_VAR_K : k_tmp);
							k_tmp = (k_tmp > (read_len - s_r_o_r)?(read_len - s_r_o_r - 1):k_tmp);
						
							ext_single_tree_alt(ref_seq, ref_seq_var, read_char, seqss2, ex_results2, c_results2, &c_result2_n, tid, posesp_r, pos_n, read_len - s_r_o_r + MAX_EDIT, read_len - s_r_o_r, 1, k_tmp, buffer_ref_seq_var, 0);

						}
#endif

						if(!c_result2_n)
						{
							continue;
						}

					}

					ex_result1 = ex_results1[tid];
					ex_result2 = ex_results2[tid];

					if(s_r_o_l == -1)
					{
						tmp_edit = ex_result2[0].edit;
						tmp_cigar_l = NULL;
						tmp_cigar_r = ex_result2[0].cigars;
						pos_tmp = pos_s;
					}
					else if(s_r_o_r == read_len)
					{
						tmp_edit = ex_result1[0].edit;
						tmp_cigar_l = ex_result1[0].cigars;
						tmp_cigar_r = NULL;
						pos_tmp = pos_s + s_r_o_l - ex_result1[0].seq_dis;

					}
					else
					{
						tmp_edit = ex_result1[0].edit + ex_result2[0].edit;
						tmp_cigar_l = ex_result1[0].cigars;
						tmp_cigar_r = ex_result2[0].cigars;
						pos_tmp = pos_s + s_r_o_l - ex_result1[0].seq_dis;
					}

					//combine cigars
					cigar_tmp = combine_cigars(tmp_cigar_l, tmp_cigar_r, s_r_o_l, s_r_o_r, read_len, tid, l_j, r_j);
				}

				ext_out[ext_out_n]->edit = tmp_edit;
				ext_out[ext_out_n]->chrno = chr_re;
				ext_out[ext_out_n]->pos = pos_tmp;
				ext_out[ext_out_n]->dir = dir_tmp;
				strcpy(ext_out[ext_out_n]->cigar, cigar_tmp);

				ext_out_n++;


				if(tmp_i < ext_out_n1)
				{
					index1 = tmp_i;
					index2 = ext_out_n - 1;
					ext_out_n2_tmp = ext_out_n;
				}
				else
				{
					index2 = tmp_i - ext_out_n1;
					index1 = ext_out_n - 1;
					ext_out_n1_tmp = ext_out_n;
				}

				mat_index1[mat_n] = index1;
				mat_index2[mat_n] = index2;

				mat_n++;

			}
		}
	
		//
		//if(!find_anchor_flag)	tol_u++;	
	}

	ext_out_n1 = ext_out_n1_tmp;
	ext_out_n2 = ext_out_n2_tmp;

	for(tmp_i = 0; tmp_i < mat_n; tmp_i++)
	{
		mat_tmp1 = mat_index1[tmp_i];
		mat_tmp2 = mat_index2[tmp_i];


		l_j = ext_out1[mat_tmp1]->dir;
		r_j = ext_out2[mat_tmp2]->dir;

		dm = ext_out1[mat_tmp1]->edit + ext_out2[mat_tmp2]->edit;

		if(dm < dm_op)
		{
			for(psp_i = 0; psp_i < v_cnt; psp_i++)
			{
				ops_dm_ex1_p[psp_i] = op_dm_ex1_p[psp_i];
				ops_dm_ex2_p[psp_i] = op_dm_ex2_p[psp_i];
				ops_rc_p1[psp_i] = op_rc_p1[psp_i];
				ops_rc_p2[psp_i] = op_rc_p2[psp_i];
			}

			op_dm_ex1_p[0] = mat_tmp1;
			op_dm_ex2_p[0] = mat_tmp2;

			op_rc_p1[0] = l_j;
			op_rc_p2[0] = r_j;

			v_cnt = 1;
			dm_op = dm;
		}
		else if(dm == dm_op)
		{
			if(v_cnt < cus_max_output_ali)
			{
				op_dm_ex1_p[v_cnt] = mat_tmp1;
				op_dm_ex2_p[v_cnt] = mat_tmp2;

				op_rc_p1[v_cnt] = l_j;
				op_rc_p2[v_cnt] = r_j;

				v_cnt++;
			}
		}
		else if(dm < dm_ops)
		{
			ops_dm_ex1_p[0] = mat_tmp1;
			ops_dm_ex2_p[0] = mat_tmp2;

			ops_rc_p1[0] = l_j;
			ops_rc_p2[0] = r_j;

			vs_cnt = 1;
			dm_ops = dm;
		}
		else if(dm == dm_ops)
		{
			if(vs_cnt < cus_max_output_ali)
			{
				ops_dm_ex1_p[vs_cnt] = mat_tmp1;
				ops_dm_ex2_p[vs_cnt] = mat_tmp2;

				ops_rc_p1[vs_cnt] = l_j;
				ops_rc_p2[vs_cnt] = r_j;

				vs_cnt++;
			}
		}
	}

	if(mat_n == 0)
	{
		for(tmp_i = 0; tmp_i < ext_out_n1; tmp_i++)
		{
			dm = ext_out1[tmp_i]->edit;
			l_j = ext_out1[tmp_i]->dir;

			if(dm < dm_op)
			{
				for(psp_i = 0; psp_i < v_cnt; psp_i++)
				{
					ops_dm_ex1_p[psp_i] = op_dm_ex1_p[psp_i];
					ops_rc_p1[psp_i] = op_rc_p1[psp_i];
				}

				op_dm_ex1_p[0] = tmp_i;
				op_rc_p1[0] = l_j;

				v_cnt = 1;
				dm_op = dm;
			}
			else if(dm == dm_op)
			{
				if(v_cnt < cus_max_output_ali)
				{
					op_dm_ex1_p[v_cnt] = tmp_i;
					op_rc_p1[v_cnt] = l_j;

					v_cnt++;
				}
			}
			else if(dm < dm_ops)
			{
				ops_dm_ex1_p[0] = tmp_i;
				ops_rc_p1[0] = l_j;

				vs_cnt = 1;
				dm_ops = dm;
			}
			else if(dm == dm_ops)
			{
				if(vs_cnt < cus_max_output_ali)
				{
					ops_dm_ex1_p[vs_cnt] = tmp_i;
					ops_rc_p1[vs_cnt] = l_j;

					vs_cnt++;
				}
			}

		}

		if(ext_out_n1)
		{
			v_cnt1 = v_cnt;
			vs_cnt1 = vs_cnt;
		}

		dm_op = MAX_OP_SCORE;
		dm_ops = MAX_OP_SCORE;

		v_cnt = 0;
		vs_cnt = 0;

		for(tmp_i = 0; tmp_i < ext_out_n2; tmp_i++)
		{
			dm = ext_out2[tmp_i]->edit;
			r_j = ext_out2[tmp_i]->dir;

			if(dm < dm_op)
			{
				for(psp_i = 0; psp_i < v_cnt; psp_i++)
				{
					ops_dm_ex2_p[psp_i] = op_dm_ex2_p[psp_i];
					ops_rc_p2[psp_i] = op_rc_p2[psp_i];
				}

				op_dm_ex2_p[0] = tmp_i;
				op_rc_p2[0] = r_j;

				v_cnt = 1;
				dm_op = dm;
			}
			else if(dm == dm_op)
			{
				if(v_cnt < cus_max_output_ali)
				{
					op_dm_ex2_p[v_cnt] = tmp_i;
					op_rc_p2[v_cnt] = r_j;

					v_cnt++;
				}
			}
			else if(dm < dm_ops)
			{
				ops_dm_ex2_p[0] = tmp_i;
				ops_rc_p2[0] = r_j;

				vs_cnt = 1;
				dm_ops = dm;
			}
			else if(dm == dm_ops)
			{
				if(vs_cnt < cus_max_output_ali)
				{
					ops_dm_ex2_p[vs_cnt] = tmp_i;
					ops_rc_p2[vs_cnt] = r_j;

					vs_cnt++;
				}
			}
		}

		if(ext_out_n2)
		{
			v_cnt2 = v_cnt;
			vs_cnt2 = vs_cnt;
		}
	}
	else
	{
		v_cnt1 = v_cnt;
		v_cnt2 = v_cnt;

		vs_cnt1 = vs_cnt;
		vs_cnt2 = vs_cnt;
	}

	if(v_cnt1)
	{
		l_j = op_rc_p1[0];
		if(l_j)
		{
			for(read_b_i = 0; read_b_i < read_len1; read_b_i++)
				sam_seq[read_b_i] = Dna5Tochar[charToDna5n[seqio[seqi].read_seq1[read_b_i]] ^ 0X3];

			sam_seq[read_b_i] = '\0';

			strrev1(sam_seq);

			strcpy(read_rev_buffer[seqi], sam_seq);
			read_rev_buffer[seqi][read_len1] = '\0';

			seqio[seqi].seq1 = read_rev_buffer[seqi];

			strrev1(qual1_buffer[seqi]);

			seqio[seqi].qual1 = qual1_buffer[seqi];

		}
		else
		{
			seqio[seqi].seq1 = seqio[seqi].read_seq1;
			seqio[seqi].qual1 = qual1_buffer[seqi];
		}
	}
	if(v_cnt2)
	{
		r_j = op_rc_p2[0];
		if(r_j)
		{
			for(read_b_i = 0; read_b_i < read_len2; read_b_i++)
				sam_seq[read_b_i] = Dna5Tochar[charToDna5n[seqio[seqi].read_seq2[read_b_i]] ^ 0X3];

			sam_seq[read_b_i] = '\0';

			strrev1(sam_seq);

			strcpy(read_rev_buffer_1[seqi], sam_seq);
			read_rev_buffer_1[seqi][read_len1] = '\0';

			seqio[seqi].seq2 = read_rev_buffer_1[seqi];

			strrev1(qual2_buffer[seqi]);

			seqio[seqi].qual2 = qual2_buffer[seqi];
		}
		else
		{
			seqio[seqi].seq2 = seqio[seqi].read_seq2;
			seqio[seqi].qual2 = qual2_buffer[seqi];
		}
	}

	if(v_cnt1 && v_cnt2)
	{
		if(mat_n)
		{
			ext_p1 = ext_out1[op_dm_ex1_p[0]];//mat_index1[op_dm_ex1_p[0]]
			ext_p2 = ext_out2[op_dm_ex2_p[0]];//mat_index2[op_dm_ex2_p[0]]
			chr_re1 = ext_p1->chrno;
			chr_re2 = chr_re1;
		}
		else
		{
			ext_p1 = ext_out1[op_dm_ex1_p[0]];
			ext_p2 = ext_out2[op_dm_ex2_p[0]];
			chr_re1 = ext_p1->chrno;
			chr_re2 = ext_p2->chrno;
		}

		sam_pos1 = ext_p1->pos - chr_end_n[chr_re1 - 1];
		sam_pos2 = ext_p2->pos - chr_end_n[chr_re2 - 1];

		if((!l_j) && r_j)
		{
			if(mat_n)
			{
				sam_flag1 = 99;
				sam_flag2 = 147;
			}
			else
			{
				sam_flag1 = 97;
				sam_flag2 = 145;
			}
			sam_cross = sam_pos2 + read_len2 - sam_pos1;
		}
		else if(l_j && (!r_j))
		{
			if(mat_n)
			{
				sam_flag1 = 83;
				sam_flag2 = 163;
			}
			else
			{
				sam_flag1 = 81;
				sam_flag2 = 161;
			}

			sam_cross = sam_pos2 - read_len1 - sam_pos1;
		}
		else if(l_j && r_j)
		{
			if(mat_n)
			{
				sam_flag1 = 67;
				sam_flag2 = 131;
			}
			else
			{
				sam_flag1 = 65;
				sam_flag2 = 129;
			}

			sam_cross = sam_pos2 - read_len1 - sam_pos1;
		}
		else
		{
			if(mat_n)
			{
				sam_flag1 = 115;
				sam_flag2 = 179;
			}
			else
			{
				sam_flag1 = 113;
				sam_flag2 = 177;
			}

			sam_cross = sam_pos2 + read_len2 - sam_pos1;
		}

		edit1 = ext_p1->edit;
		edit2 = ext_p2->edit;

		cigar_tmp = ext_p1->cigar;
		cigar_tmp_o = ext_p2->cigar;
	}
	else if(v_cnt1)
	{
		ext_p1 = ext_out1[op_dm_ex1_p[0]];//mat_index1

		chr_re1 = ext_p1->chrno;
		chr_re2 = chr_re1;

		sam_pos1 = ext_p1->pos - chr_end_n[chr_re1 - 1];
		sam_pos2 = sam_pos1;

		if(!l_j)
		{
			sam_flag1 = 73;
			sam_flag2 = 133;
		}
		else
		{
			sam_flag1 = 89;
			sam_flag2 = 181;
		}
		sam_cross = 0;

		edit1 = ext_p1->edit;
		edit2 = 0;
		cigar_tmp = ext_p1->cigar;
		cigar_tmp_o = NULL;

		seqio[seqi].seq2 = seqio[seqi].read_seq2;

	}
	else if(v_cnt2)
	{
		ext_p2 = ext_out2[op_dm_ex2_p[0]];//mat_index2

		chr_re2 = ext_p2->chrno;
		chr_re1 = chr_re2;

		sam_pos2 = ext_p2->pos - chr_end_n[chr_re1 - 1];
		sam_pos1 = sam_pos2;

		if(!r_j)
		{
			sam_flag1 = 69;
			sam_flag2 = 137;
		}
		else
		{
			sam_flag1 = 117;
			sam_flag2 = 153;
		}
		sam_cross = 0;

		edit1 = 0;
		edit2 = ext_p2->edit;

		//fprintf(stderr, "3\n");

		cigar_tmp = NULL;
		cigar_tmp_o = ext_p2->cigar;

		seqio[seqi].seq1 = seqio[seqi].read_seq1;

	}
	else
	{
		sam_pos1 = 0;
		sam_pos2 = 0;

		chr_re1 = chr_file_n;
		chr_re2 = chr_re1;

		sam_flag1 = 77;
		sam_flag2 = 141;

		cigar_tmp = NULL;
		cigar_tmp_o = NULL;

		seqio[seqi].xa_n_p1 = 0;
		seqio[seqi].xa_n_p2 = 0;

		seqio[seqi].qualc1 = 0;
		seqio[seqi].qualc2 = 0;

		seqio[seqi].seq1 = seqio[seqi].read_seq1;
		seqio[seqi].seq2 = seqio[seqi].read_seq2;

		re_flag = 0;
		
		//
		//tol_u++;
	}

#ifdef	FIX_SV
	seqio[seqi].xa_n_x1 = 0;
	seqio[seqi].xa_n_x2 = 0;
#endif
	seqio[seqi].xa_n1 = 0;
	seqio[seqi].xa_n2 = 0;

	seqio[seqi].flag1 = sam_flag1;
	seqio[seqi].flag2 = sam_flag2;
	seqio[seqi].pos1 = sam_pos1;
	seqio[seqi].pos2 = sam_pos2;
	seqio[seqi].cross = sam_cross;
	seqio[seqi].chr_re1 = chr_re1;
	seqio[seqi].chr_re2 = chr_re2;
	seqio[seqi].nm1 = edit1;
	seqio[seqi].nm2 = edit2;

	if(cigar_tmp)	strcpy(pr_cigar1_buffer[seqi], cigar_tmp);
	else	strcpy(pr_cigar1_buffer[seqi], "*");
	seqio[seqi].cigar1 = pr_cigar1_buffer[seqi];

	if(cigar_tmp_o)	strcpy(pr_cigar2_buffer[seqi], cigar_tmp_o);
	else	strcpy(pr_cigar2_buffer[seqi], "*");
	seqio[seqi].cigar2 = pr_cigar2_buffer[seqi];

	if(mat_n)
	{
		for(tmp_i = 1, xa_i = 0; tmp_i < v_cnt1; tmp_i++, xa_i++)
		{
			ext_p1 = ext_out1[op_dm_ex1_p[tmp_i]];//mat_index1[op_dm_ex1_p[tmp_i]]

			chr_res[tid][xa_i] = ext_p1->chrno;

			sam_pos1s[tid][xa_i] = ext_p1->pos - chr_end_n[ext_p1->chrno - 1];

			strcpy(cigar_p1s[tid][xa_i], ext_p1->cigar);

			xa_d1s[tid][xa_i] = "+-"[ext_p1->dir];

			lv_re1s[tid][xa_i] = ext_p1->edit;
		}

		for(tmp_i = 0; tmp_i < vs_cnt1; tmp_i++, xa_i++)
		{
			ext_p1 = ext_out1[ops_dm_ex1_p[tmp_i]];//mat_index1[ops_dm_ex1_p[tmp_i]]

			chr_res[tid][xa_i] = ext_p1->chrno;

			sam_pos1s[tid][xa_i] = ext_p1->pos - chr_end_n[ext_p1->chrno - 1];;

			strcpy(cigar_p1s[tid][xa_i], ext_p1->cigar);

			xa_d1s[tid][xa_i] = "+-"[ext_p1->dir];

			lv_re1s[tid][xa_i] = ext_p1->edit;
		}
	}
	else
	{
		for(tmp_i = 1, xa_i = 0; tmp_i < v_cnt1; tmp_i++, xa_i++)
		{
			ext_p1 = ext_out1[op_dm_ex1_p[tmp_i]];

			chr_res[tid][xa_i] = ext_p1->chrno;

			sam_pos1s[tid][xa_i] = ext_p1->pos - chr_end_n[ext_p1->chrno - 1];

			strcpy(cigar_p1s[tid][xa_i], ext_p1->cigar);

			xa_d1s[tid][xa_i] = "+-"[ext_p1->dir];

			lv_re1s[tid][xa_i] = ext_p1->edit;
		}

		for(tmp_i = 0; tmp_i < vs_cnt1; tmp_i++, xa_i++)
		{
			ext_p1 = ext_out1[ops_dm_ex1_p[tmp_i]];

			chr_res[tid][xa_i] = ext_p1->chrno;

			sam_pos1s[tid][xa_i] = ext_p1->pos - chr_end_n[ext_p1->chrno - 1];;

			strcpy(cigar_p1s[tid][xa_i], ext_p1->cigar);

			xa_d1s[tid][xa_i] = "+-"[ext_p1->dir];

			lv_re1s[tid][xa_i] = ext_p1->edit;
		}
	}


	seqio[seqi].xa_n_p1 = xa_i;


	if(xa_i > 0)
	{
		memcpy(chr_res_buffer[seqi], chr_res[tid], xa_i << 2);
		seqio[seqi].chr_res = chr_res_buffer[seqi];

		memcpy(xa_d1s_buffer[seqi], xa_d1s[tid], xa_i);
		seqio[seqi].xa_d1s = xa_d1s_buffer[seqi];

		memcpy(sam_pos1s_buffer[seqi], sam_pos1s[tid], xa_i << 2);
		seqio[seqi].sam_pos1s = sam_pos1s_buffer[seqi];

		memcpy(lv_re1s_buffer[seqi], lv_re1s[tid], xa_i << 2);
		seqio[seqi].lv_re1s = lv_re1s_buffer[seqi];

		for(v_cnt_i = 0; v_cnt_i < xa_i; v_cnt_i++)
		{
			pos_n = strlen(cigar_p1s[tid][v_cnt_i]) + 1;
			seqio[seqi].cigar_p1s[v_cnt_i] = (char* )malloc(pos_n);

			memcpy(seqio[seqi].cigar_p1s[v_cnt_i], cigar_p1s[tid][v_cnt_i], pos_n);
		}
	}


	if(mat_n)
	{
		for(tmp_i = 1, xa_i = 0; tmp_i < v_cnt2; tmp_i++, xa_i++)
		{
			ext_p2 = ext_out2[op_dm_ex2_p[tmp_i]];//mat_index2[op_dm_ex2_p[tmp_i]]

			chr_res[tid][xa_i] = ext_p2->chrno;

			sam_pos2s[tid][xa_i] = ext_p2->pos - chr_end_n[ext_p2->chrno - 1];;

			strcpy(cigar_p2s[tid][xa_i], ext_p2->cigar);

			xa_d2s[tid][xa_i] = "+-"[ext_p2->dir];

			lv_re2s[tid][xa_i] = ext_p2->edit;
		}

		for(tmp_i = 0; tmp_i < vs_cnt2; tmp_i++, xa_i++)
		{
			ext_p2 = ext_out2[ops_dm_ex2_p[tmp_i]];//mat_index2[ops_dm_ex2_p[tmp_i]]

			chr_res[tid][xa_i] = ext_p2->chrno;

			sam_pos2s[tid][xa_i] = ext_p2->pos - chr_end_n[ext_p2->chrno - 1];;

			strcpy(cigar_p2s[tid][xa_i], ext_p2->cigar);

			xa_d2s[tid][xa_i] = "+-"[ext_p2->dir];

			lv_re2s[tid][xa_i] = ext_p2->edit;
		}
	}
	else
	{
		for(tmp_i = 1, xa_i = 0; tmp_i < v_cnt2; tmp_i++, xa_i++)
		{
			ext_p2 = ext_out2[op_dm_ex2_p[tmp_i]];

			chr_res[tid][xa_i] = ext_p2->chrno;

			sam_pos2s[tid][xa_i] = ext_p2->pos - chr_end_n[ext_p2->chrno - 1];;

			strcpy(cigar_p2s[tid][xa_i], ext_p2->cigar);

			xa_d2s[tid][xa_i] = "+-"[ext_p2->dir];

			lv_re2s[tid][xa_i] = ext_p2->edit;
		}

		for(tmp_i = 0; tmp_i < vs_cnt2; tmp_i++, xa_i++)
		{
			ext_p2 = ext_out2[ops_dm_ex2_p[tmp_i]];

			chr_res[tid][xa_i] = ext_p2->chrno;

			sam_pos2s[tid][xa_i] = ext_p2->pos - chr_end_n[ext_p2->chrno - 1];;

			strcpy(cigar_p2s[tid][xa_i], ext_p2->cigar);

			xa_d2s[tid][xa_i] = "+-"[ext_p2->dir];

			lv_re2s[tid][xa_i] = ext_p2->edit;
		}
	}


	seqio[seqi].xa_n_p2 = xa_i;

	if(xa_i > 0)
	{
		if(xa_i > seqio[seqi].xa_n_p1)
		{
			memcpy(chr_res_buffer[seqi], chr_res[tid], xa_i << 2);
			seqio[seqi].chr_res = chr_res_buffer[seqi];
		}
		memcpy(xa_d2s_buffer[seqi], xa_d2s[tid], xa_i);
		seqio[seqi].xa_d2s = xa_d2s_buffer[seqi];

		memcpy(sam_pos2s_buffer[seqi], sam_pos2s[tid], xa_i << 2);
		seqio[seqi].sam_pos2s = sam_pos2s_buffer[seqi];

		memcpy(lv_re2s_buffer[seqi], lv_re2s[tid], xa_i << 2);
		seqio[seqi].lv_re2s = lv_re2s_buffer[seqi];

		for(v_cnt_i = 0; v_cnt_i < xa_i; v_cnt_i++)
		{
			pos_n = strlen(cigar_p2s[tid][v_cnt_i]) + 1;
			seqio[seqi].cigar_p2s[v_cnt_i] = (char* )malloc(pos_n);

			memcpy(seqio[seqi].cigar_p2s[v_cnt_i], cigar_p2s[tid][v_cnt_i], pos_n);
		}
	}

	if((v_cnt1 == 1) || (seqio[seqi].xa_n_p1 == 0))
	{
		if(seqio[seqi].xa_n_p1)
			seqio[seqi].qualc1 = 20;
		else	seqio[seqi].qualc1 = 60;
	}
	else	seqio[seqi].qualc1 = 0;

	if((v_cnt2 == 1) || (seqio[seqi].xa_n_p2 == 0))
	{
		if(seqio[seqi].xa_n_p2)
			seqio[seqi].qualc2 = 20;
		else	seqio[seqi].qualc2 = 60;
	}
	else	seqio[seqi].qualc2 = 0;

	return re_flag;
}

int compare_ext_out(const void * a, const void * b)
{
	ext_output** ext1 = (ext_output** )a;
	ext_output** ext2 = (ext_output** )b;

	uint32_t pos1 = (*ext1)->pos;
	uint32_t pos2 = (*ext2)->pos;

#ifdef	DEBUG_EXT_ALN
	fprintf(stderr, "com: %u %u\n", pos1, pos2);
#endif

	if(pos1 > pos2)
		return 1;
	else if(pos1 < pos2)
		return -1;
	else
		return 0;
}

int binsearch_pair_pos_ext(uint32_t x, ext_output** bp, uint16_t n, int devi)
{
	int low, high, mid;

	low = 0;
	high = n - 1;

	while ( low <= high )
	{
		mid = (low + high) >> 1;

		if(x < (bp[mid]->pos - devi))
		{
			high = mid - 1;
		}
		else if(x > (bp[mid]->pos + devi))
		{
			low = mid + 1;
		}
		else     /*found match*/
		{
			return mid;
		}
	}

	return -1;
}
uint16_t ext_single_tree_alt(uint8_t** seq_test_b, uint8_t** seq_test_b_var, uint8_t* pattern, uint8_t** seqss, ext_results** ex_results, com_result** c_results, uint16_t* re_c_result_n, uint8_t thread_i, uint32_t* extension_start_pos_array, uint16_t seq_n, uint16_t seq_length, uint16_t read_len, uint8_t dir, uint16_t k, uint64_t* buffer_ref_seq_var, uint8_t alt_flag)
{
	uint8_t seq_char = 0;
	uint8_t seq_char_d = 0;
	uint8_t seq_char_var = 0;
	uint8_t continue_f = 0;
	uint8_t next_nodes = 0;
	uint8_t vcf_pos_f = 0;
	uint8_t alt_type = 0;
	uint8_t route_rsn = 0;
	uint8_t r_node_no = 0;
	uint8_t num_tmp = 0;
	uint8_t off_tmp = 0;
	uint8_t cigar_tmp = 0;
	uint8_t cigar_tmp_pre = 0;

	uint16_t bucket[BUCKET_SIZE][EXTENSION_REF_LENGTH];
	uint16_t bucket_char_n[BUCKET_SIZE];//16
	uint16_t bucket_char_n_tmp = 0;
	uint16_t tmp_seqn_i = 0;
	uint16_t tmp_seqn_j = 0;
	uint16_t seq_i = 0;
	uint16_t ext_i = 0;
	uint16_t bucket_char_n_i = 0;
	uint16_t node_array_n = 0;
	uint16_t node_array_n_c = 0;
	uint16_t seqs_n = 0;
	uint16_t node_array_i = 0;
	uint16_t ext_n = 0;
	uint16_t seq_nums_n_tmp = 0;
	uint16_t seq_nums_tmp = 0;
	uint16_t seq_len_n_tmp = 0;
	uint16_t tmp_j = 0;
	uint16_t route_arrn_i = 0;
	uint16_t route_i = 0;
	uint16_t alt_info_cnt = 0;
	uint16_t node_offset = 0;
	uint16_t ending_node_offset = 0;
	uint16_t node_id = 0;
	uint16_t ending_node_id = 0;
	uint16_t alt_nums = 0;
	uint16_t alt_length = 0;
	uint16_t alt_info_cnt_pre = 0;
	uint16_t cnt_pre = 0, cnt = 0;
	uint16_t route_no = 0;
	uint16_t sn = 0, snt = 0;
	uint16_t pdis_tmp = 0;
	uint16_t sub_tmp = 0;
	uint16_t tail_tmp = 0;
	uint16_t seq_id = 0;
	uint16_t edit_v = 0;
	uint16_t c_ren = 0, ex_ren = 0;
	uint16_t altseq_n = 0;
	uint16_t re_cnt = 0;
	uint16_t nodeid_tmp = 0;
	uint16_t nodeid_tmp_pre = 0Xffff;
	uint16_t tmp_alt_n = 0;

	int16_t rs_i = 0;

	uint32_t qnum = 0;
	uint32_t cal_n = 0;
	uint32_t cal_arrn = 0;
	uint32_t cal_arr = 0;
	uint32_t ref_pos_end = 0;
	uint32_t route_arr_i = 0;
	uint32_t tmp_i = 0;
	uint32_t aindex_i = 0;
	uint32_t ref_index = 0;
	uint32_t alt_start = 0;
	uint32_t alt_end = 0;
	uint32_t alt_seq_begin = 0;
	uint32_t alt_info_cnt_tal = 0;
	uint32_t tmp_alt_info_cnt_tal = 0;
	uint32_t ref_pos_value = 0;
	uint32_t chr_re = 0;
	uint32_t route_id_tmp = 0;
	uint32_t route_rs[MAX_K_VAR];
	uint8_t cigar_rs[MAX_K_VAR];

	int tmp_node_i = 0;
	int tmp_node_j = 0;//uint16_t

	int64_t low = 0, high = 0, mid = 0;

	char cigars[MAX_CIGAR_LEN];
	uint8_t* seqs = NULL;
	uint8_t* aseq_array_p = NULL;
	uint8_t* type_array_p = NULL;

	uint16_t* seq_nums_p = NULL;
	uint16_t** seq_node = NULL;
	uint16_t* seq_node_n = NULL;
	uint16_t** seq_node_len = NULL;
	uint16_t* seq_node_len_n = NULL;
	uint16_t** noder_array_tree = NULL;
	uint16_t* nodern_array_tree = NULL;
	uint16_t* route_n_index_array = NULL;
	uint16_t* result_seqpos_array = NULL;
	uint16_t* route_sep = NULL;
	uint16_t* route_arr = NULL;
	uint32_t* route_index = NULL;
	uint32_t* route_arrn = NULL;
	uint32_t* route_id_mark = NULL;
	uint32_t* route_id_mark_pre = NULL;
	uint32_t* route_nodes = NULL;
	uint32_t chr_nos[MAX_EXT_ALT_POS];
	uint32_t* rpos_array_p = NULL;
	uint32_t* aindex_array_p = NULL;
	uint32_t* rpos_end_array_p = NULL;
	uint32_t* chr_rpos_a_p = NULL;

	uint64_t* chr_end_n_p = NULL;
	uint32_t* apos_array_p = NULL;

	treenode* t_node = NULL;
	treenode* t_node_c = NULL;
	treenode* tmp_node_p = NULL;
	treenode** node_array = NULL;
	alt_info_a* alt_info_array = NULL;
	ext_queue** lv_queue = NULL;
	ext_result** lv_result = NULL;
	ext_results* ex_result = NULL;
	com_result* c_result = NULL;
	route_node* tmp_r_node = NULL;
	route_node** r_node_p = NULL;
	alt_seq_a* alt_seqs = NULL;

	seqs = seqss[thread_i];
	ex_result = ex_results[thread_i];
	c_result = c_results[thread_i];
	node_array = node_arrays[thread_i];
	seq_node = seq_nodes[thread_i];
	seq_node_n = seq_node_ns[thread_i];
	seq_node_len = seq_node_lens[thread_i];
	seq_node_len_n = seq_node_len_ns[thread_i];
	alt_info_array = alt_info_arrays[thread_i];
	lv_queue = lv_queues[thread_i];
	lv_result = lv_results[thread_i];
	route_sep = route_seps[thread_i];
	alt_seqs = alt_seqss[thread_i];

	route_id_mark = route_id_marks[thread_i];
	route_id_mark_pre = route_id_marks_pre[thread_i];

	//create node
	t_node = (treenode* )calloc(1, sizeof(treenode));

	if(t_node == NULL)	exit(1);


	t_node->seq_nums_n = seq_n;
	for(bucket_char_n_i = 0; bucket_char_n_i < seq_n; bucket_char_n_i++)
	{
		t_node->seq_nums[bucket_char_n_i] = bucket_char_n_i;
		//to connect ALT
		seq_node[bucket_char_n_i][0] = 0;
		seq_node_n[bucket_char_n_i] = 1;
	}

	memset(seq_node_len_n, 0, seq_n << 1);

	seqs_n = 0;
	node_array[0] = t_node;
	node_array_n = 1;

	for(node_array_i = 0; node_array_i < node_array_n; node_array_i++)
	{
		t_node_c = node_array[node_array_i];
		t_node_c->alt_nums = 0;
		t_node_c->seqs_start = seqs_n;
		seq_nums_n_tmp = t_node_c->seq_nums_n;
		seq_nums_p = t_node_c->seq_nums;

		if(seq_nums_n_tmp == 1)
		{
			t_node_c->seqs_length = seq_length - t_node_c->seq_start;
			//memcpy(seqs + seqs_n, seq_test_b[seq_nums_p[0]] + t_node_c->seq_start, t_node_c->seqs_length);
			memcpy(seqs + seqs_n, seq_test_b_var[seq_nums_p[0]] + t_node_c->seq_start, t_node_c->seqs_length);
			seqs_n += t_node_c->seqs_length;
		}
		else
		{
			for(ext_i = t_node_c->seq_start, ext_n = 0; ext_i < seq_length; ext_i++, ext_n++)
			{
				continue_f = 1;

				seq_char_d = seq_test_b[seq_nums_p[0]][ext_i];
				seq_char_var = seq_test_b_var[seq_nums_p[0]][ext_i];

				memset(bucket_char_n, 0, BUCKET_SIZE << 1);

				for(seq_i = 0; seq_i < seq_nums_n_tmp; seq_i++)   //seq_n
				{
					seq_char = seq_test_b[seq_nums_p[seq_i]][ext_i];
					bucket[seq_char][bucket_char_n[seq_char]++] = seq_nums_p[seq_i];
					if(seq_char != seq_char_d)	continue_f = 0;
				}

				if(continue_f)
				{
					//seqs[seqs_n++] = seq_char_d;//connect all nodes' seqs together
					seqs[seqs_n++] = seq_char_var;
					continue;
				}
				else
				{
					break;
				}
			}

			t_node_c->seqs_length = ext_n;

			if(ext_i != seq_length)
			{
				node_array_n_c = node_array_n;
				for(bucket_char_n_i = 0, next_nodes = 0; bucket_char_n_i < BUCKET_SIZE; bucket_char_n_i++)
				{
					bucket_char_n_tmp = bucket_char_n[bucket_char_n_i];
					if(bucket_char_n_tmp != 0)
					{
						//create node
						t_node = (treenode* )calloc(1, sizeof(treenode));
						memcpy(t_node->seq_nums, bucket[bucket_char_n_i], bucket_char_n_tmp << 1);

						//to connect ALT
						for(tmp_i = 0; tmp_i < bucket_char_n_tmp; tmp_i++)
							seq_node[t_node->seq_nums[tmp_i]][seq_node_n[t_node->seq_nums[tmp_i]]++] = node_array_n;

						t_node->seq_start = ext_i;
						t_node->seq_nums_n = bucket_char_n_tmp;

						//
						t_node->seq_dis += (t_node_c->seq_dis + ext_n);

						node_array[node_array_n++] = t_node;

						next_nodes++;
					}
				}


				t_node_c->next_nodes_n = next_nodes;
				//add next nodes
				for(bucket_char_n_i = 0; bucket_char_n_i < next_nodes; bucket_char_n_i++)
					t_node_c->next_nodes[bucket_char_n_i] = node_array_n_c + bucket_char_n_i;
			}
		}

		//
		//t_node_c->seq_dis += t_node_c->seqs_length;

		//to connect ALT
		//seq_node_len records each seperation by each node
		for(tmp_i = 0; tmp_i < seq_nums_n_tmp; tmp_i++)
		{
			seq_nums_tmp = seq_nums_p[tmp_i];
			seq_len_n_tmp = seq_node_len_n[seq_nums_tmp];
			if(seq_len_n_tmp)	seq_node_len[seq_nums_tmp][seq_len_n_tmp] = seq_node_len[seq_nums_tmp][seq_len_n_tmp - 1] + t_node_c->seqs_length;
			else	seq_node_len[seq_nums_tmp][seq_len_n_tmp] = t_node_c->seqs_length;

			seq_node_len_n[seq_nums_tmp]++;
		}
	}

	alt_info_cnt = 1;
	route_sep[0] = 1;

	if(dir)
	{
		rpos_array_p = rpos_array;
		aindex_array_p = aindex_array;
		apos_array_p = apos_array;
		rpos_end_array_p = rpos_end_array;
		aseq_array_p = aseq_array;
		type_array_p = type_array;
		chr_end_n_p = chr_end_n;
		chr_rpos_a_p = chr_rpos_a;
	}
	else
	{
		rpos_array_p = rpos_array_rec;
		aindex_array_p = aindex_array_rec;
		apos_array_p = apos_array_rec;
		rpos_end_array_p = rpos_end_array_rec;
		aseq_array_p = aseq_array_rec;
		type_array_p = type_array_rec;
		chr_end_n_p = chr_end_n_rec;
		chr_rpos_a_p = chr_rpos_a_rec;
	}
	
	uint16_t del_tol = 0;
	uint16_t seq_n_tmp = 0;
	uint16_t seq_len_tmp = 0;
	uint16_t seq_id_index[MAX_EXT_ALT_POS];
	uint32_t read_b_i = 0;
	uint32_t ext_len_tmp = 0;
	uint32_t seqs_cp_s = 0;
	uint32_t seqs_cp_e = 0;
	
	for(tmp_seqn_i = 0; tmp_seqn_i < seq_n; tmp_seqn_i++)
		seq_id_index[tmp_seqn_i] = tmp_seqn_i;


	seq_n_tmp = seq_n;
	for(tmp_seqn_j = 0; (tmp_seqn_j < seq_n_tmp) && alt_flag; tmp_seqn_j++)
	{
		ref_pos_value = extension_start_pos_array[tmp_seqn_j];//can be not in order
		tmp_seqn_i = seq_id_index[tmp_seqn_j];

		low = 0;
		high = chr_file_n - 1;
		while ( low <= high )
		{
			mid = (low + high) >> 1;
			if(ref_pos_value < (chr_end_n_p[mid]))
			{
				high = mid - 1;
			}
			else if(ref_pos_value > (chr_end_n_p[mid]))
			{
				low = mid + 1;
			}
			else
			{
				chr_re =  mid;
				break;
			}
			chr_re = low;
		}

		chr_nos[tmp_seqn_i] = chr_re;

		vcf_pos_f = 0;
		low = chr_rpos_a_p[chr_re - 1];//0
		high = chr_rpos_a_p[chr_re] - 1;//rpos_cnt - 1
		
		if(high == -1)	continue;

		while ( low <= high )
		{
			mid = ((int64_t )(low + high)) >> 1;
			if(ref_pos_value < rpos_array_p[mid])   //absolute pos (including chr)
			{
				high = mid - 1;
			}
			else if(ref_pos_value > rpos_array_p[mid])
			{
				low = mid + 1;
			}
			else     //found match
			{
				ref_index = mid;
				vcf_pos_f = 1;
				break;
			}
		}

		if(!vcf_pos_f)	ref_index = high + 1; //what if ref pos far away from ALTs pos ?
		if(high == chr_rpos_a_p[chr_re])	continue;	//rpos_cnt break;

		del_tol = 0;

		while(1)
		{
			alt_start = rpos_array_p[ref_index];
			for(tmp_node_i = seq_node_n[tmp_seqn_i] - 1; tmp_node_i >= 0; tmp_node_i--)
			{
				if(seq_node_len[tmp_seqn_i][tmp_node_i] < alt_start - ref_pos_value + 1)   //rpos_array_p[ref_index]
				{
					tmp_node_i++;
					break;
				}
			}
			if(tmp_node_i == seq_node_n[tmp_seqn_i])
				break;

			node_id = seq_node[tmp_seqn_i][tmp_node_i];//which node the alt is located at

			//node_array[node_id]->seqs_start
			//node_offset is the offset of this alt to its previous node
			if(tmp_node_i)	node_offset = alt_start - ref_pos_value - seq_node_len[tmp_seqn_i][tmp_node_i - 1];//rpos_array_p[ref_index]
			else	node_offset = alt_start - ref_pos_value;//rpos_array_p[ref_index]

			alt_nums = aindex_array_p[ref_index + 1] - aindex_array_p[ref_index];

			node_array[node_id]->node_offset = node_offset;//find the position of alts at this node, maybe not be used?
			node_array[node_id]->alt_nums += alt_nums;//check whether this node has ALT and how many ALTs at the same position
			//node_array[node_id]->alt_info_cnt = alt_info_cnt;

			for(aindex_i = aindex_array_p[ref_index]; aindex_i < aindex_array_p[ref_index + 1]; aindex_i++)
			{
				alt_seq_begin = apos_array_p[aindex_i];
				alt_length = apos_array_p[aindex_i + 1] - apos_array_p[aindex_i];

				//find ending node ID and the offset

				alt_end = rpos_end_array_p[aindex_i];

				if(type_array_p[aindex_i] != 1)   //not I
				{
					alt_type = 0;//for debug

					for(tmp_node_j = seq_node_n[tmp_seqn_i] - 1; tmp_node_j >= 0; tmp_node_j--)
					{
						if(seq_node_len[tmp_seqn_i][tmp_node_j] < alt_end - ref_pos_value + 1)
						{
							tmp_node_j++;
							break;
						}
					}

					if(tmp_node_j != seq_node_n[tmp_seqn_i])   //found ending node
					{
						ending_node_id = seq_node[tmp_seqn_i][tmp_node_j];
						if(tmp_node_j)	ending_node_offset = alt_end - ref_pos_value - seq_node_len[tmp_seqn_i][tmp_node_j - 1];
						else	ending_node_offset = alt_end - ref_pos_value;

						del_tol += (alt_end - alt_start);

					}
					else
					{
						t_node = (treenode* )calloc(1, sizeof(treenode));
						t_node->seq_nums[0] = seq_n;
						t_node->seq_nums_n = 1;
						
						t_node->seqs_start = seqs_n;
						ext_len_tmp = seq_length - (alt_start - ref_pos_value);
						t_node->seqs_length = ext_len_tmp;//del_tol
						
						if(dir)
						{
							seqs_cp_s = alt_end;
							seqs_cp_e = alt_end + ext_len_tmp;
							for(read_b_i = seqs_cp_s, ext_i = 0; read_b_i < seqs_cp_e; read_b_i++, ext_i++)//START_POS_REF
							{
								seqs[seqs_n + ext_i] = ((buffer_ref_seq_var[read_b_i >> 4] >> ((15 - (read_b_i & 0Xf)) << 2)) & 0Xf);

							}
						}else{
							seqs_cp_s = ref_length_p - alt_end;
							seqs_cp_e = ref_length_p - (alt_end + ext_len_tmp); 							
							for(read_b_i = seqs_cp_s, ext_i = 0; read_b_i > seqs_cp_e; read_b_i--, ext_i++)//START_POS_REF
							{
								seqs[seqs_n + ext_i] = ((buffer_ref_seq_var[read_b_i >> 4] >> ((15 - (read_b_i & 0Xf)) << 2)) & 0Xf);
							}
						}
					
						seqs_n += ext_len_tmp;//t_node_c->seqs_length;

						seq_node_len[seq_n][0] = ext_len_tmp;
						seq_node_len_n[seq_n] = 1;

						seq_node[seq_n][0] = node_array_n;
						seq_node_n[seq_n] = 1;

						node_array[node_array_n] = t_node;
						
						//add new pos into extension_start_pos_array[tmp_seqn_i];
						extension_start_pos_array[seq_n_tmp] = alt_end;
						seq_id_index[seq_n_tmp] = seq_n;
						seq_n_tmp++;
						seq_n++;

						ending_node_id = node_array_n;
						ending_node_offset = 0;
						
						node_array_n++;
					}
				}
				else
				{
					alt_type = 1;//for debug
					ending_node_id = node_id;
					ending_node_offset = node_offset;
				}

				alt_info_array[alt_info_cnt].nodeid = node_id;
				alt_info_array[alt_info_cnt].node_offset = node_offset;
				alt_info_array[alt_info_cnt].alt_length = alt_length;		//0 for DEL
				alt_info_array[alt_info_cnt].ending_node_id = ending_node_id;
				alt_info_array[alt_info_cnt].ending_node_offset = ending_node_offset;
				alt_info_array[alt_info_cnt].alt_seq_begin = alt_seq_begin;
				alt_info_array[alt_info_cnt].alt_start = alt_start - ref_pos_value;
				alt_info_array[alt_info_cnt].alt_end = alt_end - ref_pos_value;
				alt_info_array[alt_info_cnt].vcf_index = ref_index;//not be used, for debug ?
				alt_info_array[alt_info_cnt].alt_type = alt_type;
				//add route ids
				alt_info_array[alt_info_cnt].ori_index = alt_info_cnt;
				//
				alt_info_array[alt_info_cnt].seq_index = tmp_seqn_i;

				//new
				alt_seqs[alt_info_cnt].ref_index = ref_index;
				alt_seqs[alt_info_cnt].aindex = aindex_i;

				alt_info_cnt++;
			}

			//node_array[node_id]->alt_p = (altnode* )calloc(1, sizeof(altnode ));
			ref_index++;
		}

		if(del_tol > MAX_EDIT) //extend last node with length of del_tol
		{
			//crate a new node with length of ... and connect it to its previous node node_array[] seq_node seq_node_n seq_node_len and add next

			t_node = (treenode* )calloc(1, sizeof(treenode));
			t_node->seq_nums[0] = tmp_seqn_i;
			t_node->seq_nums_n = 1;
			t_node->seqs_start = seqs_n;
			t_node->seqs_length = del_tol;

			seq_len_tmp = seq_node_len[tmp_seqn_i][seq_node_len_n[tmp_seqn_i] - 1];
			tmp_node_j = seq_node_n[tmp_seqn_i];

			seq_node_len[seq_n][0] = del_tol;
			seq_node_len_n[seq_n] = 1;

			seq_node[seq_n][0] = node_array_n;
			seq_node_n[seq_n] = 1;
			
			if(dir)
			{
				seqs_cp_s = ref_pos_value + seq_len_tmp;
				seqs_cp_e = ref_pos_value + seq_len_tmp + del_tol;
				for(read_b_i = seqs_cp_s, ext_i = 0; read_b_i < seqs_cp_e; read_b_i++, ext_i++)
				{
					seqs[seqs_n + ext_i] = ((buffer_ref_seq_var[read_b_i >> 4] >> ((15 - (read_b_i & 0Xf)) << 2)) & 0Xf);
				}
			}else{
				seqs_cp_s = ref_length_p - (ref_pos_value + seq_len_tmp);
				seqs_cp_e = ref_length_p - (ref_pos_value + seq_len_tmp + del_tol);
				for(read_b_i = seqs_cp_s, ext_i = 0; read_b_i < seqs_cp_e; read_b_i--, ext_i++)
				{
					seqs[seqs_n + ext_i] = ((buffer_ref_seq_var[read_b_i >> 4] >> ((15 - (read_b_i & 0Xf)) << 2)) & 0Xf);
				}				
			}

			seqs_n += del_tol;//t_node_c->seqs_length;

			node_array[seq_node[tmp_seqn_i][tmp_node_j - 1]]->next_nodes[0] = node_array_n;
			node_array[seq_node[tmp_seqn_i][tmp_node_j - 1]]->next_nodes_n = 1;
			node_array[node_array_n] = t_node;
			node_array_n++;

			//add new pos into extension_start_pos_array[tmp_seqn_i];
			extension_start_pos_array[seq_n_tmp] = ref_pos_value + seq_len_tmp;
			seq_id_index[seq_n_tmp] = seq_n;//tmp_seqn_i
			seq_n_tmp++;
			seq_n++;
		}

		route_sep[tmp_seqn_i + 1] = alt_info_cnt;
	}

	seq_n_tmp = seq_n;

	//add route ids
	for(tmp_seqn_i = 0; tmp_seqn_i < seq_n; tmp_seqn_i++)
	{
		cal_n = route_sep[tmp_seqn_i + 1] - route_sep[tmp_seqn_i];

		cal_arrn += (pow(2, cal_n) - 1);

		cal_arr += (uint32_t )(pow((double)2, (double)(cal_n - 1)) * (double)cal_n);
	}
	cal_arrn++;
	cal_arr++;

#ifdef	DEBUG_EXT_ALN
	fprintf(stderr, "there at most %u routes %u\n", cal_arrn, cal_arr);
#endif

	route_arr = (uint16_t* )calloc(cal_arr, 2);
	route_arrn = (uint32_t* )calloc(cal_arrn, 4);
	route_index = (uint32_t* )calloc(cal_arrn, 4);
	route_nodes = (uint32_t* )calloc(alt_info_cnt + 1, 4);

	//debug
	if(route_arr == NULL)	exit(1);
	if(route_arrn == NULL)	exit(1);
	if(route_index == NULL)	exit(1);

	qnum = 1;
	route_index[0] = 1;

	for(tmp_seqn_i = 0; tmp_seqn_i < seq_n; tmp_seqn_i++)
	{
		cnt_pre = route_sep[tmp_seqn_i];
		cnt = route_sep[tmp_seqn_i + 1];

		if(cnt > cnt_pre)
		{
			create_route_array(cnt_pre, cnt, &route_arr, &route_arrn, &qnum, &route_index, route_nodes, alt_info_array);
		}
	}

	//combine alt_nums and node_offset to traverse ALT seq on tree
	qsort(alt_info_array + 1, alt_info_cnt - 1, sizeof(alt_info_a), compare_alt_sort);

	alt_info_array[alt_info_cnt].nodeid = node_array_n;//to end while

	//indeed need to be sorted
	for(tmp_i = 1; tmp_i < alt_info_cnt; tmp_i++)
	{
		nodeid_tmp = alt_info_array[tmp_i].nodeid;
		if(nodeid_tmp != nodeid_tmp_pre)
		{
			node_array[nodeid_tmp]->alt_info_cnt = tmp_i;//alt_info_cnt records the separation of each node according to its ALTs
			if(nodeid_tmp_pre != 0Xffff)
				node_array[nodeid_tmp_pre]->alt_info_n = tmp_alt_n;
			tmp_alt_n = 0;
		}
		++tmp_alt_n;
		nodeid_tmp_pre = nodeid_tmp;
	}
	node_array[nodeid_tmp]->alt_info_n = tmp_alt_n;

	computeEditDistance_tree(node_array,
	                         pattern,
	                         read_len,
	                         k,
	                         lv_queue,
	                         &lv_result,
	                         &re_cnt,
	                         alt_info_array,//alt_info_array_tmp,
	                         aseq_array_p,
	                         seqs,
	                         qnum,//
							 route_id_mark,
							 route_id_mark_pre,
	                         rnode_array[thread_i],
	                         routeid_arrays[thread_i],
	                         route_arr,
	                         route_arrn,
	                         route_index,
	                         node_array_n,
	                         route_nodes,
	                         route_sep
	                        );	

	if(route_index)	free(route_index);

	r_node_p = rnode_array[thread_i];

	for(tmp_i = 0; tmp_i < re_cnt; tmp_i++)
	{
		edit_v = lv_result[tmp_i]->edit;

		r_node_no = lv_result[tmp_i]->rnode_index;
		if(r_node_no)
		{
			snt = 0;
			sub_tmp = 0;
			num_tmp = 1;
			route_rsn = 0;

			cigar_rs[0] = lv_result[tmp_i]->cigar;

			do
			{
				route_rs[route_rsn] = r_node_p[r_node_no]->pdis;
				route_rsn++;

				tmp_r_node = r_node_p[r_node_no];

				if(tmp_r_node->cigar < 4)	cigar_rs[route_rsn] = tmp_r_node->cigar;

				r_node_no = tmp_r_node->pre_index;
			}
			while(r_node_no);

			for(rs_i = route_rsn - 1; rs_i > -1; rs_i--)
			{
				pdis_tmp = route_rs[rs_i];
				cigar_tmp = cigar_rs[rs_i];
			}

			if(route_rsn > 1)
			{
				for(rs_i = route_rsn - 1; rs_i > 0; rs_i--)
				{
					pdis_tmp = route_rs[rs_i];// + (num_tmp - 1)
					cigar_tmp = cigar_rs[rs_i];

					if((cigar_tmp == cigar_rs[rs_i - 1]) && (route_rs[rs_i - 1] - pdis_tmp == 1 - (cigar_tmp >> 1)))//cigar_rs[rs_i]
						num_tmp++;
					else
					{
						off_tmp = ((cigar_tmp == 2) ? 0:num_tmp - 1);
						if(pdis_tmp - sub_tmp > off_tmp)   //
						{
							sn = sprintf(cigars + snt, "%uM", pdis_tmp - sub_tmp - off_tmp);// - (num_tmp - 1)
							snt += sn;
						}
						sn = sprintf(cigars + snt, "%u%c", num_tmp, "XID"[cigar_tmp]);
						snt += sn;

						sub_tmp = pdis_tmp + ((cigar_tmp == 2) ? 0:1);//num_tmp
						num_tmp = 1;
					}
				}
			}
			else
			{
				pdis_tmp = route_rs[0];
			}

			if(read_len - pdis_tmp)
			{
				pdis_tmp = route_rs[0];//num_tmp - 1
				cigar_tmp = cigar_rs[0];

				off_tmp = ((cigar_tmp == 2) ? 0:num_tmp - 1);
				if(pdis_tmp - sub_tmp > off_tmp)   //num_tmp - 1
				{
					sn = sprintf(cigars + snt, "%uM", pdis_tmp - sub_tmp - off_tmp);// - (num_tmp - 1)
					snt += sn;
				}
				sn = sprintf(cigars + snt, "%u%c", num_tmp, "XID"[cigar_tmp]);
				snt += sn;

				tail_tmp = pdis_tmp + ((cigar_tmp == 2) ? 0:1);
				if(tail_tmp < read_len)
				{
					sn = sprintf(cigars + snt, "%uM", read_len - tail_tmp);//read_len - 1 - pdis_tmp num_tmp
					snt += sn;
				}
			}
		}
		else
		{
			sprintf(cigars, "%uM", read_len);
			snt = read_len;
		}
		cigars[snt] = '\0';

		//For seqs
		tmp_node_p = node_array[lv_result[tmp_i]->nodeid];

		if(tmp_node_p->re_flag)	continue;

		tmp_node_p->re_flag = 1;
		for(tmp_j = 0; (tmp_j < tmp_node_p->seq_nums_n) && (ex_ren < TREE_EXT_OUTPUT); tmp_j++)
		{
			seq_id = tmp_node_p->seq_nums[tmp_j];

			ex_result[ex_ren].edit = edit_v;
			if(dir)	ex_result[ex_ren].chrno = chr_nos[seq_id];
			else	ex_result[ex_ren].chrno = chr_file_n - chr_nos[seq_id];

			route_id_tmp = lv_result[tmp_i]->routeid;

			strcpy(ex_result[ex_ren].cigars, cigars);

			ex_result[ex_ren].seq_dis = tmp_node_p->seq_dis + lv_result[tmp_i]->offset;//seq_node_n[seq_id][tmp_node_p->seq_dis]

			c_result[c_ren].seqid = seq_id;
			c_result[c_ren].resultid = ex_ren;
			c_ren++;
			ex_ren++;
		}
	}
	(*re_c_result_n) = c_ren;

	if(route_arr)	free(route_arr);

	if(route_arrn)	free(route_arrn);

	for(tmp_i = 0; tmp_i < node_array_n; tmp_i++)
		if(node_array[tmp_i])	free(node_array[tmp_i]);

	return alt_info_cnt - 1;
}

char* combine_cigars(char* cigarBuf1, char* cigarBuf2, int16_t s_r_o_l, int16_t s_r_o_r, uint16_t read_len, uint8_t tid, uint8_t l_j, uint8_t r_j)
{
	uint16_t f_cigarn = 0, s_o = 0, f_i = 0, f_c = 0, pchl = 0, sn = 0, snt = 0;
	int m_m_n = 0;

	char* pch = NULL;
	char* saveptr = NULL;
	char* cigar_p = NULL;
	char* f_cigar = NULL;
	char* b_cigar = NULL;
	char* str_o = NULL;

	cigar_p = cigar_ps[tid];
	f_cigar = f_cigars[tid];
	b_cigar = b_cigars[tid];
	str_o = str_os[tid];

	m_m_n = s_r_o_r - s_r_o_l - 1;

	if(l_j)
	{
		f_cigarn = strlen(cigarBuf1);

		strcpy(str_o, cigarBuf1);

		pch = strtok_r(cigarBuf1,"DMIX", &saveptr);

		while (pch != NULL)
		{
			pchl = strlen(pch);
			f_cigar[f_cigarn - f_c - 2] = atoi(pch);
			s_o += (pchl + 1);
			f_cigar[f_cigarn - f_c - 1] = str_o[s_o - 1];

			f_c += 2;

			pch = strtok_r(NULL, "DMIX", &saveptr);
		}
	}

	if(r_j)
	{
		strcpy(b_cigar, cigarBuf2);
		pch = strtok(cigarBuf2,"DMIX");

		if(pch)	pchl = strlen(pch);
	}

	if(l_j && r_j)
	{
		if((f_cigar[f_cigarn - 1] == 'M') && (b_cigar[pchl] == 'M'))
		{
			f_cigar[f_cigarn - 2] += (m_m_n + atoi(pch));

			for(f_i = 0; f_i < f_c; f_i += 2)
			{
				sn = sprintf(cigar_p + snt, "%u%c", f_cigar[f_cigarn - f_c + f_i], f_cigar[f_cigarn + 1 - f_c + f_i]);
				snt += sn;
			}
			sn = sprintf(cigar_p + snt, "%s", b_cigar + pchl + 1);
			snt += sn;
		}
		else if(f_cigar[f_cigarn - 1] == 'M')
		{
			f_cigar[f_cigarn - 2] += m_m_n;
			for(f_i = 0; f_i < f_c; f_i += 2)
			{
				sn = sprintf(cigar_p + snt, "%u%c", f_cigar[f_cigarn - f_c + f_i], f_cigar[f_cigarn + 1 - f_c + f_i]);
				snt += sn;
			}
			sn = sprintf(cigar_p + snt, "%s",b_cigar);
			snt += sn;
		}
		else if(b_cigar[pchl] == 'M')
		{
			for(f_i = 0; f_i < f_c; f_i += 2)
			{
				sn = sprintf(cigar_p + snt, "%u%c", f_cigar[f_cigarn - f_c + f_i], f_cigar[f_cigarn + 1 - f_c + f_i]);
				snt += sn;
			}

			sn = sprintf(cigar_p + snt, "%uM%s", m_m_n + atoi(pch), b_cigar + pchl + 1);
			snt += sn;
		}
		else
		{
			for(f_i = 0; f_i < f_c; f_i += 2)
			{
				sn = sprintf(cigar_p + snt, "%u%c", f_cigar[f_cigarn - f_c + f_i], f_cigar[f_cigarn + 1 - f_c + f_i]);
				snt += sn;
			}
			sn = sprintf(cigar_p + snt, "%uM%s", m_m_n, b_cigar);
			snt += sn;
		}
	}
	else if(l_j)
	{
		if(f_cigar[f_cigarn - 1] == 'M')
		{
			f_cigar[f_cigarn - 2] += m_m_n;
			for(f_i = 0; f_i < f_c; f_i += 2)
			{
				sn = sprintf(cigar_p + snt, "%u%c", f_cigar[f_cigarn - f_c + f_i], f_cigar[f_cigarn + 1 - f_c + f_i]);
				snt += sn;
			}
		}
		else
		{
			for(f_i = 0; f_i < f_c; f_i += 2)
			{
				sn = sprintf(cigar_p + snt, "%u%c", f_cigar[f_cigarn - f_c + f_i], f_cigar[f_cigarn + 1 - f_c + f_i]);
				snt += sn;
			}
			sn = sprintf(cigar_p + snt, "%uM", m_m_n);
			snt += sn;
		}

	}
	else if(r_j)
	{
		if(b_cigar[pchl] == 'M')
		{
			sn = sprintf(cigar_p + snt, "%uM%s", m_m_n + atoi(pch), b_cigar + pchl + 1);
			snt += sn;
		}
		else
		{
			sn = sprintf(cigar_p + snt, "%uM%s", m_m_n, b_cigar);
			snt += sn;
		}

	}
	else
	{
		sn = sprintf(cigar_p + snt, "%uM", m_m_n);
		snt += sn;
	}

	cigar_p[snt] = '\0';

	return cigar_p;
}

int compare_alt_sort(const void * a, const void * b)
{
	alt_info_a* alt1 = (alt_info_a* )a;
	alt_info_a* alt2 = (alt_info_a* )b;

	uint16_t nodeid1 = alt1->nodeid;
	uint16_t nodeid2 = alt2->nodeid;
	uint16_t node_offset1 = alt1->node_offset;
	uint16_t node_offset2 = alt2->node_offset;

	if(nodeid1 > nodeid2)
		return 1;
	else if(nodeid1 < nodeid2)
		return -1;
	else
	{
		if(node_offset1 > node_offset2)
			return 1;
		else if(node_offset1 < node_offset2)
			return -1;
		else	return 0;
	}
}
int compare_com_result(const void * a, const void * b)
{
	com_result* com_result1 = (com_result* )a;
	com_result* com_result2 = (com_result* )b;

	uint16_t seqid1 = com_result1->seqid;
	uint16_t seqid2 = com_result1->seqid;

	if(seqid1 > seqid2)
		return 1;
	else if(seqid1 < seqid2)
		return -1;
	else	return 0;
}
int binsearch_com(uint16_t x, com_result* v, int n)
{
	int low, high, mid;

	low = 0;
	high = n - 1;

	while ( low <= high )
	{
		mid = ((low + high) >> 1);
		if(x < v[mid].seqid)
		{
			high = mid - 1;
		}
		else if(x > v[mid].seqid)
		{
			low = mid + 1;
		}
		else     /*found match*/
		{
			return mid;
		}
	}

	return -1;
}

void create_route_array(uint16_t cnt_pre, uint16_t cnt, uint16_t** arr_a, uint32_t** arr_n_a, uint32_t* qnum_a, uint32_t** route_index_a, uint32_t* route_node, alt_info_a* alt_arr)
{
	uint8_t ty = 0;
	uint8_t index_f = 0;

	uint16_t i = 0;
	uint16_t n = 0;
	uint16_t ep = 0;
	uint16_t c_l = 0;

	uint32_t q_s = 0, q_e = 0, q_n = 0;
	uint32_t qhead = 0, qrear = 0;
	uint32_t arrn_tmp = 0;

	uint16_t* arr = NULL;

	uint32_t* route_index = NULL;
	uint32_t* arr_n = NULL;

	arr = (*arr_a);
	arr_n = (*arr_n_a);

	qhead = (*qnum_a);
	qrear = qhead;

	//debug
	uint16_t qrear_pre = qrear;

	route_index = (*route_index_a);

	n = cnt - cnt_pre;

	for(i = 0; i < n; i++, qrear++)
	{
		arrn_tmp = arr_n[qrear - 1];
		arr[arrn_tmp] = cnt_pre + i;
		arr_n[qrear] = arrn_tmp + 1;
		route_node[cnt_pre + i] = qrear;
	}

	//debug
	arr_tol += n;

	while(qhead != qrear)
	{
		//get the node
		q_s = arr_n[qhead - 1];
		q_e = arr_n[qhead];

		c_l = arr[q_e - 1];
		q_n = q_e - q_s;

		ep = alt_arr[c_l].alt_end;
		ty = alt_arr[c_l].alt_type;

		index_f = 0;

		for(i = c_l + 1; i < cnt; i++)
		{
			if(alt_arr[i].alt_start >= ep + ty)
			{
				arrn_tmp = arr_n[qrear - 1];
				memcpy(arr + arrn_tmp, arr + q_s, q_n << 1);
				arr[arrn_tmp + q_n] = i;

				//debug
				arr_tol += (q_n + 1);

				arr_n[qrear] = arrn_tmp + q_n + 1;

				if(!index_f)
				{
					route_index[qhead] = qrear;
					index_f = 1;
				}

				qrear++;
			}
		}

		qhead++;
	}

	(*qnum_a) = qhead;
}
