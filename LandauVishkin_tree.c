
#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <assert.h>

//#include "extension_tree_alt.h"
#include "LandauVishkin_tree.h"
//debug
#include "load_input.h"


#ifdef	ANCHOR_HASH_ALI

uint8_t k_anchor = 10;
uint8_t anchor_seed_d = 5;
uint8_t k_anchor_back = 4;
uint8_t anchor_back_mask = 0Xf;
uint8_t k_anchor_re = 22;
uint8_t k_anchor_b = 20;

uint16_t*** anchor_hash = NULL;
uint8_t*** anchor_array = NULL;
uint16_t*** anchor_point = NULL;
uint16_t*** anchor_pos = NULL;

anchor_seed** anchor_seed_buffer = NULL;

#endif

#ifdef	REDUCE_ANCHOR
uint8_t anchor_seed_length_thr = 10;//11: 1993234; 12 1992946; 15 1990***; 20 1989***; 30 1988***
#else
uint8_t anchor_seed_length_thr = 0;
#endif

uint32_t insert_dis = 500;
int64_t devi = 50 * 3;

read_h** rh = NULL;

#ifdef	ANCHOR_HASH_ALI

int compare_read_hash(const void * a, const void * b)
{
    read_h* rh1 = (read_h* )a;
    read_h* rh2 = (read_h* )b;

    if(rh1->des < rh2->des)
        return -1;
    else if(rh1->des > rh2->des)
        return 1;
    else    return 0;

}

int comepare_anchor_seed(const void * a, const void * b)
{
	anchor_seed* as1 = (anchor_seed* )a;
	anchor_seed* as2 = (anchor_seed* )b;
	
	if(as1->seed_length > as2->seed_length)
		return -1;
	else if(as1->seed_length < as2->seed_length)
		return 1;
	else	return 0;
}
#endif

//patternLen must be <= textLen
int computeEditDistance_tree(
    treenode** node_array,
    uint8_t* pattern,
    int patternLen,
    int k,
    ext_queue** lv_queue,
    ext_result*** lv_result,
    uint16_t* re_cnt,
    alt_info_a* alt_info_array,
    uint8_t* aseq_array,
    uint8_t* seqs,
    uint32_t qnum,//
	uint32_t* route_id_mark, 
	uint32_t* route_id_mark_pre,
	route_node** tmp_rnode_array,
	uint32_t** routeid_array,
	uint16_t* route_arr,
	uint32_t* route_arrn,
	uint32_t* route_index,
	uint16_t node_array_n,
	uint32_t* route_nodes,
	uint16_t* route_sep)
{
	uint8_t mis_flag = 0;
	uint8_t tmp_alt_flag = 0;
	uint8_t test_bp = 0;
	uint8_t tmp_alt_type = 0;
	uint8_t alt_search_flag = 0;
	uint8_t route_search_flag = 0;

	uint16_t tmp_pdis = 0;
	uint16_t tmp_tdis = 0;
	uint16_t pdis_i = 0;
	uint16_t re_n = 0;//max number of result
	uint16_t next_nodes_i = 0;
	uint16_t bucket_char_n_i = 0;
	uint16_t tmp_alt_info = 0;
	uint16_t tmp_alt_info_c = 0;
	uint16_t tmp_node_id_alt = 0;
	uint16_t del_p_i = 0;
	uint16_t del_n = 0;
	uint16_t tmp_del_p = 0;
	uint16_t tmp_node_i = 0;
	uint16_t tmp_offset = 0;
	uint16_t tmp_alt_length = 0;
	uint16_t ori_index_tmp = 0;
	uint16_t seq_index_tmp = 0;
	uint16_t tra_n = 0;
	uint16_t r_n = 0;
	uint16_t r_t = 0;
	
	uint32_t rnode_cnt = 0;
	uint32_t r_s = 0;
	uint32_t c_id = 0;
	uint32_t tmp_i = 0;
	uint32_t r_i = 0;
	uint32_t rm_i = 0;
	uint32_t rm_i_pre = 0;
	//uint32_t tmp_i_max_pre = 0;
	
	//to modify
	//uint32_t tmp_route = 0;
	uint32_t tmp_route_id = 0;
	uint32_t tmp_alt_seq_begin = 0;
	uint32_t qhead = 0;
	uint32_t qrear = 0;
	uint32_t route_id = 0;

	uint32_t* route_max_array = NULL;
	
	treenode* tmp_node = NULL;
	route_node* r_node = NULL;

	if((route_id_mark == NULL) || (route_id_mark_pre == NULL))
	{
		fprintf(stderr, "Fail to alloc\n");
		exit(1);
	}
	
	route_max_array = calloc((k + 1) * (k + 1), 4);
	
	rnode_cnt = 1;
	qrear = 1;
	//route_id = 1;
	
	//node 0 enter into queue
	lv_queue[0]->nodeid = 0;
	lv_queue[0]->offset = 0;
	lv_queue[0]->pdis = 0;
	lv_queue[0]->del_p = 0;
	lv_queue[0]->route_id = 0;
	lv_queue[0]->alt_flag = 0;//
	lv_queue[0]->alt_info_cnt = node_array[0]->alt_info_cnt;
	lv_queue[0]->alt_type = 0;

	while(qhead != qrear)
	{
		tmp_node_i = lv_queue[qhead]->nodeid;
		tmp_offset = lv_queue[qhead]->offset;
		tmp_alt_flag = lv_queue[qhead]->alt_flag;
		tmp_route_id = lv_queue[qhead]->route_id;
		tmp_pdis = lv_queue[qhead]->pdis;
		tmp_del_p = lv_queue[qhead]->del_p;	
		
		tmp_node = node_array[tmp_node_i];

		mis_flag = 0;

		if(tmp_alt_flag)//alt branch
		{
			tmp_alt_info_c = lv_queue[qhead]->alt_info_cnt;
			tmp_alt_seq_begin = alt_info_array[tmp_alt_info_c].alt_seq_begin;
			tmp_alt_length = alt_info_array[tmp_alt_info_c].alt_length;

			for(bucket_char_n_i = tmp_offset, pdis_i = tmp_pdis, del_p_i = tmp_del_p; (bucket_char_n_i < tmp_alt_length) && (pdis_i < patternLen); bucket_char_n_i++, del_p_i++)
			{
				if(del_p_i < del_n)	continue;

				//test_bp = ((ref_seqv_array[(bucket_char_n_i + tmp_alt_seq_begin) >> 4] >> ((15 - ((bucket_char_n_i + tmp_alt_seq_begin) & 0Xf)) << 2)) & 0Xf);
				test_bp = aseq_array[tmp_alt_seq_begin + bucket_char_n_i];

				if(test_bp != pattern[pdis_i])
				{
					r_node = tmp_rnode_array[rnode_cnt];
					r_node->nodeid = tmp_node_i;
					r_node->offset = bucket_char_n_i;
					r_node->pdis = pdis_i;
					//r_node->tdis = bucket_char_n_i - tmp_offset + tmp_tdis;
					r_node->alt_index = tmp_alt_info_c;
					r_node->alt_flag = tmp_alt_flag;
					r_node->route_id = tmp_route_id;
					r_node->pre_index = 0;
					r_node->cigar = 4;
					//r_node->del_p = 0;

					routeid_array[0][tmp_route_id] = rnode_cnt;

					route_id_mark[rm_i] = tmp_route_id;
					rm_i++;
						
					if(tmp_route_id > route_max_array[0])	route_max_array[0] = tmp_route_id;

					rnode_cnt++;

					mis_flag = 1;
					break;
				}
				pdis_i++;
			}

			if(pdis_i == patternLen)
			{
				if(re_n < MAX_EXT_ALT_POS)
				{
					((*lv_result)[re_n])->rnode_index = 0;
				((*lv_result)[re_n])->routeid = tmp_route_id;
				((*lv_result)[re_n])->nodeid = tmp_node_i;
				((*lv_result)[re_n])->edit = 0;
				((*lv_result)[re_n])->offset = tmp_offset;
				re_n++;
				}
				
			}
			else  //locate to tree
			{
				if(!mis_flag)
				{			
					tmp_node_id_alt = alt_info_array[tmp_alt_info_c].ending_node_id;

					lv_queue[qrear]->nodeid = tmp_node_id_alt;//ending node id
					lv_queue[qrear]->offset = alt_info_array[tmp_alt_info_c].ending_node_offset;//ending node offset
					lv_queue[qrear]->pdis = pdis_i;
					//lv_queue[qrear]->tdis = bucket_char_n_i - tmp_offset + tmp_tdis;
					lv_queue[qrear]->alt_flag = 0;
					lv_queue[qrear]->alt_info_cnt = 0;
					lv_queue[qrear]->del_p = del_p_i;
					lv_queue[qrear]->alt_type = alt_info_array[tmp_alt_info_c].alt_type;
					lv_queue[qrear]->route_id = tmp_route_id;///route_id
					
					qrear = ((qrear + 1) % NODES_QUEUE_MAX);
					///route_id++;								
				}
			}
		}
		else  //not alt branch
		{
			if(tmp_node->alt_nums)//have alts
			{
				tmp_alt_info_c = tmp_node->alt_info_cnt;
				tmp_alt_type = lv_queue[qhead]->alt_type;

#ifdef	ALT_TR
				for(tmp_i = 0; (tmp_i < tmp_node->alt_info_n) && (alt_info_array[tmp_alt_info_c + tmp_i].node_offset < (tmp_offset + tmp_alt_type)); tmp_i++)
					;
				tmp_alt_info_c += tmp_i;
#endif
				if(tmp_i == tmp_node->alt_info_n)
					alt_search_flag = 0;
				else{
					alt_search_flag = 1;
					tmp_alt_info = alt_info_array[tmp_alt_info_c].node_offset;
				}

				for(bucket_char_n_i = tmp_offset, pdis_i = tmp_pdis, del_p_i = tmp_del_p; (bucket_char_n_i < tmp_node->seqs_length) && (pdis_i < patternLen); bucket_char_n_i++, del_p_i++)
				{
					if(del_p_i < del_n)	continue;
					/// + tmp_alt_type
					while(alt_search_flag && (tmp_alt_info == bucket_char_n_i) && (alt_info_array[tmp_alt_info_c].nodeid == tmp_node_i))//ALT
					{
						lv_queue[qrear]->nodeid = tmp_node_i;
						lv_queue[qrear]->offset = 0;//bucket_char_n_i
						lv_queue[qrear]->pdis = pdis_i;/// + tmp_alt_type
						//lv_queue[qrear]->tdis = bucket_char_n_i - tmp_offset + tmp_tdis;
						lv_queue[qrear]->alt_flag = 1;
						lv_queue[qrear]->alt_info_cnt = tmp_alt_info_c;
						lv_queue[qrear]->del_p = del_p_i;
						lv_queue[qrear]->alt_type = 0;
						
						//route
						ori_index_tmp = alt_info_array[tmp_alt_info_c].ori_index;
						
						if(tmp_route_id > node_array_n)
						{
							c_id = route_index[tmp_route_id - node_array_n];
							if(c_id)
							{
								seq_index_tmp = alt_info_array[tmp_alt_info_c].seq_index;
								tra_n = route_sep[seq_index_tmp + 1] - route_sep[seq_index_tmp];
								//tra_n -= route_arr[route_arrn[tmp_route_id - node_array_n] - 1];
								
								r_n = route_arrn[c_id] - route_arrn[c_id - 1];
								r_s = route_arrn[c_id - 1] + r_n - 1;
								for(tmp_i = r_s, r_t = 0, route_search_flag = 0; r_t < tra_n; tmp_i += r_n, r_t++)
									if(route_arr[tmp_i] == ori_index_tmp)
									{
										route_search_flag = 1;
										break;
									}
								if(route_search_flag)	route_id = c_id + r_t + node_array_n;		
								else	route_id = tmp_route_id;
							}else	route_id = tmp_route_id;
						}else{
							route_id = route_nodes[ori_index_tmp] + node_array_n;
						}
						//

						lv_queue[qrear]->route_id = route_id;

						qrear = ((qrear + 1) % NODES_QUEUE_MAX);
						//route_id++;

						++tmp_alt_info_c;
						tmp_alt_info = alt_info_array[tmp_alt_info_c].node_offset;
					}

					test_bp = seqs[tmp_node->seqs_start + bucket_char_n_i];

					if(!((test_bp >> pattern[pdis_i]) & 1))
					{
						r_node = tmp_rnode_array[rnode_cnt];
						r_node->nodeid = tmp_node_i;
						r_node->offset = bucket_char_n_i;
						r_node->pdis = pdis_i;
						//r_node->tdis = bucket_char_n_i - tmp_offset + tmp_tdis;
						r_node->alt_index = 0;
						r_node->alt_flag = tmp_alt_flag;
						r_node->route_id = tmp_route_id;
						r_node->pre_index = 0;
						r_node->cigar = 4;
						//r_node->del_p = 0;

						routeid_array[0][tmp_route_id] = rnode_cnt;

						route_id_mark[rm_i] = tmp_route_id;
						rm_i++;
				
						if(tmp_route_id > route_max_array[0])	route_max_array[0] = tmp_route_id;

						rnode_cnt++;

						mis_flag = 1;
						break;
					}
	
					pdis_i++;
				}

				if(pdis_i == patternLen)
				{
					if(re_n < MAX_EXT_ALT_POS)
					{
							((*lv_result)[re_n])->rnode_index = 0;
					((*lv_result)[re_n])->routeid = tmp_route_id;
					((*lv_result)[re_n])->nodeid = tmp_node_i;
					((*lv_result)[re_n])->edit = 0;
					((*lv_result)[re_n])->offset = bucket_char_n_i;
					re_n++;
					}
					
				}
				else  //branch
				{
					if(!mis_flag)
					{
						for(next_nodes_i = 0; next_nodes_i < tmp_node->next_nodes_n; next_nodes_i++)
						{
							tmp_node_id_alt = tmp_node->next_nodes[next_nodes_i];

							lv_queue[qrear]->nodeid = tmp_node_id_alt;
							lv_queue[qrear]->offset = 0;
							lv_queue[qrear]->pdis = pdis_i;
							//lv_queue[qrear]->tdis = bucket_char_n_i + tmp_tdis;
							lv_queue[qrear]->alt_flag = 0;
							lv_queue[qrear]->alt_info_cnt = 0;
							lv_queue[qrear]->del_p = del_p_i;
							lv_queue[qrear]->alt_type = 0;
							
							if(tmp_route_id > node_array_n)	lv_queue[qrear]->route_id = tmp_route_id;
							else	lv_queue[qrear]->route_id = tmp_node_id_alt;

							qrear = ((qrear + 1) % NODES_QUEUE_MAX);
							//route_id++;

						}
					}
				}
			}
			else  //not have alts
			{
				for(bucket_char_n_i = tmp_offset, pdis_i = tmp_pdis, del_p_i = tmp_del_p; (bucket_char_n_i < tmp_node->seqs_length) && (pdis_i < patternLen); bucket_char_n_i++, del_p_i++)
				{
					if(del_p_i < del_n)	continue;

					test_bp = seqs[tmp_node->seqs_start + bucket_char_n_i];

					if(!((test_bp >> pattern[pdis_i]) & 1))
					{
						r_node = tmp_rnode_array[rnode_cnt];
						r_node->nodeid = tmp_node_i;
						r_node->offset = bucket_char_n_i;
						r_node->pdis = pdis_i;
						//r_node->tdis = bucket_char_n_i - tmp_offset + tmp_tdis;
						r_node->alt_index = 0;
						r_node->alt_flag = tmp_alt_flag;
						r_node->route_id = tmp_route_id;
						r_node->pre_index = 0;
						r_node->cigar = 4;
						//r_node->del_p = 0;

						routeid_array[0][tmp_route_id] = rnode_cnt;

						route_id_mark[rm_i] = tmp_route_id;
						rm_i++;
					
						if(tmp_route_id > route_max_array[0])	route_max_array[0] = tmp_route_id;

						rnode_cnt++;

						mis_flag = 1;
						break;
					}
					pdis_i++;
				}

				if(pdis_i == patternLen)
				{
					if(re_n < MAX_EXT_ALT_POS)
					{
						((*lv_result)[re_n])->rnode_index = 0;
					((*lv_result)[re_n])->routeid = tmp_route_id;
					((*lv_result)[re_n])->nodeid = tmp_node_i;
					((*lv_result)[re_n])->edit = 0;
					((*lv_result)[re_n])->offset = bucket_char_n_i;
					re_n++;	
					}
				}
				else  //branch
				{
					if(!mis_flag)
					{
						for(next_nodes_i = 0; next_nodes_i < tmp_node->next_nodes_n; next_nodes_i++)
						{
							tmp_node_id_alt = tmp_node->next_nodes[next_nodes_i];

							lv_queue[qrear]->nodeid = tmp_node_id_alt;
							lv_queue[qrear]->offset = 0;
							lv_queue[qrear]->pdis = pdis_i;
							//lv_queue[qrear]->tdis = bucket_char_n_i - tmp_offset + tmp_tdis;
							lv_queue[qrear]->alt_flag = 0;
							lv_queue[qrear]->alt_info_cnt = 0;
							lv_queue[qrear]->del_p = del_p_i;
							lv_queue[qrear]->alt_type = 0;
							
							if(tmp_route_id > node_array_n)	lv_queue[qrear]->route_id = tmp_route_id;
							else	lv_queue[qrear]->route_id = tmp_node_id_alt;

							qrear = ((qrear + 1) % NODES_QUEUE_MAX);
							//route_id++;
						}
					}
				}
			}
		}

		qhead = ((qhead + 1) % NODES_QUEUE_MAX);
	}

	int e = 0, d = 0, d_tmp1 = 0, d_tmp2 = 0, e_tmp = 0, e_tmp_c = 0;
	int best = 0, left = 0, right = 0;
	int best_index = 0, left_index = 0, right_index = 0, tmp_index = 0;
	int16_t array_index_best = 0, array_index_left = 0, array_index_right = 0, array_index_c = 0;
	uint8_t route_type = 0;
	uint32_t route_i = 0, route_i_pre = 0;

	for ( e = 1; e <= k; e++)
	{
		memcpy(route_id_mark_pre, route_id_mark, rm_i << 2);
		qsort(route_id_mark_pre, rm_i, 4, compare_route_mask);
		rm_i_pre = rm_i;
		rm_i = 0;

		e_tmp = (e - 1) * (e - 1);
		e_tmp_c = e * e;

		// Search d's in the order 0, 1, -1, 2, -2, etc to find an alignment with as few indels as possible.
		for ( d = 0; d != e+1; d = (d > 0 ? -d : -d+1))
		{
			d_tmp1 = d - 1;
			d_tmp2 = d + 1;

			//trnasform tradational two-dimension L[][] matrix into linear array in order to reduce memory allocation
			array_index_best = (e_tmp + (d > 0 ? (d << 1) - 1 : -(d << 1)));// % max_route_ids
			array_index_left = (e_tmp + (d_tmp1 > 0 ? (d_tmp1 << 1) - 1 : -(d_tmp1 << 1)));// % max_route_ids
			array_index_right = (e_tmp + (d_tmp2 > 0 ? (d_tmp2 << 1) - 1 : -(d_tmp2 << 1)));// % max_route_ids

			array_index_c = (e_tmp_c + (d > 0 ? (d << 1) - 1 : -(d << 1)));// % max_route_ids

			route_i_pre = 0Xffffffff;
			for(r_i = 0; r_i < rm_i_pre; r_i++)
			{
				route_i = route_id_mark_pre[r_i];
				
				if(route_i == route_i_pre)	continue;
				
				route_i_pre = route_i;

				if((d < 1 - e) || (d > e - 1))	best = -2;
				else
				{
					best_index = routeid_array[array_index_best][route_i];
					if(best_index)	best = tmp_rnode_array[best_index]->pdis + 1;
					else	best = -2;
				}

				if((d_tmp1 < 1 - e) || (d_tmp1 > e - 1))	left = -2;
				else
				{
					left_index = routeid_array[array_index_left][route_i];
					if(left_index)	left = tmp_rnode_array[left_index]->pdis;
					else	left = -2;
				}

				if((d_tmp2 < 1 - e) || (d_tmp2 > e - 1))	right = -2;
				else
				{
					right_index = routeid_array[array_index_right][route_i];
					if(right_index)	right = tmp_rnode_array[right_index]->pdis + 1;
					else	right = -2;
				}
			
				if((best == -2) && (left == -2) && (right == -2))	continue;

				tmp_index = best_index;

				del_n = 1;
				route_type = 0;
				if (left > best)//D
				{
					route_type = 2;
					best = left;
					tmp_index = left_index;
				}
				if (right > best)//I
				{
					route_type = 1;
					best = right;
					tmp_index = right_index;
					del_n = 0;
				}
		
				r_node = tmp_rnode_array[tmp_index];
				
				if(best == patternLen)
				{
					if(re_n < MAX_EXT_ALT_POS)
					{
						((*lv_result)[re_n])->rnode_index = tmp_index;
					((*lv_result)[re_n])->routeid = route_i;
					((*lv_result)[re_n])->nodeid = r_node->nodeid;
					((*lv_result)[re_n])->edit = e;
					((*lv_result)[re_n])->offset = r_node->offset;
					((*lv_result)[re_n])->cigar = route_type;
					re_n++;	
					}
					
					continue;
				}

				qhead = 0;
				qrear = 1;

				//node 0 enter into queue
				lv_queue[0]->nodeid = r_node->nodeid;
				lv_queue[0]->offset = r_node->offset;
				lv_queue[0]->pdis = best;
				lv_queue[0]->alt_flag = r_node->alt_flag;
				lv_queue[0]->alt_info_cnt = r_node->alt_index;
				lv_queue[0]->route_id = r_node->route_id;

				//lv_queue[0]->del_p = r_node->del_p;
				lv_queue[0]->del_p = 0;
				//lv_queue[0]->alt_type = r_node->alt_type;
				lv_queue[0]->alt_type = 0;

				//qrear = ((qrear + 1) % NODES_QUEUE_MAX);

				while(qhead != qrear)
				{
					tmp_node_i = lv_queue[qhead]->nodeid;
			
					tmp_pdis = lv_queue[qhead]->pdis;
					//tmp_tdis = lv_queue[qhead]->tdis;
					tmp_offset = lv_queue[qhead]->offset;
					tmp_alt_flag = lv_queue[qhead]->alt_flag;
					tmp_route_id = lv_queue[qhead]->route_id;
					tmp_del_p = lv_queue[qhead]->del_p;

					tmp_node = node_array[tmp_node_i];
					
					mis_flag = 0;

					if(tmp_alt_flag)//alt branch
					{
						tmp_alt_info_c = lv_queue[qhead]->alt_info_cnt;
						tmp_alt_seq_begin = alt_info_array[tmp_alt_info_c].alt_seq_begin;
						tmp_alt_length = alt_info_array[tmp_alt_info_c].alt_length;

						for(bucket_char_n_i = tmp_offset, pdis_i = tmp_pdis, del_p_i = tmp_del_p; (bucket_char_n_i < tmp_alt_length) && (pdis_i < patternLen); bucket_char_n_i++, del_p_i++)
						{
							if(del_p_i < del_n)	continue;

							//test_bp = ((ref_seqv_array[(bucket_char_n_i + tmp_alt_seq_begin) >> 4] >> ((15 - ((bucket_char_n_i + tmp_alt_seq_begin) & 0Xf)) << 2)) & 0Xf);
							test_bp = aseq_array[tmp_alt_seq_begin + bucket_char_n_i];

							if(test_bp != pattern[pdis_i])
							{
								r_node = tmp_rnode_array[rnode_cnt];
								r_node->nodeid = tmp_node_i;
								r_node->offset = bucket_char_n_i;
								r_node->pdis = pdis_i;
								//r_node->tdis = bucket_char_n_i - tmp_offset + tmp_tdis;
								r_node->alt_index = tmp_alt_info_c;
								r_node->alt_flag = tmp_alt_flag;
								r_node->route_id = tmp_route_id;
								r_node->pre_index = tmp_index;
								r_node->cigar = route_type;
								//r_node->del_p = tmp_del_p;

								routeid_array[array_index_c][tmp_route_id] = rnode_cnt;

								route_id_mark[rm_i] = tmp_route_id;
								rm_i++;

								if(tmp_route_id > route_max_array[array_index_c])	route_max_array[array_index_c] = tmp_route_id;

								rnode_cnt++;

								mis_flag = 1;
								break;
							}
							pdis_i++;
						}

						if(pdis_i == patternLen)
						{
							if(re_n < MAX_EXT_ALT_POS)
							{
							((*lv_result)[re_n])->rnode_index = tmp_index;
							((*lv_result)[re_n])->routeid = tmp_route_id;
							((*lv_result)[re_n])->nodeid = tmp_node_i;
							((*lv_result)[re_n])->edit = e;
							((*lv_result)[re_n])->offset = tmp_offset;
							((*lv_result)[re_n])->cigar = route_type;
							re_n++;
							}
							
						}
						else  //locate to tree
						{
							if(!mis_flag)
							{
								tmp_node_id_alt = alt_info_array[tmp_alt_info_c].ending_node_id;

								lv_queue[qrear]->nodeid = tmp_node_id_alt;//ending node id
								lv_queue[qrear]->offset = alt_info_array[tmp_alt_info_c].ending_node_offset;//ending node offset
								lv_queue[qrear]->pdis = pdis_i;
								//lv_queue[qrear]->tdis = bucket_char_n_i - tmp_offset + tmp_tdis;
								lv_queue[qrear]->alt_flag = 0;
								lv_queue[qrear]->alt_info_cnt = 0;
								lv_queue[qrear]->del_p = del_p_i;
								lv_queue[qrear]->alt_type = alt_info_array[tmp_alt_info_c].alt_type;
								lv_queue[qrear]->route_id = tmp_route_id;///route_id
								
								qrear = ((qrear + 1) % NODES_QUEUE_MAX);
								///route_id++;

							}
						}
					}
					else  //not alt branch
					{
						if(tmp_node->alt_nums)//have alts
						{
							tmp_alt_info_c = tmp_node->alt_info_cnt;
							tmp_alt_type = lv_queue[qhead]->alt_type;

#ifdef	ALT_TR
							for(tmp_i = 0; (tmp_i < tmp_node->alt_info_n) && (alt_info_array[tmp_alt_info_c + tmp_i].node_offset < (tmp_offset + tmp_alt_type)); tmp_i++)
								;
							tmp_alt_info_c += tmp_i;
#endif
							if(tmp_i == tmp_node->alt_info_n)
								alt_search_flag = 0;
							else{
								alt_search_flag = 1;
								tmp_alt_info = alt_info_array[tmp_alt_info_c].node_offset;
							}	

							for(bucket_char_n_i = tmp_offset, pdis_i = tmp_pdis, del_p_i = tmp_del_p; (bucket_char_n_i < tmp_node->seqs_length) && (pdis_i < patternLen); bucket_char_n_i++, del_p_i++)
							{
								if(del_p_i < del_n)	continue;
								/// + tmp_alt_type
								while(alt_search_flag && (tmp_alt_info == bucket_char_n_i) && (alt_info_array[tmp_alt_info_c].nodeid == tmp_node_i))//ALT
								{
									lv_queue[qrear]->nodeid = tmp_node_i;
									lv_queue[qrear]->offset = 0;//bucket_char_n_i
									lv_queue[qrear]->pdis = pdis_i;/// + tmp_alt_type
									//lv_queue[qrear]->tdis = bucket_char_n_i - tmp_offset + tmp_tdis;
									lv_queue[qrear]->alt_flag = 1;
									lv_queue[qrear]->alt_info_cnt = tmp_alt_info_c;
									lv_queue[qrear]->del_p = del_p_i;
									lv_queue[qrear]->alt_type = 0;

									//route
									ori_index_tmp = alt_info_array[tmp_alt_info_c].ori_index;//2
									if(tmp_route_id > node_array_n)
									{
										c_id = route_index[tmp_route_id - node_array_n];//155 - 149 = 6
										if(c_id)
										{
											seq_index_tmp = alt_info_array[tmp_alt_info_c].seq_index;//seq ID 
											tra_n = route_sep[seq_index_tmp + 1] - route_sep[seq_index_tmp];//how many routes on this seq 27 - 25 = 2
											//tra_n -= route_arr[route_arrn[tmp_route_id - node_array_n] - 1];//461 - 436 = 25
								
											r_n = route_arrn[c_id] - route_arrn[c_id - 1];
											r_s = route_arrn[c_id - 1] + r_n - 1;
											for(tmp_i = r_s, r_t = 0, route_search_flag = 0; r_t < tra_n; tmp_i += r_n, r_t++)
												if(route_arr[tmp_i] == ori_index_tmp)
												{
													route_search_flag = 1;
													break;
												}
											if(route_search_flag)	route_id = c_id + r_t + node_array_n;
											else	route_id = tmp_route_id;
										}else	route_id = tmp_route_id;
									}else{
										route_id = route_nodes[ori_index_tmp] + node_array_n;
									}

									lv_queue[qrear]->route_id = route_id;

									qrear = ((qrear + 1) % NODES_QUEUE_MAX);
									//route_id++;

									++tmp_alt_info_c;
									tmp_alt_info = alt_info_array[tmp_alt_info_c].node_offset;
								}

								test_bp = seqs[tmp_node->seqs_start + bucket_char_n_i];

								//if(test_bp != pattern[pdis_i])
								if(!((test_bp >> pattern[pdis_i]) & 1))
								{
									r_node = tmp_rnode_array[rnode_cnt];
									r_node->nodeid = tmp_node_i;
									r_node->offset = bucket_char_n_i;
									r_node->pdis = pdis_i;
									//r_node->tdis = bucket_char_n_i - tmp_offset + tmp_tdis;
									r_node->alt_index = 0;
									r_node->alt_flag = tmp_alt_flag;
									r_node->route_id = tmp_route_id;
									r_node->pre_index = tmp_index;
									r_node->cigar = route_type;
									//r_node->del_p = tmp_del_p;

									routeid_array[array_index_c][tmp_route_id] = rnode_cnt;

									route_id_mark[rm_i] = tmp_route_id;
									rm_i++;
									
									if(tmp_route_id > route_max_array[array_index_c])	route_max_array[array_index_c] = tmp_route_id;

									rnode_cnt++;

									mis_flag = 1;
									break;
								}
							
								pdis_i++;
							}
							if(pdis_i == patternLen)
							{
								if(re_n < MAX_EXT_ALT_POS)
								{
								((*lv_result)[re_n])->rnode_index = tmp_index;
								((*lv_result)[re_n])->routeid = tmp_route_id;
								((*lv_result)[re_n])->nodeid = tmp_node_i;
								((*lv_result)[re_n])->edit = e;
								((*lv_result)[re_n])->offset = bucket_char_n_i;
								((*lv_result)[re_n])->cigar = route_type;
								re_n++;
								}
								
							}
							else  //branch
							{
								if(!mis_flag)
								{
									for(next_nodes_i = 0; next_nodes_i < tmp_node->next_nodes_n; next_nodes_i++)
									{
										tmp_node_id_alt = tmp_node->next_nodes[next_nodes_i];

										lv_queue[qrear]->nodeid = tmp_node_id_alt;
										lv_queue[qrear]->offset = 0;
										lv_queue[qrear]->pdis = pdis_i;
										//lv_queue[qrear]->tdis = bucket_char_n_i + tmp_tdis;
										lv_queue[qrear]->alt_flag = 0;
										lv_queue[qrear]->alt_info_cnt = 0;
										lv_queue[qrear]->del_p = del_p_i;
										lv_queue[qrear]->alt_type = 0;
										
										if(tmp_route_id > node_array_n)	lv_queue[qrear]->route_id = tmp_route_id;
										else	lv_queue[qrear]->route_id = tmp_node_id_alt;

										qrear = ((qrear + 1) % NODES_QUEUE_MAX);
										//route_id++;

									}
								}
							}
						}
						else  //not have alts
						{
							for(bucket_char_n_i = tmp_offset, pdis_i = tmp_pdis, del_p_i = tmp_del_p; (bucket_char_n_i < tmp_node->seqs_length) && (pdis_i < patternLen); bucket_char_n_i++, del_p_i++)
							{
								if(del_p_i < del_n)	continue;

								//
								test_bp = seqs[tmp_node->seqs_start + bucket_char_n_i];

								//if(test_bp != pattern[pdis_i])
								if(!((test_bp >> pattern[pdis_i]) & 1))
								{
									r_node = tmp_rnode_array[rnode_cnt];
									r_node->nodeid = tmp_node_i;
									r_node->offset = bucket_char_n_i;
									r_node->pdis = pdis_i;
									//r_node->tdis = bucket_char_n_i - tmp_offset + tmp_tdis;
									r_node->alt_index = 0;
									r_node->alt_flag = tmp_alt_flag;
									r_node->route_id = tmp_route_id;
									r_node->pre_index = tmp_index;
									r_node->cigar = route_type;
									//r_node->del_p = tmp_del_p;

									routeid_array[array_index_c][tmp_route_id] = rnode_cnt;

									route_id_mark[rm_i] = tmp_route_id;
									rm_i++;
									
									if(tmp_route_id > route_max_array[array_index_c])	route_max_array[array_index_c] = tmp_route_id;

									rnode_cnt++;

									mis_flag = 1;
									break;
								}
								pdis_i++;
							}

							if(pdis_i == patternLen)
							{
								if(re_n < MAX_EXT_ALT_POS)
								{
								((*lv_result)[re_n])->rnode_index = tmp_index;
								((*lv_result)[re_n])->routeid = tmp_route_id;
								((*lv_result)[re_n])->nodeid = tmp_node_i;
								((*lv_result)[re_n])->edit = e;
								((*lv_result)[re_n])->offset = bucket_char_n_i;
								((*lv_result)[re_n])->cigar = route_type;
								re_n++;
								}	
							}
							else  //branch
							{
								if(!mis_flag)
								{
									for(next_nodes_i = 0; next_nodes_i < tmp_node->next_nodes_n; next_nodes_i++)
									{
										tmp_node_id_alt = tmp_node->next_nodes[next_nodes_i];

										lv_queue[qrear]->nodeid = tmp_node_id_alt;
										lv_queue[qrear]->offset = 0;
										lv_queue[qrear]->pdis = pdis_i;
										//lv_queue[qrear]->tdis = bucket_char_n_i - tmp_offset + tmp_tdis;
										lv_queue[qrear]->alt_flag = 0;
										lv_queue[qrear]->alt_info_cnt = 0;
										lv_queue[qrear]->del_p = del_p_i;
										lv_queue[qrear]->alt_type = 0;
										
										if(tmp_route_id > node_array_n)	lv_queue[qrear]->route_id = tmp_route_id;
										else	lv_queue[qrear]->route_id = tmp_node_id_alt;
										

										qrear = ((qrear + 1) % NODES_QUEUE_MAX);
										//route_id++;
									}
								}
							}
						}
					}
					
					qhead = ((qhead + 1) % NODES_QUEUE_MAX);
				}
			}
		}
	}

	for(e = 0; e < (k + 1) * (k + 1); e++)
	{
		memset(routeid_array[e], 0, (route_max_array[e] + 1) << 2);// 300000
	}

	if(route_max_array)	free(route_max_array);

	(*re_cnt) = re_n;

	return -1;
}

int compare_route_mask(const void * a, const void * b)
{
	uint32_t pos1 = (*((uint32_t* )a));
	uint32_t pos2 = (*((uint32_t* )b));

		
	if(pos1 > pos2)
		return 1;
	else if(pos1 < pos2)
		return -1;
	else	return 0;
}
	