
#include <stdlib.h>
#include <unistd.h>
#include <malloc.h>
#include <string.h>
#include <dirent.h>
#include <inttypes.h>

#include "index_build.h"
#include "bit_operation.h"
#include "load_input.h"

#define	VCF_ADD
#define	VCF_ADD_MAX	50
#define	VCF_ADD_MIN	2

#define	INS_SEQ
#define	ALT_SEQ_MAX	20000

typedef struct vcf_record_sort
{
	uint32_t vcf_buff_sort[4];
} vcf_record_s;


uint8_t bit_shift_value = 0;
uint8_t vcf_bit_shift = 0;//cannot be > 7
uint16_t chr_no_n = 0;
uint32_t vcf_seqs_offset = 0;
uint32_t chr_posend_n = 0;
uint32_t chr_file_n_vcf = 0;


char** chr_no_a = NULL;
uint32_t* chr_vcf_n = NULL;
uint32_t* chr_vcf_seq_n = NULL;


const char* vcf_seq = "vcf_seq";
const char* vcf_pos_s = "vcf_pos";


int index_build(int argc, char *argv[])
{

	if(load_input_index(argc, argv) == 1)	return 1;

	printf("kmer size to build index: %u\n", k_t);
	fflush(stdout);
	//load_reffile_kmer();
#ifndef	DEBUG_VCF
	load_reffile_kmer_fa();

	file_kmer_qsort();
#else
	/*
	chr_end_n[0] = 2049;
	chr_end_n[1] = 48131944;
	chr_end_n[2] = 99436510;
	chr_end_n[3] = 254707070;

	chr_posend_n = 4;

	ref_seq_n = 254707069;

	strcpy(chr_names[0], "chr21");
	strcpy(chr_names[1], "chr22");
	strcpy(chr_names[2], "chrX");
	*/
	
	//chr_name_n = 25;
	chr_posend_n = 26;
	chr_file_n_vcf = chr_posend_n;
	ref_seq_n = 3095696031;

strcpy(chr_names[0],"chr1");
strcpy(chr_names[1],"chr2");
strcpy(chr_names[2],"chr3");
strcpy(chr_names[3],"chr4");
strcpy(chr_names[4],"chr5");
strcpy(chr_names[5],"chr6");
strcpy(chr_names[6],"chr7");
strcpy(chr_names[7],"chr8");
strcpy(chr_names[8],"chr9");
strcpy(chr_names[9],"chr10");
strcpy(chr_names[10],"chr11");
strcpy(chr_names[11],"chr12");
strcpy(chr_names[12],"chr13");
strcpy(chr_names[13],"chr14");
strcpy(chr_names[14],"chr15");
strcpy(chr_names[15],"chr16");
strcpy(chr_names[16],"chr17");
strcpy(chr_names[17],"chr18");
strcpy(chr_names[18],"chr19");
strcpy(chr_names[19],"chr20");
strcpy(chr_names[20],"chr21");
strcpy(chr_names[21],"chr22");
strcpy(chr_names[22],"chrX");
strcpy(chr_names[23],"chrY");
strcpy(chr_names[24],"chrM");
chr_end_n[0]=2049;
chr_end_n[1]=249252670;
chr_end_n[2]=492452043;
chr_end_n[3]=690474473;
chr_end_n[4]=881628749;
chr_end_n[5]=1062544009;
chr_end_n[6]=1233659076;
chr_end_n[7]=1392797739;
chr_end_n[8]=1539161761;
chr_end_n[9]=1680375192;
chr_end_n[10]=1815909939;
chr_end_n[11]=1950916455;
chr_end_n[12]=2084768350;
chr_end_n[13]=2199938228;
chr_end_n[14]=2307287768;
chr_end_n[15]=2409819160;
chr_end_n[16]=2500173913;
chr_end_n[17]=2581369123;
chr_end_n[18]=2659446371;
chr_end_n[19]=2718575354;
chr_end_n[20]=2781600874;
chr_end_n[21]=2829730769;
chr_end_n[22]=2881035335;
chr_end_n[23]=3036305895;
chr_end_n[24]=3095679461;
chr_end_n[25]=3095696032;

#endif

	if(tree_flag)	alt_index_build();


	return 0;
}

void alt_index_build()
{
	load_vcf_chr(vcf_file_name);

	chr_vcf_index();
}

void load_vcf_chr(char* vcf_file_name)
{
	FILE* fp_in_vcf = NULL;
	//FILE* fp_out_test = fopen("./test.vcf", "w");

	char file_d[ROUTE_LENGTH_MAX];
	char one_line[VCF_LINE_MAX + 1];
	char ref_num[VCF_REF_MAX];
	char alt_num[VCF_ALT_MAX];
	char vcf_ref[VCF_REF_MAX];
	char vcf_alt[VCF_ALT_MAX];
	char chr_no_tmp[40];
	char* pch = NULL;
	char* saveptr = NULL;
	char* pch_alt = NULL;
	char* saveptr_alt = NULL;

	uint8_t field_i = 0;
	uint8_t chr_name_f = 0;
	uint8_t skip_flag = 0;
	
	uint16_t ref_l = 0;//VCF_REF_MAX
	uint16_t alt_l = 0;//VCF_ALT_MAX
	uint16_t min_l = 0;
	uint16_t file_n = 0;//how many chrs at most?
	uint16_t chr_no = 0;
	uint16_t alt_field_i = 0;
	
	uint32_t vcf_line_n = 0;
	uint32_t vcf_pos = 0;
	uint32_t ite_i = 0;
	uint32_t vcf_n_tmp = 0;
	uint32_t vcf_buff[4];

	vcf_bit_shift = 2;
	file_n = 60;

#ifdef	VCF_ADD
	uint32_t chr_start = 0;
	uint32_t num_tmp = 0;
	uint32_t num_tmp_tmp = 0;
	uint64_t a_size = 0;
	uint64_t* buffer_ref_seq = NULL;
	FILE* fp_ref_seq = NULL;
	char char_tmp[VCF_ADD_MAX];//uint8_t

	ref_length_p = ref_seq_n;

	fp_ref_seq = fopen(ref_seq, "rb");
	if (fp_ref_seq == NULL)
	{
		fputs ("File error opening the seq file\n",stderr);
		exit (1);
	}
	fprintf(stderr, "ref_seq_n: %u\n", ref_seq_n);//3095696031

	a_size = (ref_seq_n >> 5) + 1;
	buffer_ref_seq = (uint64_t* )calloc(a_size, 8);

	num_tmp = fread(buffer_ref_seq, 8, a_size, fp_ref_seq);//63676768

	fprintf(stderr, "num_tmp %u %u vcf_rpos %s\n", num_tmp, ref_seq_n, vcf_rpos);//96740501 3095696031

	if(fp_ref_seq)	fclose(fp_ref_seq);
#endif

	fp_in_vcf = fopen(vcf_file_name, "r");
	if(fp_in_vcf == NULL)	fprintf(stderr, "Error of reading input vcf file\n");

	FILE** file_p = (FILE** )malloc(file_n * sizeof(FILE* ));
	if(file_p == NULL)	fprintf(stderr, "Error of allocating memory file_p\n");

	for(ite_i = 0; ite_i < file_n; ite_i++)
		file_p[ite_i] = NULL;

	FILE** file_p_pos = (FILE** )malloc(file_n * sizeof(FILE* ));
	if(file_p_pos == NULL)	fprintf(stderr, "Error of allocating memory file_p_pos\n");

	for(ite_i = 0; ite_i < file_n; ite_i++)
		file_p_pos[ite_i] = NULL;

	chr_vcf_n = (uint32_t* )calloc(file_n, 4);
	if(chr_vcf_n == NULL)	fprintf(stderr, "Error of allocating memory chr_vcf_n\n");

	chr_vcf_seq_n = (uint32_t* )calloc(file_n, 4);
	if(chr_vcf_seq_n == NULL)	fprintf(stderr, "Error of allocating memory chr_vcf_seq_n\n");

	chr_no_a = (char** )calloc(file_n, sizeof(char* ));
	if(chr_no_a == NULL)	fprintf(stderr, "Error of allocating memory chr_no_a\n");
	for(ite_i = 0; ite_i < file_n; ite_i++)
		if((chr_no_a[ite_i] = (char* )calloc(CHR_NAME_MAX, 1)) == NULL)
			fprintf(stderr, "Error of allocating memory chr_no_a\n");

	fprintf(stderr, "begin reading and loading vcf file\n");

	//debug
	FILE* fp_debug_vcf = fopen("./debug_vcf_record","w");

	while ((!feof(fp_in_vcf)) && (fgets(one_line, VCF_LINE_MAX + 2, fp_in_vcf) != NULL))
	{
		if(one_line[0] == '#')	continue;

		field_i = 0;
		pch = strtok_r(one_line,"\t", &saveptr);

		while (pch != NULL)
		{
			//fprintf("%s ", pch);

			if(field_i == 0)
			{
				//sscanf(pch, "%[1-9]", chr_no_tmp);
				//vcf chr only has the number
				strcpy(chr_no_tmp, pch);// + 3

				skip_flag = 0;
				if(!(strncmp(chr_no_tmp,"GL",2) && strncmp(chr_no_tmp,"MT",2) && strncmp(chr_no_tmp,"NC",2) && strncmp(chr_no_tmp,"hs",2)))
				{
					skip_flag = 1;
					break;
				}
				
				chr_name_f = 0;
				for (ite_i = 0; ite_i < chr_no_n; ite_i++)
					if(strcmp(chr_no_a[ite_i], chr_no_tmp) == 0)
					{
						chr_name_f = 1;
						break;
					}

				if(chr_name_f)
				{
					chr_no = ite_i;
				}
				else
				{
					chr_no = chr_no_n;
					strcpy(chr_no_a[chr_no_n++], chr_no_tmp);
				}
			}
			else if(field_i == 1)
			{
				vcf_pos = atoi(pch);
			}
			else if(field_i == 2)
			{

			}
			else if(field_i == 3)     //ref
			{
				strcpy(vcf_ref, pch);
			}
			else if(field_i == 4)     //alt
			{
				strcpy(vcf_alt, pch);
			}
			else
			{
				break;
			}

			pch = strtok_r(NULL, "\t", &saveptr);
			field_i++;
		}

		if(skip_flag)	continue;
		
		if(file_p[chr_no] == NULL)
		{
			//should change when put into deBGA, change route
			memset(file_d, 0, ROUTE_LENGTH_MAX);
			strcpy(file_d, index_route);
			strcat(file_d, vcf_seq);
			strcat(file_d, chr_no_tmp);

			file_p[chr_no] = fopen(file_d, "wb");
			if (file_p[chr_no] == NULL)
			{
				fprintf(stderr, "File error of creating vcf seq file %s %s %s\n",file_d, vcf_seq, chr_no_tmp);
				exit(1);
			}
		}
		if(file_p_pos[chr_no] == NULL)
		{
			//should change when put into deBGA
			memset(file_d, 0, ROUTE_LENGTH_MAX);
			strcpy(file_d, index_route);
			strcat(file_d, vcf_pos_s);
			strcat(file_d, chr_no_tmp);

			file_p_pos[chr_no] = fopen(file_d, "wb");
			if (file_p_pos[chr_no] == NULL)
			{
				fputs ("File error of creating vcf pos file\n",stderr);
				exit(1);
			}
		}

		//fprintf("%u %s %s\n", vcf_pos, vcf_ref, vcf_alt);

		vcf_seqs_offset = chr_vcf_seq_n[chr_no];
		ref_l = strlen(vcf_ref);
		pch_alt = strtok_r(vcf_alt,",", &saveptr_alt);
		
		while (pch_alt != NULL)
		{
			alt_l = strlen(pch_alt);
			for (ite_i = 0; ite_i < alt_l; ite_i++)
			{
				alt_num[ite_i] = nt_table_u[(uint8_t )pch_alt[ite_i]];
				if(alt_num[ite_i] == 4)	break;
			}
			if(ite_i < alt_l)
			{
				pch_alt = strtok_r(NULL, ",", &saveptr_alt);
				continue;
			}

			for (ite_i = 0; ite_i < ref_l; ite_i++) ref_num[ite_i] = nt_table_u[(uint8_t )vcf_ref[ite_i]];

			min_l = (ref_l > alt_l) ? alt_l:ref_l;

			//fprintf("L: %u %u\n", ref_l, alt_l);

			vcf_n_tmp = 0;
			
			if((ref_l == 1) && (alt_l == 1))//SNP
			{
				vcf_buff[0] = vcf_pos - 1;//vcf 1-offset
				vcf_buff[1] = 0;
				vcf_buff[1] |= (1 << vcf_bit_shift);
				vcf_buff[2] = vcf_seqs_offset;
				vcf_buff[3] = vcf_pos;
				fwrite(vcf_buff, 4, 4, file_p_pos[chr_no]);
				fwrite(&alt_num[0], 1, 1, file_p[chr_no]);
				vcf_seqs_offset++;

				vcf_n_tmp++;
				//fprintf("one\n");
			}
			else if((ref_l < alt_l) && (memcmp(ref_num, alt_num, min_l) == 0))//I strncmp
			{
				vcf_buff[0] = vcf_pos - 1 + min_l;
				vcf_buff[1] = 1;
				vcf_buff[1] |= ((alt_l - ref_l) << vcf_bit_shift);
				vcf_buff[2] = vcf_seqs_offset;
				vcf_buff[3] = vcf_buff[0];
				fwrite(vcf_buff, 4, 4, file_p_pos[chr_no]);
				fwrite(alt_num + ref_l, 1, alt_l - ref_l, file_p[chr_no]);
				vcf_seqs_offset += (alt_l - ref_l);

				vcf_n_tmp++;
				
				fprintf(fp_debug_vcf, "I: %u\n", vcf_pos);
				//fprintf("two\n");
#ifdef	VCF_ADD				
				if(alt_l - ref_l > VCF_ADD_MIN)
				{
					for (ite_i = 0; ite_i < chr_posend_n; ite_i++)
						if(strcmp(chr_no_tmp, chr_names[ite_i] + 3) == 0)
							break;

					chr_start = chr_end_n[ite_i] - 1;
					
					for (ite_i = 0; ite_i < alt_l - ref_l; ite_i++)
					{
						num_tmp = chr_start + vcf_pos + ite_i;
						if(alt_num[ref_l + ite_i] != ((buffer_ref_seq[num_tmp >> 5] >> ((31 - (num_tmp & 0X1f)) << 1)) & 0X3))
							break;
					}
					
					fprintf(stderr, "ite_i: %u chr_names: %s pos: %u ref_l: %u alt_l: %u %u %u\n", ite_i, chr_no_tmp, vcf_pos, ref_l, alt_l, ((buffer_ref_seq[(chr_start + vcf_pos) >> 5] >> ((31 - ((chr_start + vcf_pos) & 0X1f)) << 1)) & 0X3), ((buffer_ref_seq[(chr_start + vcf_pos + 1) >> 5] >> ((31 - ((chr_start + vcf_pos + 1) & 0X1f)) << 1)) & 0X3));
										
					if(ite_i)
					{
						fprintf(stderr, "VCF_ADD INS A\n");

						vcf_buff[0] = vcf_buff[0] + ite_i;//vcf_pos - 1 + min_l
						vcf_buff[1] = 1;
						vcf_buff[1] |= ((alt_l - ref_l) << vcf_bit_shift);
						vcf_buff[2] = vcf_seqs_offset;
						vcf_buff[3] = vcf_buff[0];
						fwrite(vcf_buff, 4, 4, file_p_pos[chr_no]);
						fwrite(alt_num + ref_l, 1, alt_l - ref_l, file_p[chr_no]);
						vcf_seqs_offset += (alt_l - ref_l);
						
						vcf_n_tmp++;
					}
					
					for (ite_i = 0; ite_i < alt_l - ref_l; ite_i++)
					{
						num_tmp = chr_start + vcf_pos - 1 - ite_i;
						if(alt_num[alt_l - 1 - ite_i] != ((buffer_ref_seq[num_tmp >> 5] >> ((31 - (num_tmp & 0X1f)) << 1)) & 0X3))
							break;
					}
					
					if(ite_i)
					{
						fprintf(stderr, "VCF_ADD INS B\n");

						vcf_buff[0] = vcf_pos - 1 + min_l - ite_i;//vcf_pos - 1 + min_l
						vcf_buff[1] = 1;
						vcf_buff[1] |= ((alt_l - ref_l) << vcf_bit_shift);
						vcf_buff[2] = vcf_seqs_offset;
						vcf_buff[3] = vcf_buff[0];
						fwrite(vcf_buff, 4, 4, file_p_pos[chr_no]);
						fwrite(alt_num + ref_l, 1, alt_l - ref_l, file_p[chr_no]);
						vcf_seqs_offset += (alt_l - ref_l);
						
						vcf_n_tmp++;
					}
				}
#endif
			}
			else if((ref_l > alt_l) && (memcmp(ref_num, alt_num, min_l) == 0))//D strncmp
			{
				vcf_buff[0] = vcf_pos - 1 + min_l;
				vcf_buff[1] = 2;
				vcf_buff[1] |= ((ref_l - alt_l) << vcf_bit_shift);
				vcf_buff[2] = vcf_seqs_offset;
				vcf_buff[3] = vcf_pos - 1 + ref_l;
				fwrite(vcf_buff, 4, 4, file_p_pos[chr_no]);
				//no write into vcf seq file

				vcf_n_tmp++;
				
				//fprintf("three\n");
				fprintf(fp_debug_vcf, "D: %u\n", vcf_pos);
					
#ifdef	VCF_ADD
				if((ref_l - alt_l < VCF_ADD_MAX) && (ref_l - alt_l > VCF_ADD_MIN))
				{
					for (ite_i = 0; ite_i < chr_posend_n; ite_i++)
						if(strcmp(chr_no_tmp, chr_names[ite_i] + 3) == 0)
							break;

					chr_start = chr_end_n[ite_i] - 1;

					fprintf(stderr, "ite_i: %u chr_names[ite_i]: %s chr_start: %u\n", ite_i, chr_names[ite_i], chr_start);
					
					for (ite_i = 0; ite_i < ref_l - alt_l; ite_i++)
					{
						num_tmp = chr_start + vcf_buff[3] + ite_i;
						char_tmp[ite_i] = ((buffer_ref_seq[num_tmp >> 5] >> ((31 - (num_tmp & 0X1f)) << 1)) & 0X3);
					}

					fprintf(stderr, "DEL: chr%u %u %u %u\n", chr_no, vcf_pos, ref_l, alt_l);
					
					//taacatttttatgtgttgctt ca tccagtttgctagagtttttggagatt

					//caagtcaataaatgtgatacaccaaa taaac agaatttaaaaaaaactca
					
					for(ite_i = 0; ite_i < ref_l - alt_l; ite_i++)
						fprintf(stderr, "%c", Dna5Tochar[char_tmp[ite_i]]);
					fprintf(stderr, "\n");

					for(ite_i = 0; ite_i < ref_l - alt_l; ite_i++)
						fprintf(stderr, "%c", Dna5Tochar[*(ref_num + alt_l + ite_i)]);
					fprintf(stderr, "\n");
					
					//ref_num[ref_l] = '\0';
					//char_tmp[ref_l - alt_l] = '\0';
					if(memcmp(ref_num + alt_l, char_tmp, ref_l - alt_l) == 0)
					{
						fprintf(stderr, "VCF_ADD A\n");

						vcf_buff[0] = vcf_buff[3];//vcf_pos - 1 + min_l
						vcf_buff[1] = 2;
						vcf_buff[1] |= ((ref_l - alt_l) << vcf_bit_shift);
						vcf_buff[2] = vcf_seqs_offset;
						vcf_buff[3] = vcf_buff[0] + ref_l - alt_l;
						fwrite(vcf_buff, 4, 4, file_p_pos[chr_no]);
						
						vcf_n_tmp++;
					}
					num_tmp_tmp = vcf_pos - 1 + min_l - ref_l + alt_l;
					for (ite_i = 0; ite_i < ref_l - alt_l; ite_i++)
					{
						num_tmp = chr_start + num_tmp_tmp + ite_i;
						char_tmp[ite_i] = ((buffer_ref_seq[num_tmp >> 5] >> ((31 - (num_tmp & 0X1f)) << 1)) & 0X3);
					}
					
					for(ite_i = 0; ite_i < ref_l - alt_l; ite_i++)
						fprintf(stderr, "%c", Dna5Tochar[char_tmp[ite_i]]);
					fprintf(stderr, "\n");

					for(ite_i = 0; ite_i < ref_l - alt_l; ite_i++)
						fprintf(stderr, "%c", Dna5Tochar[*(ref_num + alt_l + ite_i)]);
					fprintf(stderr, "\n");
					
					//ref_num[ref_l] = '\0';
					//char_tmp[ref_l - alt_l] = '\0';
					if(memcmp(ref_num + alt_l, char_tmp, ref_l - alt_l) == 0)
					{
						fprintf(stderr, "VCF_ADD B\n");

						vcf_buff[0] = num_tmp_tmp;//vcf_pos - 1 + min_l
						vcf_buff[1] = 2;
						vcf_buff[1] |= ((ref_l - alt_l) << vcf_bit_shift);
						vcf_buff[2] = vcf_seqs_offset;
						vcf_buff[3] = vcf_pos - 1 + min_l;
						fwrite(vcf_buff, 4, 4, file_p_pos[chr_no]);
						
						vcf_n_tmp++;
					}
					fflush(stderr);
				}

#endif
			}
			else     //mnp
			{
				vcf_buff[0] = vcf_pos - 1;
				vcf_buff[1] = 3;
				vcf_buff[1] |= (alt_l << vcf_bit_shift);
				vcf_buff[2] = vcf_seqs_offset;
				vcf_buff[3] = vcf_pos - 1 + ref_l;
				fwrite(vcf_buff, 4, 4, file_p_pos[chr_no]);
				fwrite(alt_num, 1, alt_l, file_p[chr_no]);
				vcf_seqs_offset += alt_l;

				vcf_n_tmp++;
				
				fprintf(fp_debug_vcf, "MNP: %u\n", vcf_pos);
			}

			//chr_vcf_n[chr_no]++;
			chr_vcf_n[chr_no] += vcf_n_tmp;
			
			pch_alt = strtok_r(NULL, ",", &saveptr_alt);
		}
		chr_vcf_seq_n[chr_no] = vcf_seqs_offset;

		//fprintf(fp_out_test, "%u ", vcf_line_n);
		//fputs(one_line, fp_out_test);
		//vcf_line_n++;
	}

	//fprintf("total number of vcf lines: %u\n", vcf_line_n);

	//debug
	if(fp_debug_vcf)	fclose(fp_debug_vcf);

	fclose(fp_in_vcf);

	//fclose(fp_out_test);

	for(ite_i = 0; ite_i < file_n; ite_i++)
	{
		if(file_p[ite_i] != NULL)
			fclose(file_p[ite_i]);
	}
	free(file_p);

	for(ite_i = 0; ite_i < file_n; ite_i++)
	{
		if(file_p_pos[ite_i] != NULL)
			fclose(file_p_pos[ite_i]);
	}
	free(file_p_pos);

	//if(chr_vcf_n)	free(chr_vcf_n);

}

void chr_vcf_index()
{
	uint8_t v_type = 0;
	uint8_t chr_no_f = 0;
	uint8_t ref_c = 0;
	uint8_t ref_cv_i = 0;
	uint8_t* alt_seq_chr = NULL;
#ifdef	INS_SEQ
	uint8_t alt_seq_ins[ALT_SEQ_MAX];
#endif
	char chr_no_tmp[CHR_NAME_MAX];
	char file_d[ROUTE_LENGTH_MAX];
	
	uint16_t chr_i = 0;
	uint16_t chr_no_i = 0;
	uint16_t chr_no = 0;

	uint32_t chr_vcf_total = 0;
	uint32_t ref_cv_n = 0;
	uint32_t rpos_cnt = 0;
	uint32_t apos_cnt = 0;
	uint32_t chr_start = 0;
	uint32_t chr_end = 0;
	uint32_t v_length = 0;
	uint32_t v_apos = 0;
	uint32_t v_rpos_end = 0;
	uint32_t v_rpos = 0;
	uint32_t v_rpos_pre = 0Xffffffff;
	uint32_t w_v_rpos = 0, w_v_rpos_end = 0, ref_length_p = 0;
	int vcf_sort_i = 0;
	int64_t ite_i = 0;

	uint32_t aseq_cnt = 0;
	uint64_t ref_cv = 0;
	uint64_t ref_cv_w = 0;
	uint64_t ref_cv_w_backup = 0;
	uint64_t result_num = 0;
	uint64_t a_size = 0;
	int64_t v_rpos_pre_snp = -1;
	uint64_t* buffer_ref_seq = NULL;

	FILE* fp_in_vcf_pos = NULL;
	FILE* fp_in_vcf_seq = NULL;
	FILE* fp_ref_seq = NULL;

	vcf_record_s* vcf_rs = NULL;

	//ref_length_p = ref_seq_n - 1;
	ref_length_p = ref_seq_n;

	fp_ref_seq = fopen(ref_seq, "rb");
	if (fp_ref_seq == NULL)
	{
		fputs ("File error opening the seq file\n",stderr);
		exit (1);
	}
	fprintf(stderr, "ref_seq_n: %u\n", ref_seq_n);//3095696031

	a_size = (ref_seq_n >> 5) + 1;
	buffer_ref_seq = (uint64_t* )calloc(a_size, 8);

	result_num = fread(buffer_ref_seq, 8, a_size, fp_ref_seq);//63676768

	fprintf(stderr, "result_num %u %u vcf_rpos %s\n", result_num, ref_seq_n, vcf_rpos);//96740501 3095696031

	if(fp_ref_seq)	fclose(fp_ref_seq);

	bit_shift_value = (1 << vcf_bit_shift) - 1;

	FILE* fp_out_rpos = fopen(vcf_rpos,"wb");
	if(fp_out_rpos == NULL)
	{
		fprintf(stderr, "Error of opening vcf rpos file\n");
		exit(1);
	}

	FILE* fp_out_aindex = fopen(vcf_aindex,"wb");
	if(fp_out_aindex == NULL)
	{
		fprintf(stderr, "Error of opening vcf ALT index file\n");
		exit(1);
	}

	FILE* fp_out_apos = fopen(vcf_apos,"wb");
	if(fp_out_apos == NULL)
	{
		fprintf(stderr, "Error of opening vcf ALT pos file\n");
		exit(1);
	}

	FILE* fp_out_aseq = fopen(vcf_aseq,"wb");
	if(fp_out_aseq == NULL)
	{
		fprintf(stderr, "Error of opening vcf aseq file\n");
		exit(1);
	}

	FILE* fp_out_rpos_end = fopen(vcf_rpos_end,"wb");
	if(fp_out_rpos_end == NULL)
	{
		fprintf(stderr, "Error of opening vcf rpos end file\n");
		exit(1);
	}

	FILE* fp_out_ref_seq_v = fopen(ref_vcf,"wb");
	if(fp_out_ref_seq_v == NULL)
	{
		fprintf(stderr, "Error of opening vcf ref seq file\n");
		exit(1);
	}

	FILE* fp_out_size = fopen(vcf_size,"wb");
	if(fp_out_size == NULL)
	{
		fprintf(stderr, "Error of opening vcf size file\n");
		exit(1);
	}

	FILE* fp_out_type = fopen(vcf_type,"wb");
	if(fp_out_type == NULL)
	{
		fprintf(stderr, "Error of opening vcf size file\n");
		exit(1);
	}


	FILE* fp_out_rpos_rec = fopen(vcf_rpos_rev,"wb");
	if(fp_out_rpos_rec == NULL)
	{
		printf("Error of opening vcf rpos rev file\n");
		exit(1);
	}

	FILE* fp_out_aindex_rec = fopen(vcf_aindex_rev,"wb");
	if(fp_out_aindex_rec == NULL)
	{
		printf("Error of opening vcf ALT index rev file\n");
		exit(1);
	}

	FILE* fp_out_apos_rec = fopen(vcf_apos_rev,"wb");
	if(fp_out_apos_rec == NULL)
	{
		printf("Error of opening vcf ALT pos rev file\n");
		exit(1);
	}

	FILE* fp_out_aseq_rec = fopen(vcf_aseq_rev,"wb");
	if(fp_out_aseq_rec == NULL)
	{
		printf("Error of opening vcf aseq rev file\n");
		exit(1);
	}

	FILE* fp_out_rpos_end_rec = fopen(vcf_rpos_end_rev,"wb");
	if(fp_out_rpos_end_rec == NULL)
	{
		printf("Error of opening vcf rpos end rev file\n");
		exit(1);
	}

	FILE* fp_out_ref_seq_v_rec = fopen(ref_vcf_rev,"wb");
	if(fp_out_ref_seq_v_rec == NULL)
	{
		printf("Error of opening vcf ref seq rev file\n");
		exit(1);
	}

	FILE* fp_out_size_rec = fopen(vcf_size_rev,"wb");
	if(fp_out_size_rec == NULL)
	{
		printf("Error of opening vcf size rev file\n");
		exit(1);
	}

	FILE* fp_out_type_rec = fopen(vcf_type_rev,"wb");
	if(fp_out_type_rec == NULL)
	{
		printf("Error of opening vcf size rev file\n");
		exit(1);
	}

	fprintf(stderr, "begin oganizing and indexing forward vcf info\n");

	FILE* fp_out_chr_rpos = fopen(vcf_chr_rpos,"wb");
	if(fp_out_chr_rpos == NULL)
	{
		fprintf(stderr, "Error of opening vcf chr rpos file\n");
		exit(1);
	}

	uint32_t* chr_rpos = (uint32_t* )calloc(chr_file_n_vcf, 4);

	fprintf(stderr, "There are %u(%u) chrs to index\n", chr_file_n_vcf - 1, chr_no_n);

	//debug
	for(chr_no_i = 0; chr_no_i < chr_no_n; chr_no_i++)
		fprintf(stderr, "%s\n", chr_no_a[chr_no_i]);

	chr_rpos[0] = 0;
	for(chr_i = 1; chr_i < chr_file_n_vcf; chr_i++)
	{
		strcpy(chr_no_tmp, chr_names[chr_i - 1] + 3);//extract the number here is only for hg19

		chr_no_f = 0;
		for(chr_no_i = 0; chr_no_i < chr_no_n; chr_no_i++)
			if(strcmp(chr_no_a[chr_no_i], chr_no_tmp) == 0)//chr_file_n_vcf chr_names[] from ref; chr_no_n chr_no_a[] from vcf
			{
				chr_no_f = 1;
				break;
			}

		chr_start = chr_end_n[chr_i - 1] - 1;
		chr_end = chr_end_n[chr_i] - 1;

		if(chr_no_f)
		{
			memset(file_d, 0, ROUTE_LENGTH_MAX);
			strcpy(file_d, index_route);
			strcat(file_d, vcf_pos_s);
			strcat(file_d, chr_no_tmp);

			fp_in_vcf_pos = fopen(file_d, "rb");

			memset(file_d, 0, ROUTE_LENGTH_MAX);
			strcpy(file_d, index_route);
			strcat(file_d, vcf_seq);
			strcat(file_d, chr_no_tmp);

			fp_in_vcf_seq = fopen(file_d, "rb");

			//debug
			//fprintf(stderr, "5 %u %u %s %s %s\n", chr_no_i, chr_vcf_seq_n[chr_no_i], chr_no_a[chr_no_i], chr_no_tmp, file_d);

			alt_seq_chr = (uint8_t* )calloc(chr_vcf_seq_n[chr_no_i], 1);

			fread(alt_seq_chr, 1, chr_vcf_seq_n[chr_no_i], fp_in_vcf_seq);

			chr_vcf_total = chr_vcf_n[chr_no_i];
			vcf_rs = (vcf_record_s* )malloc(chr_vcf_total * sizeof(vcf_record_s));
			//may need to change
			for(ite_i = 0; ite_i < chr_vcf_total; ite_i++)
				fread(vcf_rs[ite_i].vcf_buff_sort, 4, 4, fp_in_vcf_pos);

			//debug
			//fprintf(stderr, "step00\n");
			//fflush(stderr);

			qsort(vcf_rs, chr_vcf_total, sizeof(vcf_record_s ), compare_vcf_sort);

			//debug
			//fprintf(stderr, "step01 %u\n", chr_vcf_total);
			//fflush(stderr);

			for(vcf_sort_i = 0; vcf_sort_i < chr_vcf_total; vcf_sort_i++)
			{
				v_rpos = vcf_rs[vcf_sort_i].vcf_buff_sort[0] + chr_start;
				v_type = (vcf_rs[vcf_sort_i].vcf_buff_sort[1]) & bit_shift_value;
				v_length = vcf_rs[vcf_sort_i].vcf_buff_sort[1] >> vcf_bit_shift;
				v_apos = vcf_rs[vcf_sort_i].vcf_buff_sort[2];
				v_rpos_end = vcf_rs[vcf_sort_i].vcf_buff_sort[3] + chr_start;

				//for test
				//fprintf(stderr, "sort: %u %u %u %u %u %u\n", v_rpos, v_type, v_length, v_apos, v_rpos_end, v_rpos_pre_snp);
				//fflush(stderr);
				//2147484023 0 1 70606 2147484024 2147483601
				//67108875 * 8 = 536871005
				//773924008
				if(v_type)
				{
					if(v_rpos_pre != v_rpos)
					{
						fwrite(&v_rpos, 4, 1, fp_out_rpos);
						rpos_cnt++;
						fwrite(&apos_cnt, 4, 1, fp_out_aindex);
					}

					fwrite(&aseq_cnt, 4, 1, fp_out_apos);
					apos_cnt++;
					fwrite(&v_rpos_end, 4, 1, fp_out_rpos_end);

					if(v_type != 2)
					{
						fwrite(alt_seq_chr + v_apos, 1, v_length, fp_out_aseq);
						aseq_cnt += v_length;
					}
					v_rpos_pre = v_rpos;

					//write type
					fwrite(&v_type, 1, 1, fp_out_type);

				}
				else //SNP write SNPs and characters between two SNPs into ref_vcf.seq
				{
					//debug
					//fprintf(stderr, "step111\n");
					//fflush(stderr);

					for(ite_i = v_rpos_pre_snp + 1; ite_i < v_rpos; ite_i++)
					{
						ref_c = ((buffer_ref_seq[ite_i >> 5] >> ((31 - (ite_i & 0X1f)) << 1)) & 0X3);
						ref_cv = (1 << ref_c);

						ref_cv_w |= (ref_cv << (((uint64_t )(15 - ref_cv_i)) << 2));
						ref_cv_i++;
						ref_cv_n++;
						if(ref_cv_i == 16)
						{
							fwrite(&ref_cv_w, 8, 1, fp_out_ref_seq_v);//ref seq between two snps
							ref_cv_w = 0;
							ref_cv_i = 0;
						}
					}

					//debug
					//fprintf(stderr, "step222\n");

					if(v_rpos_pre_snp != v_rpos)
					{
						ref_c = ((buffer_ref_seq[v_rpos >> 5] >> ((31 - (v_rpos & 0X1f)) << 1)) & 0X3);
						ref_cv = (1 << ref_c);

						//debug
						//fprintf(stderr, "%u\n", alt_seq_chr[v_apos]);
						//fflush(stderr);

						ref_cv |= (1 << alt_seq_chr[v_apos]);

						ref_cv_w |= (ref_cv << (((uint64_t )(15 - ref_cv_i)) << 2));

						ref_cv_w_backup = ref_cv_w;

						ref_cv_i++;

						ref_cv_n++;

						if(ref_cv_i == 16)
						{
							fwrite(&ref_cv_w, 8, 1, fp_out_ref_seq_v);
							ref_cv_w = 0;
							ref_cv_i = 0;
						}
					}
					else
					{
						if(ref_cv_i)
						{
							ref_cv_w |= (((uint64_t )(1 << alt_seq_chr[v_apos])) << (((uint64_t )(16 - ref_cv_i)) << 2));
						}
						else
						{
							ref_cv_w_backup |= (1 << alt_seq_chr[v_apos]);
							fseek(fp_out_ref_seq_v, -8, SEEK_CUR);
							fwrite(&ref_cv_w_backup, 8, 1, fp_out_ref_seq_v);
						}
					}
					v_rpos_pre_snp = v_rpos;
				}
			}

			//debug
			//fprintf(stderr, "step02\n");

			//free
			if(vcf_rs)	free(vcf_rs);
			if(alt_seq_chr)	free(alt_seq_chr);

			//debug
			//fprintf(stderr, "step11\n");
			//fflush(stderr);
		}
		else  //write SNP seq directly
		{
			//debug
			//fprintf(stderr, "step33\n");
			//fflush(stderr);

			for(ite_i = chr_start; ite_i < chr_end; ite_i++)
			{
				ref_c = ((buffer_ref_seq[ite_i >> 5] >> ((31 - (ite_i & 0X1f)) << 1)) & 0X3);
				ref_cv = (1 << ref_c);

				ref_cv_w |= (ref_cv << (((uint64_t )(15 - ref_cv_i)) << 2));
				ref_cv_i++;

				ref_cv_n++;

				if(ref_cv_i == 16)
				{
					fwrite(&ref_cv_w, 8, 1, fp_out_ref_seq_v);
					ref_cv_w = 0;
					ref_cv_i = 0;
				}
			}
			//debug
			//fprintf(stderr, "step44\n");
			//fflush(stderr);
		}

		chr_rpos[chr_i] = rpos_cnt;

		//debug
		//fprintf(stderr, "step55\n");
		//fflush(stderr);
	}

	if(chr_no_f)
	{
		for(ite_i = v_rpos + 1; ite_i < chr_end; ite_i++)
		{
			ref_c = ((buffer_ref_seq[ite_i >> 5] >> ((31 - (ite_i & 0X1f)) << 1)) & 0X3);
			ref_cv = (1 << ref_c);

			ref_cv_w |= (ref_cv << (((uint64_t )(15 - ref_cv_i)) << 2));
			ref_cv_i++;
			ref_cv_n++;
			if(ref_cv_i == 16)
			{
				fwrite(&ref_cv_w, 8, 1, fp_out_ref_seq_v);
				ref_cv_w = 0;
				ref_cv_i = 0;
			}
		}
	}

	fwrite(chr_rpos, chr_file_n_vcf, 4, fp_out_chr_rpos);

	fwrite(&ref_cv_w, 8, 1, fp_out_ref_seq_v);
	fwrite(&aseq_cnt, 4, 1, fp_out_apos);
	fwrite(&apos_cnt, 4, 1, fp_out_aindex);

	//ftell: 136975491 68487752
	//fprintf("ftell: %u %ld\n", ref_cv_n, ftell(fp_out_ref_seq_v));

	fwrite(&rpos_cnt, 4, 1, fp_out_size);
	fwrite(&apos_cnt, 4, 1, fp_out_size);
	fwrite(&aseq_cnt, 4, 1, fp_out_size);
	fwrite(&ref_cv_n, 4, 1, fp_out_size);

	if(fp_in_vcf_pos)	fclose(fp_in_vcf_pos);
	if(fp_in_vcf_seq)	fclose(fp_in_vcf_seq);

	if(fp_out_rpos)	fclose(fp_out_rpos);
	if(fp_out_aindex)	fclose(fp_out_aindex);
	if(fp_out_apos)	fclose(fp_out_apos);
	if(fp_out_aseq)	fclose(fp_out_aseq);
	if(fp_out_rpos_end)	fclose(fp_out_rpos_end);
	if(fp_out_ref_seq_v)	fclose(fp_out_ref_seq_v);
	if(fp_out_size)	fclose(fp_out_size);
	if(fp_out_type)	fclose(fp_out_type);


	rpos_cnt = 0;
	apos_cnt = 0;
	aseq_cnt = 0;
	ref_cv_w = 0;
	ref_cv_i = 0;
	ref_cv_n = 0;
	v_rpos_pre = 0Xffffffff;
	v_rpos_pre_snp = chr_end_n[chr_file_n_vcf - 1] - 1;
	chr_rpos[0] = 0;
	for(chr_i = chr_file_n_vcf - 1; chr_i > 0; chr_i--)
	{
		strcpy(chr_no_tmp, chr_names[chr_i - 1] + 3);//

		//fprintf(stderr, "1 %s %s\n", chr_no_tmp, chr_names[chr_i - 1]);

		chr_no_f = 0;
		for(chr_no_i = 0; chr_no_i < chr_no_n; chr_no_i++)
			if(strcmp(chr_no_a[chr_no_i], chr_no_tmp) == 0)
			{
				chr_no_f = 1;
				break;
			}

		//fprintf(stderr, "2 %u\n", chr_no_f);

		chr_start = chr_end_n[chr_i - 1] - 1;
		chr_end = chr_end_n[chr_i] - 1;

		//fprintf(stderr, "3 %s\n", index_route);

		if(chr_no_f)
		{
			memset(file_d, 0, ROUTE_LENGTH_MAX);
			strcpy(file_d, index_route);
			strcat(file_d, vcf_pos_s);
			strcat(file_d, chr_no_tmp);

			fp_in_vcf_pos = fopen(file_d, "rb");

			//fprintf(stderr, "4\n");

			memset(file_d, 0, ROUTE_LENGTH_MAX);
			strcpy(file_d, index_route);
			strcat(file_d, vcf_seq);
			strcat(file_d, chr_no_tmp);

			fp_in_vcf_seq = fopen(file_d, "rb");

			//fprintf(stderr, "5 %u %u %s %s\n", chr_no_i, chr_vcf_seq_n[chr_no_i], chr_no_a[chr_no_i], chr_no_tmp);

			alt_seq_chr = (uint8_t* )calloc(chr_vcf_seq_n[chr_no_i], 1);

			fread(alt_seq_chr, 1, chr_vcf_seq_n[chr_no_i], fp_in_vcf_seq);

			chr_vcf_total = chr_vcf_n[chr_no_i];
			vcf_rs = (vcf_record_s* )malloc(chr_vcf_total * sizeof(vcf_record_s));
			for(ite_i = 0; ite_i < chr_vcf_total; ite_i++)
				fread(vcf_rs[ite_i].vcf_buff_sort, 4, 4, fp_in_vcf_pos);

			qsort(vcf_rs, chr_vcf_total, sizeof(vcf_record_s ), compare_vcf_sort_rec);

			//fprintf(stderr, "6 %u\n", chr_vcf_total);

			for(vcf_sort_i = chr_vcf_total - 1; vcf_sort_i >= 0; vcf_sort_i--)
			{
				//fprintf(stderr, "7 %u\n", vcf_sort_i);

				v_rpos_end = vcf_rs[vcf_sort_i].vcf_buff_sort[0] + chr_start - 1;
				v_type = (vcf_rs[vcf_sort_i].vcf_buff_sort[1]) & bit_shift_value;
				v_length = vcf_rs[vcf_sort_i].vcf_buff_sort[1] >> vcf_bit_shift;
				v_apos = vcf_rs[vcf_sort_i].vcf_buff_sort[2];
				v_rpos = vcf_rs[vcf_sort_i].vcf_buff_sort[3] + chr_start - 1;

				//for test
				//printf("sorted: %u %u %u %u %u\n", v_rpos, v_type, v_length, v_apos, v_rpos_end);

				if(v_type)
				{
					//fprintf(stderr, "9\n");

					if(v_rpos_pre != v_rpos)
					{
						w_v_rpos = ref_length_p - v_rpos;
						fwrite(&w_v_rpos, 4, 1, fp_out_rpos_rec);
						rpos_cnt++;
						fwrite(&apos_cnt, 4, 1, fp_out_aindex_rec);
					}

					fwrite(&aseq_cnt, 4, 1, fp_out_apos_rec);
					apos_cnt++;
					w_v_rpos_end = ref_length_p - v_rpos_end;
					fwrite(&w_v_rpos_end, 4, 1, fp_out_rpos_end_rec);

					if(v_type != 2)
					{
#ifdef	INS_SEQ
						for(ite_i = 0; ite_i < v_length; ite_i++)
							alt_seq_ins[ite_i] = alt_seq_chr[v_apos + v_length - 1 - ite_i];
						fwrite(alt_seq_ins, 1, v_length, fp_out_aseq_rec);
#else
						fwrite(alt_seq_chr + v_apos, 1, v_length, fp_out_aseq_rec);
#endif
						aseq_cnt += v_length;
					}
					v_rpos_pre = v_rpos;

					//write type
					fwrite(&v_type, 1, 1, fp_out_type_rec);

					//fprintf(stderr, "10\n");
				}
				else //SNP
				{

					//fprintf(stderr, "11 %u %u\n", v_rpos_pre_snp, v_rpos);

					for(ite_i = v_rpos_pre_snp - 1; ite_i > v_rpos; ite_i--)
					{
						ref_c = ((buffer_ref_seq[ite_i >> 5] >> ((31 - (ite_i & 0X1f)) << 1)) & 0X3);
						ref_cv = (1 << ref_c);

						if(ref_cv == 0Xf)	fprintf(stderr, "wrong %u\n", ite_i);

						ref_cv_w |= (ref_cv << (((uint64_t )(15 - ref_cv_i)) << 2));
						ref_cv_i++;
						ref_cv_n++;
						if(ref_cv_i == 16)
						{
							fwrite(&ref_cv_w, 8, 1, fp_out_ref_seq_v_rec);
							ref_cv_w = 0;
							ref_cv_i = 0;
						}
					}

					//fprintf(stderr, "12\n");

					if(v_rpos_pre_snp != v_rpos)
					{
						ref_c = ((buffer_ref_seq[v_rpos >> 5] >> ((31 - (v_rpos & 0X1f)) << 1)) & 0X3);
						ref_cv = (1 << ref_c);
						ref_cv |= (1 << alt_seq_chr[v_apos]);

						//debug
						if(ref_cv == 0Xf)	fprintf(stderr, "wrong %u\n", ite_i, ref_c, alt_seq_chr[v_apos]);

						ref_cv_w |= (ref_cv << (((uint64_t )(15 - ref_cv_i)) << 2));

						ref_cv_w_backup = ref_cv_w;

						ref_cv_i++;

						ref_cv_n++;

						if(ref_cv_i == 16)
						{
							fwrite(&ref_cv_w, 8, 1, fp_out_ref_seq_v_rec);
							ref_cv_w = 0;
							ref_cv_i = 0;
						}
					}
					else
					{
						if(ref_cv_i)
						{
							ref_cv_w |= (((uint64_t )(1 << alt_seq_chr[v_apos])) << (((uint64_t )(16 - ref_cv_i)) << 2));
						}
						else
						{
							ref_cv_w_backup |= (1 << alt_seq_chr[v_apos]);
							fseek(fp_out_ref_seq_v_rec, -8, SEEK_CUR);
							fwrite(&ref_cv_w_backup, 8, 1, fp_out_ref_seq_v_rec);
						}
					}
					v_rpos_pre_snp = v_rpos;

					//fprintf(stderr, "13\n");
				}
			}

			//fprintf(stderr, "7\n");

			//free
			if(vcf_rs)	free(vcf_rs);
			if(alt_seq_chr)	free(alt_seq_chr);

		}
		else
		{
			for(ite_i = chr_end - 1; ite_i >= chr_start; ite_i--)
			{
				ref_c = ((buffer_ref_seq[ite_i >> 5] >> ((31 - (ite_i & 0X1f)) << 1)) & 0X3);
				ref_cv = (1 << ref_c);

				//debug
				if(ref_cv == 0Xf)	fprintf(stderr, "wrong %u\n", ite_i, ref_c, alt_seq_chr[v_apos]);

				ref_cv_w |= (ref_cv << (((uint64_t )(15 - ref_cv_i)) << 2));
				ref_cv_i++;

				ref_cv_n++;

				if(ref_cv_i == 16)
				{
					fwrite(&ref_cv_w, 8, 1, fp_out_ref_seq_v_rec);
					ref_cv_w = 0;
					ref_cv_i = 0;
				}
			}

		}

		chr_rpos[chr_file_n_vcf - chr_i] = rpos_cnt;

	}

	if(chr_no_f)
	{
		for(ite_i = v_rpos - 1; ite_i >= 0; ite_i--)
		{
			ref_c = ((buffer_ref_seq[ite_i >> 5] >> ((31 - (ite_i & 0X1f)) << 1)) & 0X3);
			ref_cv = (1 << ref_c);

			//debug
			if(ref_cv == 0Xf)	fprintf(stderr, "wrong %u\n", ite_i, ref_c, alt_seq_chr[v_apos]);

			ref_cv_w |= (ref_cv << (((uint64_t )(15 - ref_cv_i)) << 2));
			ref_cv_i++;
			ref_cv_n++;
			if(ref_cv_i == 16)
			{
				fwrite(&ref_cv_w, 8, 1, fp_out_ref_seq_v_rec);
				ref_cv_w = 0;
				ref_cv_i = 0;
			}
		}
	}
//fprintf(stderr, "4\n");
	fwrite(chr_rpos, chr_file_n_vcf, 4, fp_out_chr_rpos);

	fwrite(&ref_cv_w, 8, 1, fp_out_ref_seq_v_rec);
	fwrite(&aseq_cnt, 4, 1, fp_out_apos_rec);
	fwrite(&apos_cnt, 4, 1, fp_out_aindex_rec);

	fwrite(&rpos_cnt, 4, 1, fp_out_size_rec);
	fwrite(&apos_cnt, 4, 1, fp_out_size_rec);
	fwrite(&aseq_cnt, 4, 1, fp_out_size_rec);
	fwrite(&ref_cv_n, 4, 1, fp_out_size_rec);

	if(buffer_ref_seq)	free(buffer_ref_seq);
	if(fp_out_rpos_rec)	fclose(fp_out_rpos_rec);
	if(fp_out_aindex_rec)	fclose(fp_out_aindex_rec);
	if(fp_out_apos_rec)	fclose(fp_out_apos_rec);
	if(fp_out_aseq_rec)	fclose(fp_out_aseq_rec);
	if(fp_out_rpos_end_rec)	fclose(fp_out_rpos_end_rec);
	if(fp_out_ref_seq_v_rec)	fclose(fp_out_ref_seq_v_rec);
	if(fp_out_size_rec)	fclose(fp_out_size_rec);
	if(fp_out_type_rec)	fclose(fp_out_type_rec);
	if(fp_out_chr_rpos)	fclose(fp_out_chr_rpos);

}

int compare_vcf_sort(const void * a, const void * b)
{
	vcf_record_s* vs1 = (vcf_record_s* )a;
	vcf_record_s* vs2 = (vcf_record_s* )b;

	uint32_t pos1 = vs1->vcf_buff_sort[0];
	uint32_t pos2 = vs2->vcf_buff_sort[0];

	if(pos1 > pos2)
		return 1;
	else if(pos1 < pos2)
		return -1;
	else	return 0;
}

int compare_vcf_sort_rec(const void * a, const void * b)
{
	vcf_record_s* vs1 = (vcf_record_s* )a;
	vcf_record_s* vs2 = (vcf_record_s* )b;

	uint32_t pos1 = vs1->vcf_buff_sort[3];
	uint32_t pos2 = vs2->vcf_buff_sort[3];

	if(pos1 > pos2)
		return 1;
	else if(pos1 < pos2)
		return -1;
	else	return 0;
}

void upper_string(char *string)
{
	while(*string)
	{
		if ( *string >= 'a' && *string <= 'z' )
		{
			*string = *string - 32;
		}
		string++;
	}
}

void load_reffile_kmer_fa()
{
	fp_n = fopen(N_route, "w");
	if(fp_n == NULL)    printf("cannot open N statistical file\n");

	uint8_t n_end_f = 0;
	uint8_t n_w_f = 0;
	uint8_t line_i = 0;
	uint8_t last_kmer_f = 0;
	uint8_t fisrt_chr = 0;
	uint8_t line_l = 0;
	uint8_t l_p = 0;

	char input = charN;
	char output;

	char one_line[REF_FASTA_LINE + 1] = "";
	char kmer_pre[KT_LENGTH_MAX + 1] = "";
	char kmer_com[KT_LENGTH_MAX + REF_FASTA_LINE + 1] = "";
	char kmer[KT_LENGTH_MAX + 1];
	char kmer_f[KF_LENGTH_MAX] = "";
	char kmer_s[KT_LENGTH_MAX] = "";
	char file_d[KF_LENGTH_MAX + ROUTE_LENGTH_MAX] = "";

	uint32_t n_cnt = 0;
	uint32_t n_tr = 0;
	uint32_t kmer_i = 0;
	uint32_t line_t_c = 0;
	uint32_t i = 0;
	uint32_t chr_i = 0;
	uint32_t chr_name_n = 0;
	uint32_t file_n = 0;
	uint32_t add = 0;

	uint64_t ref_seq_buffer[128];

	FILE* fpin_ref = NULL;


	file_n = ((uint32_t )1 << (f << 1));
	//modification
	//file_n = ((1 << (f << 1)) << 1);
	fprintf(stderr, "%u files to be written in folder div\n", file_n);
	fflush(stderr);

	FILE** file_p = (FILE** )malloc(file_n * sizeof(FILE* ));
	if(file_p == NULL)
	{
		printf("Fail to allocate memory for file_p\n");
		fflush(stdout);
		exit(1);
	}
	for(i = 0; i < file_n; i++)
		file_p[i] = NULL;

	FILE* file_sta = NULL;
	file_sta = fopen(filename_sta,"w");
	if (file_sta == NULL)
	{
		fputs ("File error of creating statistical file\n",stderr);
		fflush(stdout);
		exit(1);
	}

#ifdef	UNPIPATH_OFF_K20
	uint32_t write_buff[4];
	uint64_t pos = START_POS_REF;
	uint64_t pos_tr = 0;
	uint64_t line_tol = 0;
	uint64_t* file_line_cnt = (uint64_t* )calloc(file_n, 8);
#else
	uint32_t write_buff[3];
	uint32_t pos = START_POS_REF;
	uint32_t pos_tr = 0;
	uint32_t line_tol = 0;
	uint32_t* file_line_cnt = (uint32_t* )calloc(file_n, 4);
#endif

	if(file_line_cnt == NULL)
	{
		printf("Fail to allocate memory for file_line_cnt\n");
		fflush(stdout);
		exit(1);
	}

	fp_ref_seq = fopen(ref_seq, "wb");

	//begin loading file and loading kmer into graph
	fpin_ref = fopen(filename_ref, "r");
	if (fpin_ref == NULL)
	{
		fputs ("File error of reading ref file and wrong ref file name or route\n",stderr);
		exit(1);
	}

	l_p = 0;

	fprintf(stderr,"3\n");

#ifdef	CHR_NAME_SPLIT
	char* pch = NULL;
	char* saveptr = NULL;
	char tmp_chr_name[REF_FASTA_LINE + 1];
#endif

	fisrt_chr = 1;
	chr_posend_n = 1;
	chr_end_n[0] = START_POS_REF + 1;
	ref_seq_n = START_POS_REF;

	while ((!feof(fpin_ref)) && (fgets(one_line, REF_FASTA_LINE + 2, fpin_ref) != NULL))
	{
		if(strstr(one_line,identifier) != NULL)
		{
			one_line[strlen(one_line) - 1] = '\0';

#ifdef	CHR_NAME_SPLIT
			strcpy(tmp_chr_name, one_line + 1);
			pch = strtok_r(tmp_chr_name, " ", &saveptr);

			strcpy(chr_names[chr_name_n++], pch);
#else
			strcpy(chr_names[chr_name_n++], one_line + 1);
#endif

			if(fisrt_chr)
			{
				fisrt_chr = 0;
				continue;
			}

			//write the last kmer
			if(kmer_pre[0] != '\0')
			{
				seq_exact(kmer_pre,kmer,kmer_i,0,k_t)

				last_kmer_f = 1;

				if((strstr(kmer, char_N) != NULL) || (strstr(kmer, char_n) != NULL))
				{
					input = charN;
					last_kmer_f = 0;

					//N number
					++n_cnt;
					n_w_f = 1;
					n_end_f = 1;
				}

				if(last_kmer_f)
				{
					output = charN;

					//do something here, kmer[], pos + 1
					strncpy(kmer_f, (const char* )kmer, f);
					kmer_address(kmer_f, f, kmer_i, add)

					upper_string(kmer_f);

					memset(file_d, 0, KF_LENGTH_MAX + ROUTE_LENGTH_MAX);
					strcpy(file_d, (const char* )filename_div);
					strcat(file_d, (const char* )kmer_f);
					strcat(file_d, suff);

					if(file_p[add] == NULL)
					{
						file_p[add] = fopen(file_d, "wb");
						if (file_p[add] == NULL)
						{
							fputs ("File error of creating divisional file\n",stderr);
							exit(1);
						}
					}

					seq_exact(kmer, kmer_s, kmer_i, f, strlen(kmer));

					//binary file
					//must assign in this order
					kmer_bit32a(kmer_s, k_t - f, kmer_i, write_buff, 3)
					write_buff[1] |= (input << 24);
					write_buff[1] |= (output << 16);

#ifdef	UNPIPATH_OFF_K20
					write_buff[3] = (pos + 1) >> 32;
					write_buff[0] = (pos + 1) & 0Xffffffff;
					fwrite(write_buff, 4, 4, file_p[add]);
#else
					write_buff[0] = pos + 1;
					fwrite(write_buff, 4, 3, file_p[add]);
					//fwrite(write_buff, sizeof(write_buff ), 1, file_p[add]);
#endif

#ifdef	DEBUG_INPUT
					if(write_buff[2] == 8258096)
						printf("\n\n2. %s:\nline: %s\n%s\n%u\n\n\n", tmp_chr_name, one_line, kmer, (uint16_t )write_buff[1]);
#endif

					//
					file_line_cnt[add]++;
					line_tol++;
				}
			}


			memset(kmer_pre, 0, KT_LENGTH_MAX + 1);
			kmer_pre[0] = '\0';

			input = charN;

			pos += k_t;

			//end last kmer


			chr_end_n[chr_posend_n++] = pos + 1;

			fprintf(stderr, "Has finished loading %s\n", chr_names[chr_name_n - 2]);
			fflush(stderr);

			l_p = 0;
		}
		else
		{
			line_t_c++;

			line_l = strlen(one_line);

			if(one_line[line_l - 1] == '\n')
			{
				one_line[line_l - 1] = '\0';
				line_l--;
			}

			//ref_load
			for(line_i = 0; line_i < line_l; line_i++)
			{
				//modification
				one_line[line_i] = charTochar[(uint8_t )one_line[line_i]];

				ref_seq_buffer[(ref_seq_n & 0Xfff) >> 5] |= (((uint64_t )charToDna5_N2[(uint8_t )one_line[line_i]]) << ((31 - ((ref_seq_n & 0Xfff) & 0X1f)) << 1));

				++ref_seq_n;
				if((ref_seq_n & 0Xfff) == 0)
				{
					fwrite(ref_seq_buffer, 8, 128, fp_ref_seq);
					memset(ref_seq_buffer, 0, 1024);
				}
			}

			memset(kmer_com, 0, KT_LENGTH_MAX + REF_FASTA_LINE);

			seq_exact_2(kmer_pre,kmer_com,kmer_i,0)
			seq_exact_2(one_line,kmer_com,kmer_i,l_p)

			for(line_i = 0; line_i < line_l - k_t + l_p; line_i++,pos++)
			{
				seq_exact(kmer_com,kmer,kmer_i,line_i,line_i + k_t)

				if((strstr(kmer, char_N) != NULL) || (strstr(kmer, char_n) != NULL))
				{
					memset(kmer_pre, 0, KT_LENGTH_MAX + 1);
					kmer_pre[0] = '\0';

					input = charN;

					//N number
					++n_cnt;
					n_w_f = 1;

					continue;
				}

				//N position and number
				if(n_w_f)
				{
					pos_tr = pos + 1;

					if(n_end_f)
					{
						n_tr = n_cnt;
						n_end_f = 0;
					}
					else	n_tr = n_cnt - (k_t - 1);


					if(n_tr <= 0)
					{
						printf("wrong n %u %u\n", n_tr, n_cnt);
						exit(1);
					}
					fwrite(&pos_tr, 4, 1, fp_n);
					fwrite(&n_tr, 4, 1, fp_n);


					n_cnt = 0;
					n_w_f = 0;
				}

				output = kmer_com[line_i + k_t];

				//do something here, kmer[], pos + 1
				strncpy(kmer_f, (const char* )kmer, f);
				kmer_address(kmer_f, f, kmer_i, add)

				upper_string(kmer_f);

				memset(file_d, 0, KF_LENGTH_MAX + ROUTE_LENGTH_MAX);
				strcpy(file_d, (const char* )filename_div);
				strcat(file_d, (const char* )kmer_f);
				strcat(file_d, suff);

				if(file_p[add] == NULL)
				{
					file_p[add] = fopen(file_d, "wb");
					if (file_p[add] == NULL)
					{
						fputs ("File error of creating divisional file\n",stderr);
						exit(1);
					}
				}

				seq_exact(kmer, kmer_s, kmer_i, f, strlen(kmer));

				//binary file
				//must assign in this order
				kmer_bit32a(kmer_s, k_t - f, kmer_i, write_buff, 3)
				write_buff[1] |= (input << 24);
				write_buff[1] |= (output << 16);

#ifdef	UNPIPATH_OFF_K20
				write_buff[3] = (pos + 1) >> 32;
				write_buff[0] = (pos + 1) & 0Xffffffff;
				fwrite(write_buff, 4, 4, file_p[add]);
#else
				write_buff[0] = pos + 1;
				fwrite(write_buff, 4, 3, file_p[add]);
				//fwrite(write_buff, sizeof(write_buff ), 1, file_p[add]);
#endif

#ifdef	DEBUG_INPUT
				if(write_buff[2] == 8258096)
					printf("\n\n%2. s:\nline: %s\n%s\n%u\n\n\n", tmp_chr_name, one_line, kmer, write_buff[1]);
#endif
				//

				file_line_cnt[add]++;
				line_tol++;
				//

				input = kmer_com[line_i];
			}

			seq_exact(kmer_com, kmer_pre, kmer_i, line_l - k_t + l_p, line_l+l_p)

			l_p = k_t;
		}
	}

	//write the last kmer
	if(kmer_pre[0] != '\0')
	{
		seq_exact(kmer_pre,kmer,kmer_i,0,k_t)

		last_kmer_f = 1;

		if((strstr(kmer, char_N) != NULL) || (strstr(kmer, char_n) != NULL))
		{
			input = charN;
			last_kmer_f = 0;

			//N number
			++n_cnt;
			n_w_f = 1;
			n_end_f = 1;
		}

		if(last_kmer_f)
		{
			output = charN;

			//do something here, kmer[], pos + 1
			strncpy(kmer_f, (const char* )kmer, f);
			kmer_address(kmer_f, f, kmer_i, add)

			upper_string(kmer_f);

			memset(file_d, 0, KF_LENGTH_MAX + ROUTE_LENGTH_MAX);
			strcpy(file_d, (const char* )filename_div);
			strcat(file_d, (const char* )kmer_f);
			strcat(file_d, suff);

			if(file_p[add] == NULL)
			{
				file_p[add] = fopen(file_d, "wb");
				if (file_p[add] == NULL)
				{
					fputs ("File error of creating divisional file\n",stderr);
					exit(1);
				}
			}

			seq_exact(kmer, kmer_s, kmer_i, f, strlen(kmer));


			//binary file
			//must assign in this order
			kmer_bit32a(kmer_s, k_t - f, kmer_i, write_buff, 3)
			write_buff[1] |= (input << 24);
			write_buff[1] |= (output << 16);

#ifdef	UNPIPATH_OFF_K20
			write_buff[3] = (pos + 1) >> 32;
			write_buff[0] = (pos + 1) & 0Xffffffff;
			fwrite(write_buff, 4, 4, file_p[add]);
#else
			write_buff[0] = pos + 1;
			fwrite(write_buff, 4, 3, file_p[add]);
			//fwrite(write_buff, sizeof(write_buff ), 1, file_p[add]);
#endif

#ifdef	DEBUG_INPUT
			if(write_buff[2] == 8258096)
				printf("\n\n2. %s:\nline: %s\n%s\n%u\n\n\n", tmp_chr_name, one_line, kmer, write_buff[1]);
#endif
			//

			file_line_cnt[add]++;
			line_tol++;
		}
	}


	pos += k_t;
	chr_end_n[chr_posend_n++] = pos + 1;

	fprintf(stderr, "chr_end_n: ");
	for(i = 0; i < chr_posend_n; i++)
		fprintf(stderr, "%u ", chr_end_n[i]);

	fprintf(stderr, "\n");

	fflush(stderr);

	printf("Has finished loading %s\n", chr_names[chr_name_n - 1]);
#ifdef	UNPIPATH_OFF_K20
	printf("total number of lines: %"PRId64"\n",line_tol);
	printf("total number of positions on ref seq: %"PRId64"\n", pos+1);
#else
	printf("total number of lines: %u\n",line_tol);
	printf("total number of positions on ref seq: %u\n", pos+1);
#endif
	fclose(fpin_ref);

	fwrite(ref_seq_buffer, 8, ((ref_seq_n & 0Xfff) >> 5) + 1, fp_ref_seq);
	fclose(fp_ref_seq);

	fprintf(stderr, "ref_seq_n: %"PRId64"\n", ref_seq_n);
	fflush(stderr);
	
	printf("total number of chars to allocate for ref seq: %"PRId64"\n", (((ref_seq_n >> 5) + 1) << 3));

	fflush(stdout);

	FILE* fp_chr = fopen(unichr,"w");
	if (fp_chr == NULL)
	{
		fputs ("File error of creating chr info file\n",stderr);
		exit(1);
	}

#ifdef	UNPIPATH_OFF_K20
	for(chr_i = 1; chr_i < chr_posend_n; chr_i++)
	{
		fprintf(fp_chr,"%s\n%"PRId64"\n",chr_names[chr_i - 1], chr_end_n[chr_i]);
	}
	fclose(fp_chr);

	for(i = 0; i < file_n; i++)
		fprintf(file_sta,"%"PRId64"\n",file_line_cnt[i]);
	fprintf(file_sta,"\n%"PRId64"",line_tol);
	fflush(file_sta);
#else
	for(chr_i = 1; chr_i < chr_posend_n; chr_i++)
	{
		fprintf(fp_chr,"%s\n%u\n",chr_names[chr_i - 1], chr_end_n[chr_i]);
	}
	fclose(fp_chr);

	for(i = 0; i < file_n; i++)
		fprintf(file_sta,"%u\n",file_line_cnt[i]);
	fprintf(file_sta,"\n%u",line_tol);
	fflush(file_sta);
#endif

	for(chr_i = 1; chr_i < chr_posend_n; chr_i++)
	{
		fprintf(stderr,"chr_names[%u]=\"%s\"\n",chr_i - 1, chr_names[chr_i - 1]);
	}
	for(chr_i = 0; chr_i < chr_posend_n; chr_i++)
	{
		fprintf(stderr,"chr_end_n[%u]=%"PRId64"\n",chr_i, chr_end_n[chr_i]);
	}
	
	for(i = 0; i < file_n; i++)
	{
		if(file_p[i] != NULL)
			fclose(file_p[i]);
	}

	free(file_p);

	free(file_line_cnt);

	fclose(file_sta);

	fclose(fp_n);

	if(tree_flag)	chr_file_n_vcf = chr_posend_n;

	fprintf(stderr, "debug: %u %u %u\n", chr_posend_n, chr_name_n, chr_file_n_vcf);

}


int compare (const void * a, const void * b)
{
	k_p* kp1 = (k_p* )a;
	k_p* kp2 = (k_p* )b;

	uint16_t kmerf1 = (uint16_t )((kp1->kp)[1]);
	uint16_t kmerf2 = (uint16_t )((kp2->kp)[1]);
	uint32_t kmers1 = (uint32_t )((kp1->kp)[2]);
	uint32_t kmers2 = (uint32_t )((kp2->kp)[2]);

	if (kmerf1 > kmerf2)
		return 1;
	else if (kmerf1 < kmerf2)
		return -1;
	else
	{
		if (kmers1 > kmers2)
			return 1;
		else if (kmers1 < kmers2)
			return -1;
		else
		{
#ifdef	UNPIPATH_OFF_K20
			/*
			uint64_t tmp_kp_pos1 = (kp1->kp)[3];
			uint64_t tmp_kp_pos2 = (kp2->kp)[3];
			tmp_kp_pos1 = ((tmp_kp_pos1 << 32) | ((kp1->kp)[0]));
			tmp_kp_pos2 = ((tmp_kp_pos2 << 32) | ((kp2->kp)[0]));

			if (tmp_kp_pos1 > tmp_kp_pos2)
			    return 1;
			else if (tmp_kp_pos1 < tmp_kp_pos2)
			    return -1;
			else
			{
			    printf("same pos\n");
				exit(2);
			    return 0;
			}
			*/

			if ((kp1->kp)[3] > (kp2->kp)[3])
				return 1;
			else if ((kp1->kp)[3] < (kp2->kp)[3])
				return -1;
			else
			{
				if ((kp1->kp)[0] > (kp2->kp)[0])
					return 1;
				else if ((kp1->kp)[0] < (kp2->kp)[0])
					return -1;
				else
				{
					printf("same pos\n");
					fflush(stdout);
					return 0;
				}
			}
#else
			if ((kp1->kp)[0] > (kp2->kp)[0])
				return 1;
			else if ((kp1->kp)[0] < (kp2->kp)[0])
				return -1;
			else
			{
				printf("same pos\n");
				fflush(stdout);
				return 0;
			}
#endif
		}
	}
}

int compare_pu (const void * a, const void * b)
{
	p_u* pu1 = (p_u *)a;
	p_u* pu2 = (p_u *)b;

	if(pu1->pos > pu2->pos)
		return 1;
	else if(pu1->pos < pu2->pos)
		return -1;
	else	return 0;
}

int compare_us (const void * a, const void * b)
{
	k_u* ku1 = (k_u *)a;
	k_u* ku2 = (k_u *)b;

	if(ku1->kmer > ku2->kmer)
		return 1;
	else if(ku1->kmer < ku2->kmer)
		return -1;
	else
	{
		return 0;
	}
}

uint32_t file_kmer_qsort()
{
#ifdef	LAST_DEBUG
	uint8_t k_i = 0;
	char char_s[32] = "";
	char char_l[64] = "";
	//uint8_t i_kmer = 0;

	uint32_t i = 0;
	uint8_t kmer_i = 0;

	uint32_t file_n = 0;
	file_n = (1 << (f << 1));
	char file_sta[ROUTE_LENGTH_MAX + KF_LENGTH_MAX];

	FILE* file_s = NULL;
	file_s = fopen(filename_sta,"r");
	if(file_s == NULL)	printf("Error of opening graph file\n");

	long int line_cnt = 0;
	uint32_t i_line = 0;

#ifdef	UNPIPATH_OFF_K20
	uint64_t line_tol = 0;
	uint64_t* file_line_cnt = (uint64_t* )calloc(file_n, 8);
	if(file_line_cnt == NULL)
	{
		printf("fail to allocate memory\n");
		exit(1);
	}
#else
	uint32_t line_tol = 0;
	uint32_t* file_line_cnt = (uint32_t* )calloc(file_n, 4);
	if(file_line_cnt == NULL)
	{
		printf("fail to allocate memory\n");
		exit(1);
	}
#endif

	while(!feof(file_s))
	{
		fgets(file_sta, 32, file_s);
		line_cnt = atol(file_sta);
#ifdef	UNPIPATH_OFF_K20
		if(i_line < file_n)
			file_line_cnt[i_line] = (uint64_t )line_cnt;
		if(i_line == file_n + 1)
			line_tol = (uint64_t )line_cnt;
#else
		if(i_line < file_n)
			file_line_cnt[i_line] = (uint32_t )line_cnt;
		if(i_line == file_n + 1)
			line_tol = (uint32_t )line_cnt;
#endif
		i_line++;
	}

	//begin qsorting
	char file_div[KF_LENGTH_MAX];
	char file_open[ROUTE_LENGTH_MAX + KF_LENGTH_MAX];
	FILE** file_p = (FILE** )malloc(file_n * sizeof(FILE* ));
	if(file_p == NULL)
	{
		printf("fail to allocate memory\n");
		exit(1);
	}

	for(i = 0; i < file_n; i++)
		file_p[i] = NULL;

	k_p* k_p_file = NULL;

	//load into hash graph declaration
	uint32_t k_p_i = 0;
	uint64_t kmer_l = 0;
	uint32_t kmer_f = 0;
	uint32_t kmer_s = 0;
	uint64_t kmer_l_p = 0Xffffffffffffffff;
	uint32_t kmer_f_p = 0Xffffffff;
	uint64_t hash_n = ((uint64_t )1 << (k << 1));

	printf("hash_n: %"PRId64"\n", hash_n);//268435456

	//alloc the memory used to store hash graph
#ifdef	UNPIPATH_OFF_K20
	uint64_t* pos_array = (uint64_t* )calloc(line_tol, 8);
	if(pos_array == NULL)
	{
		printf("fail to allocate memory pos_array\n");
		fflush(stdout);
		exit(1);
	}
	uint64_t* kmer_point = (uint64_t* )calloc(line_tol + 1, 8);
	if(kmer_point == NULL)
	{
		printf("fail to allocate memory kmer_point\n");
		fflush(stdout);
		exit(1);
	}
	uint64_t* hash_array = (uint64_t* )calloc(hash_n + 1, 8);
	if(hash_array == NULL)
	{
		printf("fail to allocate memory\n");
		fflush(stdout);
		exit(1);
	}
#else
	uint32_t* pos_array = (uint32_t* )calloc(line_tol, 4);
	if(pos_array == NULL)
	{
		printf("fail to allocate memory pos_array\n");
		exit(1);
	}
	uint32_t* kmer_point = (uint32_t* )calloc(line_tol + 1, 4);
	if(kmer_point == NULL)
	{
		printf("fail to allocate memory kmer_point\n");
		exit(1);
	}


	uint32_t* hash_array = (uint32_t* )calloc(hash_n + 1, 4);
	if(hash_array == NULL)
	{
		printf("fail to allocate memory\n");
		exit(1);
	}
#endif
	//actually smaller than this amount
	uint32_t* kmer_array = (uint32_t* )calloc(line_tol, 4);
	if(kmer_array == NULL)
	{
		printf("fail to allocate memory\n");
		fflush(stdout);
		exit(1);
	}

	uint8_t* edge_array = (uint8_t* )calloc(line_tol, 1);
	if(edge_array == NULL)
	{
		printf("fail to allocate memory edge_array\n");
		fflush(stdout);
		exit(1);
	}

	uint8_t* edge_flag = (uint8_t* )calloc(line_tol, 1);
	if(edge_flag == NULL)
	{
		printf("fail to allocate memory edge_flag\n");
		fflush(stdout);
		exit(1);
	}

	printf("total number of chars to allocate: %" PRIu64 "\n", (((uint64_t )line_tol) << 2));//36215262392
	fflush(stdout);

	uint64_t kmerf_cnt = 0;
	uint32_t hash_add = 0;

	uint8_t input = 0;
	uint8_t output = 0;

#ifdef	UNPIPATH_OFF_K20
	uint64_t tmp_pos64 = 0;
	uint64_t file_line_cnt_p = 0;
#else
	uint32_t file_line_cnt_p = 0;
#endif

	const uint32_t edge_bit = 1;
	uint32_t hash_first = 0;

	//end declaration
	for(i = 0; i < (1 << (f << 1)); i++)
	{
		bit32_kmer(i, f, kmer_i, file_div)
		strcpy(file_open, (const char* )filename_div);
		strcat(file_open, (const char* )file_div);
		strcat(file_open, suff);

		if(file_p[i] == NULL)
		{
			file_p[i] = fopen(file_open, "rb"); //all of div files are existing by default
			if(file_p[i] == NULL)
			{
				continue;
			}

			k_p_file = (k_p* )calloc(file_line_cnt[i], sizeof(k_p));
			fread(k_p_file, sizeof(k_p), file_line_cnt[i], file_p[i]);

			qsort(k_p_file, file_line_cnt[i], sizeof(k_p), compare);

			printf("has finished %u file sorting: %u\n", i, file_line_cnt[i]);

			//load into hash graph
			kmer_l_p = 0Xffffffffffffffff;
			kmer_f_p = 0Xffffffff;

			for(k_p_i = 0; k_p_i < file_line_cnt[i]; k_p_i++)
			{
				//add position
				//also can write to file directly
#ifdef	UNPIPATH_OFF_K20
				tmp_pos64 = (k_p_file[k_p_i].kp)[3];

				pos_array[k_p_i + file_line_cnt_p] = ((tmp_pos64 << 32) | ((k_p_file[k_p_i].kp)[0]));

#else
				pos_array[k_p_i + file_line_cnt_p] = (k_p_file[k_p_i].kp)[0];
#endif
				//get first and whole, also get input and output char
				bit32a_bit64((k_p_file[k_p_i].kp), 3, k_t - f, k_t - k, kmer_i, kmer_l, kmer_f, kmer_s, input, output)

#ifdef	DEBUG_NOT_FOUND
				//0 8258096 126 8258096
				printf("\nreadin: k_p_i %u %u %u %u %u\n", k_p_i, (uint16_t )((k_p_file[k_p_i].kp)[1]), ((k_p_file[k_p_i].kp)[2]), kmer_f, kmer_l);
#endif
				bit32_kmer(kmer_l, k_t - f, k_i, char_l)

				if(kmer_l != kmer_l_p)
				{
					bit32_kmer(kmer_s, k_t - k, k_i, char_s)
					bit32_kmer(kmer_l, k_t - f, k_i, char_l)

					//new kmer begins and add pos number to kmer
					kmer_array[kmerf_cnt] = kmer_s;
#ifdef	DEBUG_NOT_FOUND
					printf("kmer: %u\n\n", kmer_s);
					if(kmerf_cnt > 0)
						printf("pre edge: %u\n", edge_array[kmerf_cnt - 1]);
#endif
					kmer_point[kmerf_cnt] = k_p_i + file_line_cnt_p;

					if(kmer_f != kmer_f_p)
					{
						//new hash address and add number to hash
						hash_add = (i << ((k - f) << 1)) + kmer_f;
						hash_array[hash_add] = kmerf_cnt;
#ifdef	DEBUG_NOT_FOUND
						printf("hash: %u %u\n\n\n", hash_add,kmerf_cnt );
#endif
						if(kmerf_cnt == 0)	hash_first = hash_add;
					}
					++kmerf_cnt;


				}
#ifdef	DEBUG_NOT_FOUND
				printf("edge: %u %u %u\n", kmerf_cnt, input, output);
#endif
				//add kmer edge
				if(charToDna5[input] <= 3)
					edge_array[kmerf_cnt - 1] |= (edge_bit << (7 - charToDna5[input]));

				if(charToDna5[output] <= 3)
					edge_array[kmerf_cnt - 1] |= (edge_bit << (3 - charToDna5[output]));
				else	edge_flag[kmerf_cnt - 1] = 1;

				kmer_l_p = kmer_l;
				kmer_f_p = kmer_f;

				//
				//kmer_s_p = kmer_s;
			}

			printf("has finished %u file loading\n",i);
			fflush(stdout);

			//end loading

			//end output

			file_line_cnt_p += file_line_cnt[i];

			free(k_p_file);
		}
	}
	kmer_point[kmerf_cnt] = line_tol;
	hash_array[hash_n] = kmerf_cnt;

	for(i = 0; i < file_n; i++)
	{
		if(file_p[i] != NULL)
			fclose(file_p[i]);
	}
	free(file_p);


	//for debug from here///////////////
	/*
		FILE* fp_kmer_point_w = fopen("/home/ghz/viruses_index_all/tmp_kmer_point", "wb");
	    if(fp_kmer_point_w == NULL)
		{
			printf("cannot open the kmer point tmp file\n");
			//exit(1);
		}

		fwrite(kmer_point, 8, kmerf_cnt + 1, fp_kmer_point_w);//line_tol

		fclose(fp_kmer_point_w);

		fp_kmer_point_w = fopen("/home/ghz/viruses_index_all/tmp_kmer_pos", "wb");
	    if(fp_kmer_point_w == NULL)
		{
			printf("cannot open the kmer pos tmp file\n");
			//exit(1);
		}

		fwrite(pos_array, 8, file_line_cnt_p, fp_kmer_point_w);//line_tol

		fclose(fp_kmer_point_w);
	*/

	////////////////

	//traverse the graph and get supernode

	char ws[KT_LENGTH_MAX+2] = "";
	char wsn[KT_LENGTH_MAX+2] = "";
	char fs[KF_LENGTH_MAX+2] = "";
	char ss[KF_LENGTH_MAX+2] = "";

	uint8_t l_flag = 0;
	uint8_t n_flag = 0;
	uint8_t d_flag = 0;
	uint8_t fy_flag = 0;
	uint8_t ry_flag = 0;

	//unipath seq file
	fp_us = fopen(uniseq, "w");
	if(fp_us == NULL)
	{
		printf("cannot open the unipath seq file\n");
		fflush(stdout);
		exit(1);
	}

	//unipath offset binary file
	fp_usf_b = fopen(uniseqf_b, "wb");
	if(fp_usf_b == NULL)
	{
		printf("cannot open the unipath seq offset binary file\n");
		fflush(stdout);
		exit(1);
	}

	/*
	//unipath branch node file
	fp_ue = fopen(uniedge, "wb");
	if(fp_ue == NULL)
	{
		printf("cannot open the unipath branch node file\n");
		exit(1);
	}
	*/

	uint32_t ue_n = 0;
	uint32_t uni_offset_ori = 0;
	uint64_t kmer_ws = 0;
	uint64_t kmer_wsli = 0;
	uint64_t kmer_fs = 0;
	uint64_t kmer_ss = 0;
	uint64_t kmer_wsn = 0;
	uint32_t t_hash_i = 0;

	uint8_t node_i = 0;
	uint8_t edge_i = 0;
	uint8_t edge_out = 0;
	int64_t binary_r = 0;
	uint64_t one_bit64 = 1;
#ifdef	UNPIPATH_OFF_K20
	uint64_t hash_v_p = 0;
	uint64_t t_kmer_i = 0;
#else
	uint32_t hash_v_p = 0;
	uint32_t t_kmer_i = 0;
#endif
	uint32_t super_node_id = 0;
	uint8_t uni_edge = 0;

	//for statistic
	uint32_t fy_branch_cnt = 0;
	uint32_t ry_branch_cnt = 0;
	uint32_t x_branch_cnt = 0;
	uint32_t l_blunt_cnt = 0;
	uint32_t b_blunt_cnt = 0;
	uint32_t linear_cnt = 0;
	uint32_t e_blunt_cnt = 0;
	uint32_t es_blunt_cnt = 0;
	uint8_t li_cnt = 0;
	uint32_t kmer_intervali = 0;

	//fill the blank of hash array
	for(t_hash_i = hash_n; t_hash_i > 0; t_hash_i--)
	{
		if((hash_array[t_hash_i] == 0) && (t_hash_i > hash_first))
			hash_array[t_hash_i] = hash_v_p;
		else	hash_v_p = hash_array[t_hash_i];
	}

	printf("finish filling the hash blank\n");
	fflush(stdout);

	uint64_t usf_n = 1;

#ifdef UNPIPATH_OFF_K20
	uint64_t uni_max_l = 0;
	uint64_t uni_l = 0;
	uint64_t uni_tmp = 0;
	uint64_t unioff = 0;
	uint64_t uni_kmer_off = 0;
	fwrite(&unioff, 8, 1, fp_usf_b);
	uint64_t* uni_offset = (uint64_t* )calloc(kmerf_cnt, 8);
#else
	uint32_t uni_max_l = 0;
	uint32_t uni_l = 0;
	uint32_t uni_tmp = 0;
	uint32_t unioff = 0;
	uint32_t uni_kmer_off = 0;
	fwrite(&unioff, 4, 1, fp_usf_b);
	uint32_t* uni_offset = (uint32_t* )calloc(kmerf_cnt, 4);
#endif

	if(uni_offset == NULL)
	{
		printf("Failed to allocate the unipath offset array\n");
		fflush(stdout);
		exit(1);
	}

	printf("begin creating unipath: %"PRId64"\n", hash_n);
	fflush(stdout);

	uint32_t uni_node_cnt = 0;

	for(t_hash_i = 0; t_hash_i < hash_n; t_hash_i++)//268435456
	{
		for(t_kmer_i = hash_array[t_hash_i]; t_kmer_i < hash_array[t_hash_i + 1]; t_kmer_i++)
		{
			//get the whole kmer
			//kmer_array[t_kmer_i]: get last kmer
			kmer_ws = (((uint64_t )t_hash_i) << ((k_t - k) << 1)) + kmer_array[t_kmer_i];

			bit32_kmer(kmer_ws, k_t, kmer_i, ws)

			node_i = node_indentity(edge_array[t_kmer_i], edge_flag[t_kmer_i]);

			if(node_i == 1)	linear_cnt++;

			if(node_i != 1)
			{
				//statistic
				if(node_i == 2)	fy_branch_cnt++;
				if(node_i == 3)	ry_branch_cnt++;
				if(node_i == 4)	x_branch_cnt++;
				if(node_i == 5)	l_blunt_cnt++;
				if(node_i == 6)	b_blunt_cnt++;
				if(node_i == 7)	e_blunt_cnt++;
				if(node_i == 8)	es_blunt_cnt++;

				//unipath edge

				if((node_i != 2) && (node_i != 7))
				{
					uni_edge |= ((edge_array[t_kmer_i] >> 4) << 4);
				}

				//trace to record super node ID
				uni_l = 0;
				fy_flag = 0;
				ry_flag = 0;
				if(node_i == 2)	fy_flag = 1;
				if((node_i == 3) || (node_i == 5))
				{
					if(strlen(ws) != k_t)
					{
						printf("kmer length error\n");
						fflush(stdout);
						exit(1);
					}

					uni_offset[t_kmer_i] = uni_kmer_off;

					++uni_kmer_off;

					unioff += strlen(ws);
					uni_l += strlen(ws);

					ry_flag = 1;

					fprintf(fp_us, "\n%s",ws);
				}

				d_flag = 0;
				if((node_i == 4) || (node_i == 6))
				{
					if(strlen(ws) != k_t)
					{
						printf("kmer length error\n");
						fflush(stdout);
						exit(1);
					}

					uni_offset[t_kmer_i] = uni_kmer_off;

					uni_kmer_off += k_t;

					unioff += strlen(ws);
					uni_l += strlen(ws);

					d_flag = 1;

					++super_node_id;

#ifdef UNPIPATH_OFF_K20
					fwrite(&unioff, 8, 1, fp_usf_b);
#else
					fwrite(&unioff, 4, 1, fp_usf_b);
#endif

					usf_n++;

					if(uni_l > uni_max_l)	uni_max_l = uni_l;

					fprintf(fp_us, "\n%s",ws);

					uni_edge |= (edge_array[t_kmer_i] & 0Xf);
					//fwrite(&uni_edge, 1, 1, fp_ue);
					ue_n++;

					uni_edge = 0;
				}

				if(node_i == 8)
				{
					if(strlen(ws) != k_t)
					{
						printf("kmer length error\n");
						fflush(stdout);
						exit(1);
					}

					uni_offset[t_kmer_i] = uni_kmer_off;

					uni_kmer_off += k_t;

					unioff += strlen(ws);
					uni_l += strlen(ws);

					++super_node_id;
#ifdef UNPIPATH_OFF_K20
					fwrite(&unioff, 8, 1, fp_usf_b);
#else
					fwrite(&unioff, 4, 1, fp_usf_b);
#endif
					usf_n++;

					if(uni_l > uni_max_l)	uni_max_l = uni_l;

					fprintf(fp_us, "\n%s",ws);

					//fwrite(&uni_edge, 1, 1, fp_ue);
					ue_n++;

					uni_edge = 0;

					continue;
				}

				uni_tmp = uni_l;
				for(edge_i = 0; edge_i < 4; edge_i++) //A,C,G,T
				{
					edge_out = (uint8_t )((edge_array[t_kmer_i] >> (3 - edge_i)) & 0x1);
					if(edge_out == 0)
					{
						continue;
					}

					uni_l = uni_tmp;

					next_kmer(kmer_ws, kmer_wsn, k_t, edge_i)

					bit32_kmer(kmer_wsn, k_t, kmer_i, wsn)

					uint64_div(kmer_wsn, kmer_fs, kmer_ss, k_t - k)

					bit32_kmer(kmer_fs, k, kmer_i, fs)
					bit32_kmer(kmer_ss, k_t-k, kmer_i, ss)
#ifdef	UNPIPATH_OFF_K20
					binary_r = binsearch_offset_index64(kmer_ss, kmer_array, hash_array[kmer_fs + 1] - hash_array[kmer_fs], hash_array[kmer_fs]);
#else
					binary_r = binsearch_offset_index(kmer_ss, kmer_array, hash_array[kmer_fs + 1] - hash_array[kmer_fs], hash_array[kmer_fs]);
#endif

					if(binary_r == -1)
					{
						printf("kmer search error\n");
#ifdef	DEBUG_64BIT
						//504 2240 3090 3088
						printf("binsearch_offset_index1, not found %u %u %"PRId64" %"PRId64"\n", kmer_fs, kmer_ss, hash_array[kmer_fs + 1], hash_array[kmer_fs]);
						uint64_t tmp_j = 0;
						printf("\n");
						for(tmp_j = hash_array[kmer_fs]; tmp_j < hash_array[kmer_fs + 1]; tmp_j++)
							printf("%u\n", kmer_array[tmp_j]);
#endif
						fflush(stdout);
						exit(1);
					}

					kmer_wsli = kmer_wsn;

					//record seq
					node_i = node_indentity(edge_array[binary_r], edge_flag[binary_r]);
					n_flag = 0;	//2 4 6

					if(((d_flag == 1) || (fy_flag == 1)) && ((node_i == 1) || (node_i == 2) || (node_i == 7)))
					{
						if(strlen(wsn) != k_t)
						{
							printf("kmer length error\n");
#ifdef	DEBUG_64BIT
							printf("wrong kmer length %u %u\n", strlen(wsn), k_t);
#endif
							fflush(stdout);
							exit(1);
						}

						uni_offset[binary_r] = uni_kmer_off;

						++uni_kmer_off;

						unioff += strlen(wsn);
						uni_l += strlen(wsn);

						fprintf(fp_us, "\n%s",wsn);

						//unipath input edge
						uni_edge |= ((edge_array[binary_r] >> 4) << 4);

						n_flag = 1;
					}

					//s_flag = 0; //3 5
					if((ry_flag == 1) && ((node_i == 1) || (node_i == 2) || (node_i == 7)))
					{
						fprintf(fp_us, "%c", Dna5Tochar[edge_i]);

						uni_offset[binary_r] = uni_kmer_off;

						++uni_kmer_off;

						++unioff;
						++uni_l;

						//s_flag = 1;
					}

					l_flag = 0;
					li_cnt = 0;
					while(node_indentity(edge_array[binary_r], edge_flag[binary_r]) == 1)
					{
						++li_cnt;
						edge_out = edgeout_node(edge_array[binary_r]);
						if(edge_out == 4)
						{
							printf("edge identification error \n");
							fflush(stdout);
							exit(1);
						}

						++unioff;
						++uni_l;

						fprintf(fp_us, "%c", Dna5Tochar[edge_out]);

						l_flag = 1;

						next_kmer(kmer_wsli, kmer_wsn, k_t, edge_out)

						uint64_div(kmer_wsn, kmer_fs, kmer_ss, k_t - k)

						bit32_kmer(kmer_fs, k, kmer_i, fs)
						bit32_kmer(kmer_ss, k_t-k, kmer_i, ss)
#ifdef	UNPIPATH_OFF_K20
						binary_r = binsearch_offset_index64(kmer_ss, kmer_array, hash_array[kmer_fs + 1] - hash_array[kmer_fs], hash_array[kmer_fs]);
#else
						binary_r = binsearch_offset_index(kmer_ss, kmer_array, hash_array[kmer_fs + 1] - hash_array[kmer_fs], hash_array[kmer_fs]);
#endif
						if(binary_r == -1)
						{
							printf("kmer search error\n");
#ifdef	DEBUG_64BIT
							printf("binsearch_offset_index2, not found %u %u %"PRId64" %"PRId64"\n", kmer_fs, kmer_ss, hash_array[kmer_fs + 1], hash_array[kmer_fs]);
#endif
							fflush(stdout);
							exit(1);
						}
						uni_offset_ori = uni_offset[binary_r];
						uni_offset[binary_r] = uni_kmer_off;

						++uni_kmer_off;

						kmer_wsli = kmer_wsn;
					}

					node_i = node_indentity(edge_array[binary_r], edge_flag[binary_r]);

					if(((node_i == 3) || (node_i == 4) || (node_i == 8)) && (l_flag == 1))
					{
						fseek(fp_us, -1, SEEK_CUR);
						--unioff;
						--uni_l;

						uni_offset[binary_r] = uni_offset_ori;

						--uni_kmer_off;
					}

					if(n_flag || l_flag || ry_flag)
					{
						++super_node_id;

#ifdef UNPIPATH_OFF_K20
						fwrite(&unioff, 8, 1, fp_usf_b);
#else
						fwrite(&unioff, 4, 1, fp_usf_b);
#endif
						usf_n++;

						if(uni_l > uni_max_l)	uni_max_l = uni_l;

						//unipath edge
						if((node_i == 2) || (node_i == 7))
							uni_edge |= (edge_array[binary_r] & 0Xf);
						else
						{
							if(l_flag)	uni_edge |= (1 << (3 - edge_out));
							else	uni_edge |= (edge_array[t_kmer_i] & 0Xf);
						}

						//fwrite(&uni_edge, 1, 1, fp_ue);
						ue_n++;

						uni_edge = 0;

						uni_kmer_off += (k_t - 1);
					}
				}
			}

			uni_node_cnt++;
		}

	}

	if(edge_array)	free(edge_array);
	if(edge_flag)	free(edge_flag);

	fprintf(fp_us,"\n");

	printf("number of supernode: %u\n", super_node_id);
	printf("number of linear node: %u \n", linear_cnt);
	printf("number of forward Y branch node: %u \n", fy_branch_cnt);
	printf("number of reverse Y branch node: %u \n", ry_branch_cnt);
	printf("number of X branch node: %u \n", x_branch_cnt);
	printf("number of linear blunt node: %u \n", l_blunt_cnt);
	printf("number of branch blunt node: %u \n", b_blunt_cnt);
	printf("number of end blunt node: %u \n", e_blunt_cnt);
	printf("number of ends blunt node: %u \n", es_blunt_cnt);

	fflush(stdout);

#ifdef	UNPIPATH_OFF_K20
	printf("the max length of supernode: %"PRId64"\n", uni_max_l);
#else
	printf("the max length of supernode: %u \n", uni_max_l);
#endif

	//end traverse
	//fclose(fp_ue);
	fclose(fp_us);
	fclose(fp_usf_b);

	printf("begin wrtting hash and kmer array and offset on unipath into file %"PRId64"\n", kmerf_cnt);
	fflush(stdout);

	fp_hash = fopen(unihash_g,"wb");
	if(fp_hash == NULL)
	{
		printf("wrong file route of hash to write\n");
		fflush(stdout);
		exit(1);
	}

	fp_kmer = fopen(unikmer_g,"wb");
	if(fp_kmer == NULL)
	{
		printf("wrong file route of kmer to write\n");
		fflush(stdout);
		exit(1);
	}

	fp_off = fopen(unioff_g,"wb");
	if(fp_off == NULL)
	{
		printf("wrong file route of off_set to write\n");
		fflush(stdout);
		exit(1);
	}

	fwrite(kmer_array, 4, kmerf_cnt, fp_kmer);

#ifdef UNPIPATH_OFF_K20
	fwrite(hash_array, 8, hash_n+1, fp_hash);
	fwrite(uni_offset, 8, kmerf_cnt, fp_off);
#else
	fwrite(hash_array, 4, hash_n+1, fp_hash);
	fwrite(uni_offset, 4, kmerf_cnt, fp_off);
#endif
	if(uni_offset)	free(uni_offset);

	fclose(fp_hash);
	fclose(fp_kmer);
	fclose(fp_off);
#else
	//create the unipath position file
	printf("Begin creating unipath position and its point file\n");
	fflush(stdout);

	//for debug from here //////////////////////////////////
	uint64_t one_bit64 = 1;
	uint64_t kmer_fs = 0;
	uint64_t kmer_ss = 0;
	uint64_t kmerf_cnt = 0;
	uint64_t usf_n = 0;
	uint64_t ue_n = 0;
	int64_t binary_r = 0;

	char fs[KF_LENGTH_MAX+2] = "";
	char ss[KF_LENGTH_MAX+2] = "";

	//FILE* fp_hash = NULL;
	//FILE* fp_kmer = NULL;
	FILE* fp_kmer_point = NULL;

	uint64_t a_size = 0;
	int64_t result_hash_g = 0;
	int64_t result_kmer_g = 0;
	uint64_t hash_n = 268435456;
	uint64_t kmer_n = 0;

	hash_n = ((hash_n + 1) << 3);

	printf("Load unipath hash\n");

	fp_hash = fopen("/home/ghz/viruses_index_all/unipath_g.hash", "rb");
	if(fp_hash == NULL)
	{
		fputs ("File error opening the graph hash file\n",stderr);
		exit (1);
	}

	a_size = (hash_n >> 3);
	uint64_t* hash_array = (uint64_t* ) malloc (hash_n);
	if (hash_array == NULL)
	{
		fputs("Memory error hash_array",stderr);
		exit(2);
	}

	// copy the file into the buffer:
	result_hash_g = fread (hash_array, 8, a_size, fp_hash);

	if (result_hash_g != a_size)
	{
		fputs("Reading error",stderr);
		exit(3);
	}

	fclose(fp_hash);

	//read input graph kmer file
	printf("Load unipath kmer\n");

	fp_kmer = fopen("/home/ghz/viruses_index_all/unipath_g.kmer", "rb");
	if(fp_kmer == NULL)
	{
		fputs ("File error opening the graph hash file\n",stderr);
		exit (1);
	}

	fseek(fp_kmer, 0, SEEK_END);// non-portable
	kmer_n = ftell(fp_kmer);
	rewind(fp_kmer);

	a_size = (kmer_n >> 2);

	//a_size = kmer_num;

	uint32_t* kmer_array = (uint32_t* ) malloc (kmer_n);
	if (kmer_array == NULL)
	{
		fputs ("Memory error kmer_array", stderr);
		exit (2);
	}

	// copy the file into the buffer:
	result_kmer_g = fread (kmer_array, 4, a_size, fp_kmer);
	if (result_kmer_g != a_size)
	{
		fputs ("Reading error kmer_array", stderr);
		exit (3);
	}

	fclose(fp_kmer);


	printf("Load unipath kmer point\n");

	fp_kmer_point = fopen("/home/ghz/viruses_index_all/tmp_kmer_point", "rb");
	if(fp_kmer_point == NULL)
	{
		fputs ("File error opening the graph hash file\n",stderr);
		exit (1);
	}

	fseek(fp_kmer_point, 0, SEEK_END);// non-portable
	kmer_n = ftell(fp_kmer_point);
	rewind(fp_kmer_point);

	a_size = (kmer_n >> 3);

	//a_size = kmer_num;

	uint64_t* kmer_point = (uint64_t* ) malloc (kmer_n);
	if (kmer_point == NULL)
	{
		fputs ("Memory error kmer_point",stderr);
		exit (2);
	}

	// copy the file into the buffer:
	result_kmer_g = fread (kmer_point, 8, a_size, fp_kmer_point);
	if (result_kmer_g != a_size)
	{
		fputs ("Reading error kmer_point",stderr);
		exit (3);
	}

	fclose(fp_kmer_point);

	printf("Load unipath kmer pos\n");

	fp_kmer_point = fopen("/home/ghz/viruses_index_all/tmp_kmer_pos", "rb");
	if(fp_kmer_point == NULL)
	{
		fputs ("File error opening the graph hash file\n",stderr);
		exit (1);
	}

	fseek(fp_kmer_point, 0, SEEK_END);// non-portable
	kmer_n = ftell(fp_kmer_point);
	rewind(fp_kmer_point);

	a_size = (kmer_n >> 3);

	//a_size = kmer_num;

	uint64_t* pos_array = (uint64_t* ) malloc (kmer_n);
	if (pos_array == NULL)
	{
		fputs ("Memory error pos_array",stderr);
		exit (2);
	}

	// copy the file into the buffer:
	result_kmer_g = fread (pos_array, 8, a_size, fp_kmer_point);
	if (result_kmer_g != a_size)
	{
		fputs ("Reading error pos_array",stderr);
		exit (3);
	}

	fclose(fp_kmer_point);

	int uni_max_l = 600000;
	uint32_t kmer_i = 0;
	////////////////////////////////
#endif

	fp_us = fopen(uniseq, "r");
	if(fp_us == NULL)
	{
		printf("cannot read the unipath seq file\n");
		fflush(stdout);
		exit(1);
	}

	//unipath statistics
	/*
	fp_sta = fopen(unista, "wb");
	if(fp_sta == NULL)
	{
		printf("cannot open the unipath statistical file\n");
		fflush(stdout);
		exit(1);
	}
	*/
	//unipath position file
	fp_up = fopen(unipos, "wb");
	if(fp_up == NULL)
	{
		printf("cannot open the unipath position file\n");
		fflush(stdout);
		exit(1);
	}

	//unipath position point file
	fp_upp = fopen(uniposp, "wb");
	if(fp_upp == NULL)
	{
		printf("cannot open the unipath position point file\n");
		fflush(stdout);
		exit(1);
	}

	uint64_t up_n = 0;
	uint64_t upp_n = 0;

	uint8_t uni_de = 0;
	uint32_t pos_l_i = 0;

	char* uni_in_seq = (char* )malloc(uni_max_l + 2);
	if(uni_in_seq == NULL)
	{
		printf("fail to allocate memory for unipath seq input array\n");
		fflush(stdout);
		exit(1);
	}

	char uni_fkmer[KT_LENGTH_MAX];

	uint32_t uni_in_i = 0;
	uint32_t line_l = 0;
	uint32_t first_i = 0;
	uint32_t line_cntu = 0;
	uint64_t cnt_cl_cnt = 0;
	//the number of kmers in unipath seq hash kmer array
	uint64_t us_h_n = 0;

#ifdef	UNPIPATH_OFF_K20
	fwrite(&first_cnt_cl, 8, 1, fp_upp);
	uint64_t* first = NULL;
	uint64_t* first_pos = NULL;
	uint64_t first_cnt = 0;
	uint64_t r_cnt = 0;
	uint64_t uni_pos_i = 0;
	uint64_t pos_first = 0;
	uint64_t* r = NULL;
#else
	fwrite(&first_cnt_cl, 4, 1, fp_upp);
	uint32_t* first = NULL;
	uint32_t* first_pos = NULL;
	uint32_t first_cnt = 0;
	uint32_t r_cnt = 0;
	uint32_t uni_pos_i = 0;
	uint32_t pos_first = 0;
	uint32_t* r = NULL;
#endif
	int64_t second = 0;
	uint64_t uni_kmer = 0;

	upp_n++;

	//transfer to binary seq file
#ifdef UNI_SEQ64
	uint64_t uni_arr[UNI_SEQ_WRI_ARR];
#else
	uint8_t uni_arr[UNI_SEQ_WRI_ARR];
#endif
	uint32_t s_b_tmp = 0;
	uint64_t a_cnt = 0;
	uint32_t seq_i = 0;
	uint32_t w_offset = 0;


	fp_us_b = fopen(uniseq_b, "wb");
	if(fp_us_b == NULL)
	{
		printf("cannot open the unipath seq file\n");
		fflush(stdout);
		exit(1);
	}
#ifdef UNI_SEQ64
	memset(uni_arr, 0, (UNI_SEQ_WRI_ARR << 3));
#else
	memset(uni_arr, 0, (UNI_SEQ_WRI_ARR));
#endif
	//repetitive read the last line

	fgets(uni_in_seq, uni_max_l + 2, fp_us);
	while((!feof(fp_us)) && (fgets(uni_in_seq, uni_max_l + 2, fp_us) > 0))
	{
		line_cntu++;
		line_l = strlen(uni_in_seq) - 1;

		if(!line_l)	continue;

		if(line_l < k_t)
		{
			printf("Error of unipath seq file: %u %u %s\n", line_cntu, line_l, uni_in_seq);
			fflush(stdout);
			exit(1);
		}
		else
		{
			us_h_n += (line_l - k_t + 1);

			//transfer to binary seq file
			for(seq_i = 0; seq_i < line_l; seq_i++)
			{
#ifdef UNI_SEQ64
				binary64_array64_f(uni_arr, w_offset, s_b_tmp, charToDna5[(uint8_t )uni_in_seq[seq_i]], UNI_SEQ_WRI_ARR, a_cnt, fp_us_b)
#else
				binary64_array(uni_arr, w_offset, s_b_tmp, charToDna5[(uint8_t )uni_in_seq[seq_i]], UNI_SEQ_WRI_ARR, a_cnt, fp_us_b)
#endif
				w_offset++;
			}

			strncpy(uni_fkmer, uni_in_seq, k_t);

			kmer_bit64(uni_kmer, k_t, kmer_i, uni_fkmer)

			uint64_div(uni_kmer, kmer_fs, kmer_ss, k_t - k)

			//debug
			bit32_kmer(kmer_fs, k, kmer_i, fs)
			bit32_kmer(kmer_ss, k_t-k, kmer_i, ss)
#ifdef	UNPIPATH_OFF_K20
			binary_r = binsearch_offset_index64(kmer_ss, kmer_array, hash_array[kmer_fs + 1] - hash_array[kmer_fs], hash_array[kmer_fs]);
#else
			binary_r = binsearch_offset_index(kmer_ss, kmer_array, hash_array[kmer_fs + 1] - hash_array[kmer_fs], hash_array[kmer_fs]);
#endif
			if(binary_r == -1)
			{
				printf("kmer search error\n");
#ifdef	DEBUG_64BIT
				printf("binsearch_offset_index, not found %u %u %"PRId64" %"PRId64"\n", kmer_fs, kmer_ss, hash_array[kmer_fs + 1], hash_array[kmer_fs]);
#endif
				fflush(stdout);
				exit(1);
			}
			first_cnt = kmer_point[binary_r + 1] - kmer_point[binary_r];

#ifdef	UNPIPATH_OFF_K20
			first = (uint64_t* )calloc(first_cnt, 8);
			first_pos = (uint64_t* )calloc(first_cnt, 8);
#else
			first = (uint32_t* )calloc(first_cnt, 4);
			first_pos = (uint32_t* )calloc(first_cnt, 4);
#endif


			first_i = 0;
			for(uni_pos_i = kmer_point[binary_r]; uni_pos_i < kmer_point[binary_r + 1]; uni_pos_i++)
			{
				first[first_i] = pos_array[uni_pos_i];


				first_pos[first_i] = pos_array[uni_pos_i];

				first_i++;
			}

			uni_de = 0;

			for(uni_in_i = 1; uni_in_i < line_l - k_t + 1; uni_in_i++)
			{
				strncpy(uni_fkmer, uni_in_seq + uni_in_i, k_t);
				kmer_bit64(uni_kmer, k_t, kmer_i, uni_fkmer)

				uint64_div(uni_kmer, kmer_fs, kmer_ss, k_t - k)

#ifdef	UNPIPATH_OFF_K20
				second = binsearch_offset_index64(kmer_ss, kmer_array, hash_array[kmer_fs + 1] - hash_array[kmer_fs], hash_array[kmer_fs]);
#else
				second = binsearch_offset_index(kmer_ss, kmer_array, hash_array[kmer_fs + 1] - hash_array[kmer_fs], hash_array[kmer_fs]);
#endif

				if(second != -1)
				{
#ifdef	UNPIPATH_OFF_K20
					r_cnt = uni_pos_combine64(kmer_point, pos_array, first, first_cnt, second, &r);
#else
					r_cnt = uni_pos_combine(kmer_point, pos_array, first, first_cnt, second, &r);
#endif
				}
				else
				{
					printf("kmer search error\n");
#ifdef	DEBUG_64BIT
					printf("error -1 binsearch_offset_index, not found %u %u %"PRId64" %"PRId64"\n", kmer_fs, kmer_ss, hash_array[kmer_fs + 1], hash_array[kmer_fs]);
#endif
					fflush(stdout);
					exit(1);
				}

				if(!r_cnt)
				{
					if(!uni_de)
					{
						for(uni_pos_i = 0; uni_pos_i < first_cnt; uni_pos_i++)
						{
							for(pos_l_i = 0; pos_l_i < chr_posend_n; pos_l_i++)//chr_file_n_vcf
							{
								if(first_pos[uni_pos_i] < chr_end_n[pos_l_i + 1])
									break;
							}
						}

						uni_de = 1;
					}
				}

				first_cnt = r_cnt;
				first = r;

				if(!r_cnt)	break;

			}

			for(uni_pos_i = 0; uni_pos_i < first_cnt; uni_pos_i++)
			{
				pos_first = first[uni_pos_i] + 1 - uni_in_i;

				if(pos_first <= START_POS_REF)
				{
					printf("combine pos error\n");
					fflush(stdout);
					exit(1);
				}
#ifdef	UNPIPATH_OFF_K20
				fwrite(&pos_first, 8, 1, fp_up);
#else
				fwrite(&pos_first, 4, 1, fp_up);
#endif
				up_n++;
			}

			first_cnt_cl += first_cnt;
#ifdef	UNPIPATH_OFF_K20
			fwrite(&first_cnt_cl, 8, 1, fp_upp);
#else
			fwrite(&first_cnt_cl, 4, 1, fp_upp);
#endif
			upp_n++;

			cnt_cl_cnt++;

			free(first);

			free(first_pos);
		}
	}

#ifdef UNI_SEQ64
	fwrite(uni_arr, 8, UNI_SEQ_WRI_ARR, fp_us_b);
#else
	fwrite(uni_arr, 1, UNI_SEQ_WRI_ARR, fp_us_b);
#endif

	//free hash graph memory
	if(pos_array)	free(pos_array);
	if(kmer_array)	free(kmer_array);
	if(kmer_point)	free(kmer_point);
	if(hash_array)	free(hash_array);
	if(uni_in_seq)	free(uni_in_seq);
	//for index file close
	fclose(fp_us_b);
	fclose(fp_up);
	fclose(fp_upp);
	fclose(fp_us);

	printf("begin writing info about size of index files\n");
	fflush(stdout);

	fp_num = fopen(unisize, "wb");
	if(fp_num == NULL)	printf("error of creating index size file\n");

#ifdef UNI_SEQ64
	sta_num_write = ((a_cnt + 1) * (UNI_SEQ_WRI_ARR << 3));
#else
	sta_num_write = ((a_cnt + 1) * (UNI_SEQ_WRI_ARR));
#endif

	fwrite(&sta_num_write, 8, 1, fp_num);

#ifdef UNPIPATH_OFF_K20
	sta_num_write = (usf_n << 3);
#else
	sta_num_write = (usf_n << 2);
#endif

	fwrite(&sta_num_write, 8, 1, fp_num);

	sta_num_write = ue_n;
	fwrite(&sta_num_write, 8, 1, fp_num);

#ifdef	UNPIPATH_OFF_K20
	sta_num_write = (up_n << 3);
#else
	sta_num_write = (up_n << 2);
#endif
	fwrite(&sta_num_write, 8, 1, fp_num);

#ifdef	UNPIPATH_OFF_K20
	sta_num_write = (upp_n << 3);
#else
	sta_num_write = (upp_n << 2);
#endif
	fwrite(&sta_num_write, 8, 1, fp_num);

#ifdef	UNPIPATH_OFF_K20
	sta_num_write = ((hash_n + 1) << 3);
#else
	sta_num_write = ((hash_n + 1) << 2);
#endif
	fwrite(&sta_num_write, 8, 1, fp_num);

	sta_num_write = (kmerf_cnt << 2);
	fwrite(&sta_num_write, 8, 1, fp_num);

#ifdef UNPIPATH_OFF_K20
	sta_num_write = (kmerf_cnt << 3);
#else
	sta_num_write = (kmerf_cnt << 2);
#endif
	fwrite(&sta_num_write, 8, 1, fp_num);

	sta_num_write = 0;
	fwrite(&sta_num_write, 8, 1, fp_num);

	sta_num_write = (((ref_seq_n >> 5) + 1) << 3);
	fwrite(&sta_num_write, 8, 1, fp_num);

	if(tree_flag)	fwrite(&ref_seq_n, 8, 1, fp_num);

	fflush(fp_num);

	fclose(fp_num);


#ifdef	HANDLE_DIR
	//delete the div tmp folder
	printf("deleting div temporary index file folder\n");

	char rm_route[ROUTE_LENGTH_MAX];
	DIR* directory_pointer = NULL;

	if((directory_pointer = opendir(filename_div)) != NULL)
	{
		memset(rm_route, 0, ROUTE_LENGTH_MAX);
		strcpy(rm_route, sys_c_rm);
		strcat(rm_route, filename_div);
		system(rm_route);
	}

	memset(rm_route, 0, ROUTE_LENGTH_MAX);
	strcpy(rm_route, sys_rm);
	strcat(rm_route, uniseq);
	system(rm_route);
#endif

	printf("finish index building\n");

	return w_offset;
}

//create pos-unipath table
void build_pos_unipath()
{
	/*
	printf("begin creating position unipath file\n");
	#ifdef	UNPIPATH_OFF_K20
	printf("number of unipath positions: %"PRId64"\n",first_cnt_cl);
	#else
	printf("number of unipath positions: %u\n",first_cnt_cl);
	#endif
	//unipath position file
	fp_up = fopen(unipos, "rb");
	if(fp_up == NULL)	printf("cannot read the unipath position file\n");

	//unipath position point file
	fp_upp = fopen(uniposp, "rb");
	if(fp_upp == NULL)	printf("cannot read the unipath position point file\n");

	//position unipath file
	fp_pu = fopen(unipu, "wb");
	if(fp_pu == NULL)	printf("cannot open the  position unipath file\n");

	uint32_t p_cnt = 0;
	uint32_t p_cnt_i = 0;
	uint32_t p_cnt_pre = 0;
	uint32_t p_i = 0;

	p_u* p_u_file = NULL;

	p_u_file = (p_u* )calloc(first_cnt_cl, sizeof(p_u));

	while((!feof(fp_upp)) && (fread(&p_cnt, 4, 1, fp_upp) > 0))
	{
	    for(p_i = 0; p_i < p_cnt - p_cnt_pre; p_i++)
	    {
	        fread(&(p_u_file[p_u_c].pos), 4, 1, fp_up);
	        p_u_file[p_u_c].uniid = p_cnt_i;

	        p_u_c++;
	    }
	    p_cnt_pre = p_cnt;
	    p_cnt_i++;
	}

	qsort(p_u_file, p_u_c, sizeof(p_u), compare_pu);

	for(p_i = 0; p_i < p_u_c; p_i++)
	{
	    fwrite(&(p_u_file[p_i].pos),4, 1, fp_pu);
	    fwrite(&(p_u_file[p_i].uniid),4, 1, fp_pu);
	}

	free(p_u_file);
	fclose(fp_up);
	fclose(fp_upp);
	fclose(fp_pu);


	uint64_t p_u_c = 0;
	sta_num_write = (p_u_c << 3);
	fwrite(&sta_num_write, 8, 1, fp_num);
	*/
}

#ifdef	UNPIPATH_OFF_K20
uint64_t uni_pos_combine64(uint64_t kmer_point[], uint64_t pos_array[], uint64_t first[], uint64_t first_cnt, int64_t second, uint64_t** r)
{
	uint64_t pos_i = 0;
	uint64_t r_cnt = 0;
	uint64_t alloc_cnt = 0;
	uint64_t r_i = 0;
	int64_t b_r = 0;

	alloc_cnt = kmer_point[second + 1] - kmer_point[second];

	r_cnt = (first_cnt > alloc_cnt ? first_cnt : alloc_cnt);
	(*r) = (uint64_t* )calloc(r_cnt, 8);

	//can change to less time usage
	for(pos_i = 0; pos_i < first_cnt; pos_i++)
	{
		if((b_r = binsearch_offset_index_com64(first[pos_i] + 1, pos_array, alloc_cnt, kmer_point[second])) != -1)
		{
			(*r)[r_i] = first[pos_i] + 1;
			r_i++;
		}
	}

#ifdef	DEBUG_COMBINE
	if(r_i == 0)
	{
		printf("new\nposes:\n");
		for(pos_i = 0; pos_i < first_cnt; pos_i++)
			printf("%"PRId64"\n", first[pos_i]);
		printf("\n");

		printf("poses:\n");
		uint64_t tmp_se = 0;
		for(tmp_se = kmer_point[second]; tmp_se < kmer_point[second + 1]; tmp_se++)
			printf("%"PRId64"\n", pos_array[tmp_se]);
	}
#endif

	free(first);


	return r_i;
}
#else
uint32_t uni_pos_combine(uint32_t kmer_point[], uint32_t pos_array[], uint32_t first[], uint32_t first_cnt, uint32_t second, uint32_t** r)
{
	uint32_t pos_i = 0;
	uint32_t r_cnt = 0;
	uint32_t alloc_cnt = 0;
	uint32_t r_i = 0;
	int64_t b_r = 0;

	alloc_cnt = kmer_point[second + 1] - kmer_point[second];
	r_cnt = (first_cnt > alloc_cnt ? first_cnt : alloc_cnt);
	(*r) = (uint32_t* )calloc(r_cnt, 4);

	//can change to less time usage
	for(pos_i = 0; pos_i < first_cnt; pos_i++)
	{
		if((b_r = binsearch_offset_index(first[pos_i] + 1, pos_array, alloc_cnt, kmer_point[second])) != -1)
		{
			(*r)[r_i] = first[pos_i] + 1;
			r_i++;
		}
	}

	free(first);

	return r_i;
}
#endif
//below is some operations on graph

//get the number of 1 in one value
uint8_t one_number_uint(uint8_t n)
{
	n=(n & 0x55555555)+((n>>1) & 0x55555555);
	n=(n & 0x33333333)+((n>>2) & 0x33333333);
	n=(n & 0x0F0F0F0F)+((n>>4) & 0x0F0F0F0F);
	n=(n & 0x00FF00FF)+((n>>8) & 0x00FF00FF);
	n=(n & 0x0000FFFF)+((n>>16) & 0x0000FFFF);

	return n;
}

//identify branching node
uint8_t node_indentity(uint8_t edge, uint8_t edge_flag)
{
	uint8_t node_i = 0;

	if(edge_flag)	node_i =  2;
	else
	{
		uint8_t edge_out = 0;
		uint8_t edge_in = 0;
		edge_out = (edge & 0xf);
		edge_in = (edge >> 4);

		if((one_number_uint(edge_in) == 1) && (one_number_uint(edge_out) == 1))	node_i =  1; //linear
		if((one_number_uint(edge_in) == 1) && (one_number_uint(edge_out) > 1))	node_i =  2; //forward Y branch
		if((one_number_uint(edge_in) > 1) && (one_number_uint(edge_out) == 1))	node_i =  3; //reverse Y branch
		if((one_number_uint(edge_in) > 1) && (one_number_uint(edge_out) > 1))	node_i =  4; //X branch
		if((one_number_uint(edge_in) == 0) && (one_number_uint(edge_out) == 1))	node_i =  5; //linear blunt
		if((one_number_uint(edge_in) == 0) && (one_number_uint(edge_out) > 1))	node_i =  6; //branch blunt
		if((one_number_uint(edge_in) == 1) && (one_number_uint(edge_out) == 0))	node_i =  7; //end blunt
		if((one_number_uint(edge_in) > 1) && (one_number_uint(edge_out) == 0))	node_i =  8; //ends blunt
	}

	return node_i;
}

//get the unique char of out edge
uint8_t edgeout_node(uint8_t edge)
{
	uint8_t outedge = 4;
	uint8_t edge_out = 0;
	edge_out = (edge & 0xf);

	if(edge_out == 1)	outedge = 3;
	if(edge_out == 2)	outedge = 2;
	if(edge_out == 4)	outedge = 1;
	if(edge_out == 8)	outedge = 0;

	return outedge;
}

//binary search offset
int64_t binsearch_offset_index_com64(uint64_t x, uint64_t v[], uint64_t n, uint64_t offset)
{
	int64_t low, high, mid;

	low = 0;
	high = n - 1;

	while ( low <= high )
	{
		mid = (low + high) >> 1;
		if(x < v[mid + offset])
		{
			high = mid - 1;
		}
		else if(x > v[mid + offset])
		{
			low = mid + 1;
		}
		else  /*found match*/
		{
			return mid + offset;
		}
	}

	return -1;
}
int64_t binsearch_offset_index64(uint32_t x, uint32_t v[], uint64_t n, uint64_t offset)
{
	int64_t low, high, mid;

	low = 0;
	high = n - 1;

	while ( low <= high )
	{
		mid = (low + high) >> 1;
		if(x < v[mid + offset])
		{
			high = mid - 1;
		}
		else if(x > v[mid + offset])
		{
			low = mid + 1;
		}
		else  /*found match*/
		{
			return mid + offset;
		}
	}

	return -1;
}
int64_t binsearch_offset_index(uint32_t x, uint32_t v[], uint32_t n, uint32_t offset)
{
	int64_t low, high, mid;

	low = 0;
	high = n - 1;

	while ( low <= high )
	{
		mid = (low + high) >> 1;
		if(x < v[mid + offset])
		{
			high = mid - 1;
		}
		else if(x > v[mid + offset])
		{
			low = mid + 1;
		}
		else  /*found match*/
		{
			return mid + offset;
		}
	}

	return -1;
}


