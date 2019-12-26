/*
 * da_idx.hpp
 *
 *  Created on: 2018-10-26
 *      Author: Gaoyang Li
 */

#ifndef SRC_IDX_H_
#define SRC_IDX_H_
#include "bwt.h"
#include "lib/kvec.h"

typedef struct{
	char ref_name[128] ;
	uint64_t seq_l;
	uint64_t seq_offset;
}REF_INFO;

typedef struct {
	uint32_t ref_list;
	uint32_t length;
	//uint64_t offset;
}UNITIG;

typedef struct{
	uint32_t direction:1, ref_ID:31;
	uint32_t UNITIG_ID;
	uint32_t ref_offset;
}REF_UNI;

#define MIN_UNI_L 35//the min length for storing unitig when classify

typedef struct//UINTIG to REFRENCE mapping; one uintig may come from multiple reference positions
{
	uint64_t
	global_offset:40,
	ref_ID:23,
	direction:1;
}REF_POS;

kvec_T(REF_POS, REF_POS_V)

kvec_T(uint8_t, uint8_V)
kvec_T(uint64_t, uint64_V)
kvec_T(UNITIG, UNITIG_V)
kvec_T(REF_UNI, REF_UNI_V)
kvec_T(REF_INFO, REF_INFO_V)

#define L_PRE_IDX 13
#define PRE_IDX_MASK 0x3FFFFFF//26 bit 1

#define Q_MEM_MAX 2000
#define MAX_LV_WRONG 20
#define MAX_LV_R_LEN 20
#define Q_CONF_MAX 100
typedef struct
{
	int *Q_MEM;//[Q_MEM_MAX];
	int (*Q_LV)[MAX_LV_R_LEN];//[MAX_LV_WRONG][MAX_LV_R_LEN]
}MAPQ;

typedef struct
{
	uint8_t 	len_e_kmer;
	uint64_t 	e_kmer_hash_mask;
	uint64_t 	e_kmer_size;//size of ekmer table
	uint8_t* 	e_kmer0;
	uint8_t* 	e_kmer1;
	int 		single_base_max;//the max number of a single base in a kmer; this value are used to filter less complex kmers in pacbio reads
}E_KMER;

typedef struct
{
	//USED ONlY in building
	uint64_t 	n_kmer;
	REF_UNI_V 	r_u_v;//replace by r_p_v in classifying

	//USED in building and classifying
	bwt 		bt;
	E_KMER 		ek;
	UNITIG_V 	unitig_v;
	uint8_V		uni_str;
	REF_POS_V 	r_p_v;//replace r_u_v in classifying

	REF_INFO_V	r_i_v;
	uint8_V		ref_bin;

	MAPQ 		mapQ;
	int 		filter_min_length;
	int 		filter_min_score;
	int 		filter_min_score_LV3;
	//int			strain_mode;
}DA_IDX;

typedef struct
{
	uint8_t 	value[8];
	int8_t  	sp_pos;
	int8_t 		last_char;
}SP_KMER;

void set_ekmer_par(E_KMER *ek);
void load_idx(DA_IDX *idx, const char *dirPath);
#endif /* SRC_IDX_H_ */
